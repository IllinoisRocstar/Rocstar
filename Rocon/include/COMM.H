/// \file
/// \ingroup irad_group
/// \brief Communication utilities.
///
/// Encapsulates the entire communication substrate.  Nothing
/// outside these objects should know about MPI.
#ifndef _COMM_H_
#define _COMM_H_
#include "mpi.h"
#include <cstdlib>

#include "primitive_utilities.H"

namespace IRAD { 

  /// MPI Interface. 
  namespace Comm {
    /// Supported data types.
    enum DataTypes {DTDOUBLE,DTFLOAT,DTINT,DTUINT,DTSIZET,DTCHAR,DTUCHAR,DTBYTE,DTUBYTE};
    /// Operations for collectives.
    enum Ops {MAXOP, MINOP, SUMOP, PRODOP, MINLOCOP,MAXLOCOP};

    /// \brief provides communication for complex objects.
    ///
    /// The MobileObject provides buffering and an interface
    /// for complex objects to use communication routines.
    class MobileObject 
    {
    protected:
      void *_buf;
      bool  _mine;
      int   _bsize;
    public:
      MobileObject(): _buf(NULL), _mine(false), _bsize(0) {};
      void **GetBufPtr() { return(_buf ? &_buf : NULL); };
      const void *GetBuffer() const { return (_buf); };
      void *GetBuffer() {return (_buf); };
      int BufSize() { return(_bsize); };
      virtual int PrepareBuffer(size_t bsize);
      virtual int Pack(void **inbuf = NULL);
      virtual int UnPack(const void *outbuf = NULL);
      virtual void DestroyBuffer() 
      { 
	if(_buf && _mine) 
	  delete [] (char*)_buf; 
	_buf = NULL;
	_mine = false;
	_bsize = 0;
      };
      virtual ~MobileObject() 
      { 
	DestroyBuffer();
      };
    };

    /// \brief Main encapsulation of MPI.
    ///
    /// The CommunicatorObject provides an interface to MPI.
    /// This is for convenience, consistency, and to provide an
    /// "easy" access point for replacing the communication substrate
    /// with something that's not MPI.
    class CommunicatorObject
    {
    private:
      int _rank;
      MPI_Comm _comm;
      bool _master;
      bool _own_comm;
      bool _initd;
      int _error;
      int _rc;
      int _nproc;
      std::vector<MPI_Request> _send_requests;
      std::vector<MPI_Request> _recv_requests;
      std::vector<int> _send_tags;
      std::vector<int> _recv_tags;
      std::vector<MPI_Status> _status;
    public:
      MPI_Comm GetCommunicator(){ return (_comm); };
      MPI_Datatype IntegerTypeID;// = MPI_INTEGER;
      CommunicatorObject();
      CommunicatorObject(MPI_Comm &incomm);
      CommunicatorObject(int* narg,char*** args);
      MPI_Datatype ResolveDataType(const Comm::DataTypes &dt);
      MPI_Op ResolveOp(const Comm::Ops &op);
      int Split(int color,int key,CommunicatorObject &newcomm);
      int WaitRecv(int recvid);
      // wait on any and all messages
      int WaitAll();
      // clear up any persistent requests
      void ClearRequests();
      int NOpenRequests() { return(_send_requests.size() + _recv_requests.size()); };
      int Initialize(CommunicatorObject &incomm);
      int Initialize(int* narg = NULL,char*** args = NULL);
      int SetExit(int errin){return(_error = errin);};
      int SetErr(int errin){return(_error = errin);};
      void ClearErr(){_error = 0;};
      int Check(Comm::Ops op=Comm::MAXOP);
      int Rank();
      MPI_Comm World(){return(_comm);};
      int Finalize();
      int Barrier(){return(MPI_Barrier(_comm));};
      int Size();
      ~CommunicatorObject();
      int StartSend(unsigned int rid);
      int SendAll();
      int StartRecv(unsigned int rid);
      int RecvAll();
      int BroadCast(std::string &sval,int root_rank = 0);
      int BroadCast(MobileObject *mo,int root_rank = 0);
      int _ASend(void *buf,int sendsize,unsigned int remote_rank,int tag = 0);
      int _SetSend(void *buf,int sendsize,unsigned int remote_rank,int tag = 0);
      int _ARecv(void *buf,int recvsize,unsigned int remote_rank,int tag = 0);
      int _SetRecv(void *buf,int recvsize,unsigned int remote_rank,int tag = 0);
      int _AllGatherv(void *sendbuf,int mysendcnt,int datasize,void *recvbuf);
      int _Send(void *buf,int sendsize,unsigned int remote_rank,int tag = 0);
      int _Recv(void *buf,int recvsize,unsigned int remote_rank,int tag = 0);
      // Native support for Mobile Objects
      int _BroadCastMOV(std::vector<MobileObject *> &mos,int root_rank=0);
      int _GatherMO(MobileObject *sPtr,std::vector<MobileObject *> &rVec,int sndcnt,int root = 0);
      int _GatherMOV(std::vector<MobileObject *> &sVec,std::vector<MobileObject *> &rVec,
		     std::vector<int> &nsend_all,int root = 0);
      int _AllGatherMO(MobileObject *sPtr,std::vector<MobileObject *> &rVec,int sndcnt=1);
      int _AllGatherMOV(std::vector<MobileObject *> &sVec,std::vector<MobileObject *> &rVec,
			std::vector<int> &nsend_all);


      template<typename DataType>
      int ASend(std::vector<DataType> &sendbuf,unsigned int remote_rank,
		int tag = 0)
      {
	int sizeofdata = sizeof(DataType);
	int sendcnt = sendbuf.size();
	return (_ASend(&sendbuf[0],sizeofdata*sendcnt,remote_rank,tag));
      };

      template<typename DataType>
      int SetSend(std::vector<DataType> &sendbuf,unsigned int remote_rank,
		  int tag = 0)
      {
	int sizeofdata = sizeof(DataType);
	int sendcnt = sendbuf.size();
	return(_SetSend(&sendbuf[0],sendcnt*sizeofdata,remote_rank,tag));
      };

      template<typename DataType>
      int SetRecv(std::vector<DataType> &recvbuf,unsigned int remote_rank,
		  int tag = 0)
      {
	int sizeofdata = sizeof(DataType);
	int recvcnt = recvbuf.size();
	return(_SetRecv(&recvbuf[0],recvcnt*sizeofdata,remote_rank,tag));
      };

      template<typename DataType>
      int ARecv(std::vector<DataType> &recvbuf,unsigned int remote_rank,
		int tag=0)
      {
	int sizeofdata = sizeof(DataType);
	int recvcnt = recvbuf.size();
	return(_ARecv(&recvbuf[0],sizeofdata*recvcnt,remote_rank,tag));
      };

      template<typename DataType>
      int BroadCast(DataType &buf,int root_rank)
      {
	int sizeofdata = sizeof(DataType);
	return((_rc = MPI_Bcast(&buf,sizeofdata,MPI_CHAR,root_rank,_comm)));
      };

      template<typename DataType>
      int BroadCast(std::vector<DataType> &buf,int root_rank)
      {
	int sizeofdata = sizeof(DataType);
	int bufsize = buf.size();
	_rc = MPI_Bcast(&bufsize,1,MPI_INT,root_rank,_comm);
	if(_rank != root_rank)
	  buf.resize(bufsize);
	_rc = MPI_Bcast(&buf[0],buf.size()*sizeofdata,MPI_CHAR,root_rank,_comm);
	return(_rc);
      };
    
      // doesn't work, duh
      template<typename DataType>
      int Reduce(DataType &send,DataType &recv,
		 const Comm::DataTypes &dt,const Comm::Ops &op,int root)
      {
	_rc = MPI_Reduce(&send,&recv,1,ResolveDataType(dt),
			 ResolveOp(op),root,_comm);
	assert(_rc == 0);
	return(_rc);
      };

      // doesn't work, duh - can't reduce arbitrary DataType
      template<typename DataType>
      int Reduce(std::vector<DataType> &send,std::vector<DataType> &recv,
		 const Comm::DataTypes &dt,const Comm::Ops &op,int root)
      {
	int count = send.size();
	//    size_t datasize = sizeof(DataType);
	//    MPI_Datatype mpi_data_type = MPI_DOUBLE;
	//    if(datasize == sizeof(int))
	//      mpi_data_type = MPI_INTEGER;
	if(_rank == root)
	  recv.resize(count);
	_rc = MPI_Reduce(&send[0],&recv[0],count,ResolveDataType(dt),
			 ResolveOp(op),root,_comm);
	assert(_rc == 0);
	return(_rc);
      };

      // doesn't work, duh
      template<typename DataType>
      int AllReduce(std::vector<DataType> &send,std::vector<DataType> &recv,
		    const Comm::DataTypes &dt,const Comm::Ops &op)
      {
	int count = send.size();
	recv.resize(count);
	_rc = MPI_Allreduce(&send[0],&recv[0],count,ResolveDataType(dt),
			    ResolveOp(op),_comm);
	assert(_rc == 0);
	return(_rc);
      };


      // doesn't work, duh
      template<typename DataType>
      int AllReduce(DataType &send,DataType &recv,
		    const Comm::DataTypes &dt,const Comm::Ops &op)
      {
	_rc = MPI_Allreduce(&send,&recv,1,ResolveDataType(dt),
			    ResolveOp(op),_comm);
	assert(_rc == 0);
	return(_rc);
      };

      template<typename DataType>
      int AllGather(std::vector<DataType> &sendvec,std::vector<DataType> &recvvec,
		    int sndcnt=0,int recvcnt=0)
      {
	size_t datasize = sizeof(DataType);
	if(sndcnt == 0)
	  sndcnt = sendvec.size();
	if(recvcnt == 0)
	  recvcnt = sndcnt;
	std::cout << "Sendcount = " << sndcnt << " ReceiveCount = " << recvcnt << std::endl;
	_rc = MPI_Allgather((void *)(&(sendvec[0])),sndcnt*datasize,MPI_CHAR,
			    (void *)(&(recvvec[0])),recvcnt*datasize,MPI_CHAR,_comm);
	assert(_rc == 0);
	return(_rc);
      };


      //     template<typename T>
      //     int BroadCastMOVector(std::vector<T> &mov,int root_rank = 0)
      //     {
      //       std::vector<MobileObject *> moc;
      //       moc.resize(mov.size());
      //       std::vector<MobileObject *>::iterator moci = moc.begin();
      //       typename std::vector<T>::iterator oi = mov.begin();
      //       while(moci != moc.end())
      // 	*moci++ = dynamic_cast<MobileObject *>(&(*oi++));
      //       return(_BroadCastMOV(moc,root_rank));
      //     };

      template<typename DataType>
      int AllGather(DataType &sendval,std::vector<DataType> &recvvec)
      {
	size_t messagesize = sizeof(DataType);
	recvvec.resize(_nproc);
	_rc = MPI_Allgather((void *)&sendval,messagesize,MPI_CHAR,
			    (void *)&recvvec[0],messagesize,MPI_CHAR,_comm);
	assert(_rc == 0);
	return(_rc);
      };

      template<typename DataType>
      int AllGatherv(std::vector<DataType> &sendvec,std::vector<DataType> &recvvec,std::vector<int> &nsend_all)
      {
	int datasize = sizeof(DataType);
	int mysendcnt = sendvec.size();
	int totalcnt = 0;
	if(nsend_all.empty()){
	  nsend_all.resize(_nproc,0);
	  this->AllGather(mysendcnt,nsend_all);
	}
	this->AllReduce(mysendcnt,totalcnt,DTINT,SUMOP);
	recvvec.resize(totalcnt);
	return(_AllGatherv(&sendvec[0],mysendcnt,datasize,&recvvec[0]));
      };

      template<typename DataType>
      int AllGatherv(std::vector<DataType> &sendvec,std::vector<DataType> &recvvec)
      {
	int datasize = sizeof(DataType);
	int mysendcnt = sendvec.size();
	int totalcnt = 0;
	std::vector<int> nsend_all;
	nsend_all.resize(_nproc,0);
	this->AllGather(mysendcnt,nsend_all);
	this->AllReduce(mysendcnt,totalcnt,DTINT,SUMOP);
	recvvec.resize(totalcnt);
	return(_AllGatherv(&sendvec[0],mysendcnt,datasize,&recvvec[0]));
      };
    
      template<typename DataType>
      int Gather(DataType &sendval,std::vector<DataType> &recvvec,int root=0)
      {
	size_t messagesize = sizeof(DataType);
	if(this->Rank() == root){
	  recvvec.resize(this->Size());
	}
	_rc = MPI_Gather((void *)&sendval,messagesize,MPI_CHAR,
			 (void *)&recvvec[0],messagesize,MPI_CHAR,
			 root,_comm);
	assert(_rc == 0);
	return(_rc);
      };

      template<typename DataType>
      int Gather(std::vector<DataType> &sendvec,std::vector<DataType> &recvvec,
		 int sndcnt=0,int recvcnt=0,int root=0)
      {
	size_t datasize = sizeof(DataType);
	if(sndcnt == 0)
	  sndcnt = sendvec.size();
	if(recvcnt == 0)
	  recvcnt = sndcnt;
	if(_rank == root)
	  recvvec.resize(recvcnt*_nproc);
	_rc = MPI_Gather((void *)(&(sendvec[0])),sndcnt*datasize,MPI_CHAR,
			 (void *)(&(recvvec[0])),recvcnt*datasize,MPI_CHAR,
			 root,_comm);
	assert(_rc == 0);
	return(_rc);
      };

      template<typename DataType>
      int Gatherv(std::vector<DataType> &sendvec,std::vector<DataType> &recvvec,
		  std::vector<int> &nsend_all,int nsend = 0,int root=0)
      {
	int datasize = sizeof(DataType);
	if(nsend == 0)
	  nsend = sendvec.size();
	int nrecv = 0;
	if(nsend_all.empty()){
	  nsend_all.resize(_nproc);
	  nsend_all[_rank] = nsend;
	  this->Gather(nsend,nsend_all,root);
	}
	for(int i = 0;i < _nproc;i++)
	  nrecv += nsend_all[i];
	std::vector<int> allsizes(_nproc,0);
	std::vector<int> disps(_nproc,0);
        // MS : Debug: the following lines should be executed by all of the
        //      ranks, therefore if is not needed. Another problem is the
        //      datasize. MPI_CHAR datasize is 8, so 8-byte datatypes can be
        //      gathered from the client code by this method reliably.
	//if(_rank == root){
	  recvvec.resize(nrecv);
	  for(int i = 0; i < _nproc;i++){
	    allsizes[i] = nsend_all[i] * datasize;
	    if(i > 0)
	      disps[i] = disps[i-1]+allsizes[i-1];
	  }
	//}

	_rc = MPI_Gatherv((void *)(&(sendvec[0])),allsizes[_rank],MPI_CHAR,
			  (void *)(&(recvvec[0])),&allsizes[0],&disps[0],MPI_CHAR,
			  root,_comm);
        assert(_rc == 0);
	return(_rc);
      };

      // MS: Scatterv implementation for IRAD COMM object
      template<typename DataType>
      int Scatterv(std::vector<DataType> &sendvec, std::vector<int> &sendcounts,
                   std::vector<DataType> &recvvec, int root = 0)
      {
	int datasize = sizeof(DataType);
        if (sendcounts.size() == 0)
           sendcounts.resize(_nproc);
        // sending number of nrecv to each process
        this->BroadCast(sendcounts, root);
        // initializing recvec
        if (recvvec.size() == 0)
           recvvec.resize(sendcounts[_rank]);
        // calculating displs
	std::vector<int> allsizes(_nproc,0);
	std::vector<int> disps(_nproc,0);        
	for(int i = 0; i < _nproc;i++){
	    allsizes[i] = sendcounts[i] * datasize;
	    if(i > 0)
	      disps[i] = disps[i-1]+allsizes[i-1];
	}
        // scattering data
	_rc = MPI_Scatterv((void *)(&(sendvec[0])), &allsizes[0], &disps[0], MPI_CHAR,
			  (void *)(&(recvvec[0])), allsizes[_rank], MPI_CHAR,
			  root, _comm);
        assert(_rc == 0);
	return(_rc);
      };
      // MS: End

      // Native Mobile Object support below
      // ------------------------------------

      template<typename MOType>
      int BroadCastMobileObject(MOType &mo,int root_rank = 0)
      {
	return(this->BroadCast(dynamic_cast<MobileObject *>(&mo),root_rank));
      };

      template<typename MOType>
      int BroadCastMO(MOType &mo,int root_rank = 0)
      {
	return(this->BroadCast(dynamic_cast<MobileObject *>(&mo),root_rank));
      };
    
      template<typename MOType>
      int BroadCastMO(std::vector<MOType> &mov,int root_rank = 0)
      {
	std::vector<MobileObject *> moc;
	moc.resize(mov.size());
	std::vector<MobileObject *>::iterator moci = moc.begin();
	typename std::vector<MOType>::iterator oi = mov.begin();
	while(moci != moc.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	return(_BroadCastMOV(moc,root_rank));
      };

      template<typename MOType>
      int BroadCastMOVector(std::vector<MOType> &mov,int root_rank = 0)
      {
	std::vector<MobileObject *> moc;
	moc.resize(mov.size());
	std::vector<MobileObject *>::iterator moci = moc.begin();
	typename std::vector<MOType>::iterator oi = mov.begin();
	while(moci != moc.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	return(_BroadCastMOV(moc,root_rank));
      };

      // To do this properly, we need a GATHERV underneath the covers.  The
      // mobile object on each processor can be of a different size depending 
      // on the actual object's local implementation of the pack/unpack
      // routines. 
      template<typename MOType>
      int GatherMO(MOType &sendval,std::vector<MOType> &recvvec,int root = 0)
      {
	recvvec.resize(0);
	MobileObject *sendPtr = dynamic_cast<MobileObject *>(&sendval);
	std::vector<MobileObject *> recv_v;
	recv_v.resize(0);
	if(_rank == root){
	  recvvec.resize(_nproc);
	  recv_v.resize(_nproc,NULL);
	  std::vector<MobileObject *>::iterator moci = recv_v.begin();
	  typename std::vector<MOType>::iterator oi  = recvvec.begin();
	  while(moci != recv_v.end())
	    *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	}
	_rc = this->_GatherMO(sendPtr,recv_v,1,root);
	//      assert(_rc == 0);
	return(_rc);
      };

      template<typename MOType>
      int GatherMO(std::vector<MOType> &sendvec,std::vector<MOType> &recvvec,
		   std::vector<int> &nsend_all,int root = 0)
      {
	int nsend = sendvec.size();
	if(nsend_all.empty()){
	  nsend_all.resize(_nproc,0);
	  nsend_all[_rank] = nsend;
	  this->Gather(nsend,nsend_all,root);
	}
	int nrecv = 0;
	if(_rank == root){
	  std::vector<int>::iterator nsi = nsend_all.begin();
	  while(nsi != nsend_all.end())
	    nrecv += nsi++;
	  recvvec.resize(nrecv);
	}
	std::vector<MobileObject *> send_v(nsend,NULL);
	std::vector<MobileObject *> recv_v(nrecv,NULL);
	std::vector<MobileObject *>::iterator moci = recv_v.begin();
	typename std::vector<MOType>::iterator oi  = recvvec.begin();
	while(moci != recv_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	moci = send_v.begin();
	oi  = sendvec.begin();
	while(moci != send_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	_rc = this->_GatherMOV(send_v,recv_v,nsend_all,root);
	//      assert(_rc == 0);
	return(_rc);
      };

      template<typename MOType>
      int AllGatherMO(MOType &sendval,std::vector<MOType> &recvvec)
      {
	recvvec.resize(_nproc);
	MobileObject *sendPtr = dynamic_cast<MobileObject *>(&sendval);
	std::vector<MobileObject *> recv_v(_nproc,NULL);
	std::vector<MobileObject *>::iterator moci = recv_v.begin();
	typename std::vector<MOType>::iterator oi = recvvec.begin();
	while(moci != recv_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	_rc = this->_AllGatherMO(sendPtr,recv_v);
	assert(_rc == 0);
	return(_rc);
      };

      template<typename MOType>
      int AllGatherMO(std::vector<MOType> &sendvec,std::vector<MOType> &recvvec,
		      std::vector<int> &nsend_all,int root = 0)
      {
	int nsend = sendvec.size();
	if(nsend_all.empty()){
	  nsend_all.resize(_nproc,0);
	  nsend_all[_rank] = nsend;
	  this->Gather(nsend,nsend_all,root);
	}
	int nrecv = 0;
	if(_rank == root){
	  std::vector<int>::iterator nsi = nsend_all.begin();
	  while(nsi != nsend_all.end())
	    nrecv += nsi++;
	  recvvec.resize(nrecv);
	}
	std::vector<MobileObject *> send_v(nsend,NULL);
	std::vector<MobileObject *> recv_v(nrecv,NULL);
	std::vector<MobileObject *>::iterator moci = recv_v.begin();
	typename std::vector<MOType>::iterator oi  = recvvec.begin();
	while(moci != recv_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	moci = send_v.begin();
	oi  = sendvec.begin();
	while(moci != send_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	_rc = this->_GatherMOV(send_v,recv_v,nsend_all,root);
	//      assert(_rc == 0);
	return(_rc);
      };

      template<typename MOType>
      int GatherMOV(std::vector<MOType> &sendvec,std::vector<MOType> &recvvec,
		    std::vector<int> &nsend_all,int root = 0)
      {
	int nsend = sendvec.size();
	if(nsend_all.empty()){
	  nsend_all.resize(_nproc,0);
	  nsend_all[_rank] = nsend;
	  this->Gather(nsend,nsend_all,root);
	}
	int nrecv = 0;
	if(_rank == root){
	  std::vector<int>::iterator nsi = nsend_all.begin();
	  while(nsi != nsend_all.end())
	    nrecv += nsi++;
	  recvvec.resize(nrecv);
	}
	std::vector<MobileObject *> send_v(nsend,NULL);
	std::vector<MobileObject *> recv_v(nrecv,NULL);
	std::vector<MobileObject *>::iterator moci = recv_v.begin();
	typename std::vector<MOType>::iterator oi  = recvvec.begin();
	while(moci != recv_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	moci = send_v.begin();
	oi  = sendvec.begin();
	while(moci != send_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	_rc = this->_GatherMOV(send_v,recv_v,nsend_all,root);
	//      assert(_rc == 0);
	return(_rc);
      };

      template<typename MOType>
      int AllGatherMOV(std::vector<MOType> &sendvec,std::vector<MOType> &recvvec,
		       std::vector<int> &nsend_all,int root = 0)
      {
	int nsend = sendvec.size();
	if(nsend_all.empty()){
	  nsend_all.resize(_nproc,0);
	  nsend_all[_rank] = nsend;
	  this->Gather(nsend,nsend_all,root);
	}
	int nrecv = 0;
	if(_rank == root){
	  std::vector<int>::iterator nsi = nsend_all.begin();
	  while(nsi != nsend_all.end())
	    nrecv += nsi++;
	  recvvec.resize(nrecv);
	}
	std::vector<MobileObject *> send_v(nsend,NULL);
	std::vector<MobileObject *> recv_v(nrecv,NULL);
	std::vector<MobileObject *>::iterator moci = recv_v.begin();
	typename std::vector<MOType>::iterator oi  = recvvec.begin();
	while(moci != recv_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	moci = send_v.begin();
	oi  = sendvec.begin();
	while(moci != send_v.end())
	  *moci++ = dynamic_cast<MobileObject *>(&(*oi++));
	_rc = this->_AllGatherMOV(send_v,recv_v,nsend_all);
	//      assert(_rc == 0);
	return(_rc);
      };

      template < typename T >
      int Union(std::vector< T > &input_data,std::vector< T > &output_data)
      {
	std::vector< T > all_items;
	int err = AllGatherv< T >(input_data,all_items);
	std::sort(all_items.begin(),all_items.end());
	typename std::vector< T >::iterator ui = std::unique(all_items.begin(),all_items.end());
	output_data.resize(ui - all_items.begin());
	std::copy(all_items.begin(),ui,output_data.begin());
	return err;
      };

    

    };
  
  
    /// 
    /// Utility class for creating derived objects that are parallel.
    ///
    /// Inheriting from this class provides parallel capabilities.
    ///
    class ParallelObject {
    protected:
      CommunicatorObject _communicator;
    public:
      ParallelObject(){};
      ParallelObject(CommunicatorObject &incomm){
	_communicator.Initialize(incomm);
      };
      virtual ~ParallelObject(){};
      virtual CommunicatorObject &Communicator(){return(_communicator);};
    };
  };
};
#endif
