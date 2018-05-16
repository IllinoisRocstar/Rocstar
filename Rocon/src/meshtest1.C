#include <ctime>

#include "Profiler.H"
#include "Global.H"
#include "FEM.H"
#ifdef PMETIS
#include "parmetis.h"
#else
#include "metis.h"
#endif

///
/// \brief ComLineObject for testing app
///
class MeshTestComLine : public ComLineObject
{
public:
  MeshTestComLine(const char *args[])
    : ComLineObject(args)
  {};
  void Initialize(){
    AddOption('h',"help");
    AddOption('p',"partition",2,"number");
    AddOption('v',"verb",1,"level");
    AddOption('c',"clone",2,"number");
    AddOption('a',"assembly");
    AddOption('m',"mesh");
    AddOption('d',"debug");
    AddOption('k',"checking");
    AddOption('r',"reorient");
    AddOption('g',"generate");
    AddOption('s',"sparse");
    AddOption('t',"metis");
    AddOption('n',"renumber");
    AddArgument("input",1);
    AddHelp("metis","Metis testing stub.");
    AddHelp("sparse","Write out the sparse matrix for visualization.");
    AddHelp("clone","Generate <number> of partitions identical to the input mesh.");
    AddHelp("generate","Generate a uniform mesh with N nodes and quit.");
    AddHelp("checking","Paranoid and insanely verbose dumping of all important arrays to Log.");
    AddHelp("partition","Performs Metis partitioning of input mesh into <number> partitions.");
    AddHelp("assembly","Performs assembly test.");
    AddHelp("mesh","Performs mesh tests.");
    AddHelp("help","Prints this long version of help.");
    AddHelp("verb","Makes the test more verbose. Default level is 1.");
    AddHelp("config","Specifies the name of the configuration file.");
    AddHelp("out","Specifies the name of the output file.");
    AddHelp("renumber","Uses ParMETIS to do optimal graph reordering.");
    AddArgHelp("input","Mode dependent arguments");
    std::ostringstream Ostr;
    //    Ostr << "Use fixed problem size in scalability analysis.  Only makes"
    //	 << "\n\t\tsense when scalability mode is enabled.";
    //    Ostr.str("");
    Ostr << "Test tool for exercising the mesh library.";
    _description.assign(Ostr.str());
  };
};


void GetStatistics(std::ostream &Ostr,Mesh::Connectivity &con)
{
  Ostr << "Size = " << con.size() << std::endl
       << "Nelem = " << con.Nelem() << std::endl;
  double mean = 0.0;
  Mesh::IndexType max = 0;
  Mesh::IndexType min = 100000000;
  Mesh::IndexType minvalue = 10000000;
  Mesh::IndexType maxvalue = 0;
  Mesh::IndexType nvalues = 0;
  double meanvalue = 0.0;
  for(Mesh::IndexType iii = 0;iii < con.Nelem();iii++){
    mean += con.Esize(iii+1);
    if(max < con.Esize(iii+1))
      max = con.Esize(iii+1);
    if(min > con.Esize(iii+1))
      min = con.Esize(iii+1);
    for(Mesh::IndexType jjj = 0;jjj < con.Esize(iii+1);jjj++){
      if(minvalue > con[iii][jjj])
	minvalue = con[iii][jjj];
      if(maxvalue < con[iii][jjj])
	maxvalue = con[iii][jjj];
      meanvalue += con[iii][jjj];
      nvalues++;
      
    }
  }
  Ostr << "Esizes(min/max/mean): ("
       << min << "/" << max << "/"
       << mean/con.Nelem() << ")"
       << std::endl
       << "Flattened Size = " << nvalues
       << std::endl
       << "Values(min/max/mean): (" 
       << minvalue << "/" 
       << maxvalue << "/"
       << meanvalue/nvalues << ")"
       << std::endl;
}

// *** Warning ***
// This is called *after* the local dof id's have been mapped to the new dof id's and the 
// info.doffset has been reset to an appropriate value.  The NodalDofs array should be the new 
// global numbers with the doffset applied.   We need local dof id's (i.e. 1 - ndof_local) in 
// order to do the stiffness matrix properly.  So we need to determine the new "rank ordering" 
// (i.e. determined the doffset).
void GlobalDofReMapExchange(Mesh::Connectivity &ec,
			    Mesh::Connectivity &NodalDofs,
			    Mesh::Connectivity &ElementDofs,
			    Mesh::Connectivity &RemoteNodalDofs,
			    std::vector<Mesh::Border> &borders,
			    Mesh::PartInfo &info,
			    std::vector<Mesh::IndexType> &order,
			    Comm::CommunicatorObject &comm,
			    std::ostream *Out,bool verbose)
{
  // First, build and post the receives for the new remote nodal dof id's
  std::vector<std::vector<Mesh::IndexType> > RcvBuf;
  unsigned int nborders = borders.size();
  RcvBuf.resize(nborders);
  for(unsigned int i = 0;i < nborders;i++){
    RcvBuf[i].resize(borders[i].remote_dofcount);
    // Go ahead and post the receives for the new dof id's
    comm.ARecv<Mesh::IndexType>(RcvBuf[i],borders[i].rpart-1); // async recv
    if(Out && verbose)
      *Out << "Posting receive for " << RcvBuf[i].size() << " dofs from " 
	   << borders[i].rpart << std::endl;
  }
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs receive posted in GlobalDofReMapExchange." << std::endl;

  // Now form and send the new global dofs for nodes I own
  // *** Warning ***
  // This routine assumes that the dofs passed in the NodalDofs array need
  // to be offset by info.doffset.  If this is not the case, make sure that 
  // info.doffset is set to 0 upon entry.
  std::vector<std::vector<Mesh::IndexType> > SndBuf;
  SndBuf.resize(nborders);
  for(unsigned int i = 0;i < nborders;i++){
    SndBuf[i].resize(0);
    std::vector<Mesh::IndexType>::iterator rni = borders[i].nrecv.begin();
    while(rni != borders[i].nrecv.end()){
      std::vector<Mesh::IndexType>::iterator ldi = NodalDofs[*rni-1].begin();
      while(ldi != NodalDofs[*rni-1].end()){
	SndBuf[i].push_back(order[*ldi-1]+1); // + info.doffset?  nah
	ldi++;
      }
      rni++;
    }
    if(Out && verbose){
      *Out << "Posting send of " << SndBuf[i].size() << " nodal dofs to " 
	   << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,SndBuf[i]," ");
      *Out << std::endl;
    }
    //    borders[i].recvsize = SndBuf[i].size();
    comm.ASend<Mesh::IndexType>(SndBuf[i],borders[i].rpart-1); // asynchronous send
  }
  comm.WaitAll();
  comm.ClearRequests();
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs send posted 2." << std::endl;

  // Unpack the received global dofs into the RemoteNodalDofs array
  for(unsigned int i = 0;i < nborders;i++){
    std::vector<Mesh::IndexType>::iterator rdi = RcvBuf[i].begin();
    std::vector<Mesh::IndexType>::iterator sni = borders[i].nsend.begin();
    if(Out && verbose){
      *Out << "Unpacking recv'd dofs:" << std::endl;
      DumpContents(*Out,RcvBuf[i]," ");
      *Out << std::endl;
    }
    while(sni != borders[i].nsend.end()){
      Mesh::IndexType rnindex = *sni - info.nlocal - 1;
      std::vector<Mesh::IndexType>::iterator ndi = RemoteNodalDofs[rnindex].begin();
      while(ndi != RemoteNodalDofs[rnindex].end())
	*ndi++ = *rdi++;
      sni++;
    }
  }
  comm.Barrier();
  //
  // All Nodal Dofs have been updated.  Now need to update the communication buffers

  // Now calculate the indices of the stiffness entries I have for the remote 
  // nodal dofs. 

  Mesh::Connectivity dc;
  Mesh::IndexType number_of_nodes = NodalDofs.size();
  if(Out && verbose)
    *Out << "Forming dual connectivity for " << number_of_nodes << " nodes." << std::endl;
  ec.Inverse(dc,number_of_nodes);
  if(Out && verbose)
    *Out << "dual done." << std::endl;
  //  std::vector<std::vector<Mesh::IndexType> > SndBuf(nborders);
  std::vector<Mesh::IndexType> ssizes(nborders,0);
  for(unsigned int i = 0;i < nborders;i++){
    //    borders[i].data.SendAp.resize(0);
    borders[i].data.SendAi.resize(0);
    //    borders[i].data.SendAp.push_back(0);
    std::vector<Mesh::IndexType>::iterator sni = borders[i].nsend.begin();
    while(sni != borders[i].nsend.end()){
      // The list is the same for each border *node*
      std::list<Mesh::IndexType> doflist;
      // Count up the dof contributions from every local
      // element that touches the remote node
      std::vector<Mesh::IndexType>::iterator ei = dc[*sni-1].begin();
      while(ei != dc[*sni-1].end()){ // for every element touching the remote node
	std::vector<Mesh::IndexType>::iterator eni = ec[*ei-1].begin();
	while(eni != ec[*ei-1].end()){ // loop thru the element's nodes
	  std::vector<Mesh::IndexType>::iterator ndi;
	  std::vector<Mesh::IndexType>::iterator ndend;
	  if(*eni > info.nlocal){  // then it's a remote node
	    ndi   = RemoteNodalDofs[*eni - info.nlocal - 1].begin();
	    ndend = RemoteNodalDofs[*eni - info.nlocal -1].end();
	    while(ndi != ndend)
	      doflist.push_back(*ndi++);
	  }
	  else{
	    //	    offset = info.doffset;
	    ndi   = NodalDofs[*eni-1].begin();
	    ndend = NodalDofs[*eni-1].end();
	    while(ndi != ndend){
	      doflist.push_back(order[*ndi-1]+1);
	      ndi++;
	    }
	  }
	  eni++;
	}
	eni = ElementDofs[*ei-1].begin();
	while(eni != ElementDofs[*ei-1].end()){
	  doflist.push_back(order[*eni-1]+1);
	  eni++;
	}
	ei++;					
      }
      doflist.sort();
      doflist.unique();
      std::list<Mesh::IndexType>::iterator dli = doflist.begin();
      while(dli != doflist.end())
	borders[i].data.SendAi.push_back(*dli++);
      sni++;
    }
    if(Out && verbose){
      *Out << "SendAp for color " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.SendAp," ");
      *Out << std::endl
	   << "SendAi for color " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.SendAi," ");
      *Out << std::endl;
    }
  }
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs have remapped SendAi." << std::endl;
  // By now, the remote processor has sent us an array of sizes for the nz entries
  // relevant to each local *node*.  To get the total size of the stiffness block
  // that will be sent, we need to multiply by the number of actual dofs living 
  // on the local node.  Right now, we are just interested in the dof id's that 
  // we will be receiving/node.
  for(unsigned int i = 0;i < nborders;i++){
    //    borders[i].data.RecvAp.resize(RcvBuf[i].size()+1);
    //    std::vector<Mesh::IndexType>::iterator rbi = RcvBuf[i].begin();
    //    std::vector<Mesh::IndexType>::iterator rai = borders[i].data.RecvAp.begin();
    //    *rai = 0;
    //    while(rbi != RcvBuf[i].end()){
    //      Mesh::IndexType pv = *rai++;
    //      *rai = pv + *rbi++;
    //    }
    if(Out && verbose){
      *Out << "RecvAP from partition " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.RecvAp," ");
      *Out << std::endl;
    }
    // We can go ahead and post receives for the Ai
    std::vector<Mesh::IndexType>::reverse_iterator rsi = borders[i].data.RecvAp.rbegin();
    unsigned int totalsize = *rsi;
    borders[i].data.RecvAi.resize(totalsize);
    if(Out && verbose)
      *Out << "Posting RecvAi (" << borders[i].data.RecvAi.size()
	   << ")  from partition" << borders[i].rpart << std::endl;
    comm.ARecv<Mesh::IndexType>(borders[i].data.RecvAi,borders[i].rpart-1); // async recv
    // Note how we've already conveniently filled the SendAi, send it out now
    comm.ASend<Mesh::IndexType>(borders[i].data.SendAi,borders[i].rpart-1); // async recv
    if(Out && verbose){
      *Out << "Posting SendAi (" << borders[i].data.SendAi.size() 
	   << ") to partition" << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.SendAi," ");
      *Out << std::endl;
    }
  }
  comm.WaitAll();
  comm.ClearRequests();
  comm.Barrier();
  if(Out && verbose){
    for(unsigned int i = 0;i < nborders;i++){
      *Out << "Recieved RecvAi (" << borders[i].data.RecvAi.size() 
	   << ") from partition" << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.RecvAi," ");
      *Out << std::endl;
    }
  } 
}

void TestMPI(std::vector<Mesh::Border> &borders,Mesh::Connectivity &NodalDofs,
	     Mesh::Connectivity &RemoteNodalDofs,Mesh::PartInfo &info,
	     Comm::CommunicatorObject &comm,std::ostream *Out,bool verbose)
{
  unsigned int okcount = 0;
  unsigned int bcount = 0;
  for(unsigned int ntrial = 0;ntrial < 100;ntrial++){
    comm.Barrier();
    if(Out && verbose)
      *Out << "Inside MPI Test(" << ntrial << ")" << std::endl;
    // First just exchange the sizes (in number of dofs) for each node
    //
    // Build the receive buffers and post the receives
    //    Mesh::IndexType doffset = info.doffset;
    std::vector<std::vector<Mesh::IndexType> > RcvBuf;
    unsigned int nborders = borders.size();
    RcvBuf.resize(nborders);
    for(unsigned int i = 0;i < nborders;i++){
      // For each remotely owned node, we will receive one integer indicating the
      // number of dofs on that remotely owned node
      RcvBuf[i].resize(borders[i].nsend.size());
      if(Out && verbose)
	*Out << "Posting receive for " << borders[i].nsend.size() << " nodes from " 
	     << borders[i].rpart << std::endl << std::flush;
      comm.ARecv<Mesh::IndexType>(RcvBuf[i],borders[i].rpart-1); // async recv
    }
    comm.Barrier();
    if(Out && verbose)
      *Out << "All procs receives posted." << std::endl;
    // Build the send buffer for sending the number of dofs on each 
    // node on the partition boundary that I own
    std::vector<std::vector<Mesh::IndexType> > SndBuf;
    SndBuf.resize(nborders);
    for(unsigned int i = 0;i < nborders;i++){
      SndBuf[i].resize(borders[i].nrecv.size());
      unsigned int j = 0;
      std::vector<Mesh::IndexType>::iterator rni = borders[i].nrecv.begin();
      while(rni != borders[i].nrecv.end()){
	SndBuf[i][j] = NodalDofs[*rni-1].size();
	rni++;
	j++;
      }
      if(Out && verbose){
	*Out << "Posting send of " << borders[i].nrecv.size() << " node's dof sizes to " 
	     << borders[i].rpart << std::endl;
	*Out << "Sending: ";
	DumpContents(*Out,SndBuf[i]," ");
	*Out << std::endl << std::flush;
      }
      comm.ASend<Mesh::IndexType>(SndBuf[i],borders[i].rpart-1); // asynchronous send
    }
    comm.WaitAll();
    comm.ClearRequests();
    comm.Barrier();
    for(unsigned int i = 0;i < nborders;i++){
      *Out << "Received(" << i+1 << "): ";
      DumpContents(*Out,RcvBuf[i]," ");
      *Out << std::endl << std::flush;
      if(RcvBuf[i][0] == 1)
	okcount++;
      else
	bcount++;
    }
  }
  *Out << "Bad Receives: " << bcount << std::endl << std::flush;
  comm.Barrier();
}

void GlobalDofExchange(std::vector<Mesh::Border> &borders,Mesh::Connectivity &NodalDofs,
		       Mesh::Connectivity &RemoteNodalDofs,Mesh::PartInfo &info,
		       Comm::CommunicatorObject &comm,std::ostream *Out,bool verbose)
{
  if(Out && verbose)
    *Out << "Inside Global Dof Exchange" << std::endl;
  // First just exchange the sizes (in number of dofs) for each node
  //
  // Build the receive buffers and post the receives
  Mesh::IndexType doffset = info.doffset;
  std::vector<std::vector<Mesh::IndexType> > RcvBuf;
  unsigned int nborders = borders.size();
  RcvBuf.resize(nborders);
  for(unsigned int i = 0;i < nborders;i++){
    // For each remotely owned node, we will receive one integer indicating the
    // number of dofs on that remotely owned node
    RcvBuf[i].resize(borders[i].nsend.size());
    if(Out && verbose)
      *Out << "Posting receive for " << borders[i].nsend.size() << " nodes from " 
	   << borders[i].rpart << std::endl;
    comm.ARecv<Mesh::IndexType>(RcvBuf[i],borders[i].rpart-1); // async recv
  }
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs receives posted." << std::endl;
  // Build the send buffer for sending the number of dofs on each 
  // node on the partition boundary that I own
  std::vector<std::vector<Mesh::IndexType> > SndBuf;
  SndBuf.resize(nborders);
  for(unsigned int i = 0;i < nborders;i++){
    SndBuf[i].resize(borders[i].nrecv.size());
    unsigned int j = 0;
    std::vector<Mesh::IndexType>::iterator rni = borders[i].nrecv.begin();
    while(rni != borders[i].nrecv.end()){
      SndBuf[i][j] = NodalDofs[*rni-1].size();
      rni++;
      j++;
    }
    if(Out && verbose)
      *Out << "Posting send of " << borders[i].nrecv.size() << " node's dof sizes to " 
	   << borders[i].rpart << std::endl;
    comm.ASend<Mesh::IndexType>(SndBuf[i],borders[i].rpart-1); // asynchronous send
  }
  comm.WaitAll();
  comm.ClearRequests();
  comm.Barrier();
  
  if(Out && verbose)
    *Out << "All procs sends posted." << std::endl;
  // Unpack the receive buffer and appropriately size the NodalDof array for 
  // each remote node.
  std::vector<unsigned int> totalsizes(nborders,0);
  for(unsigned int i = 0;i < nborders;i++){
    std::vector<Mesh::IndexType>::iterator sni = borders[i].nsend.begin();
    unsigned int j = 0;
    if(Out && verbose){
      *Out << "Received node sizes (in dofs) from " << borders[i].rpart
	   << ":" << std::endl;
      DumpContents(*Out,RcvBuf[i]);
      *Out << std::endl;
    }
    while(sni != borders[i].nsend.end()){
      totalsizes[i] += RcvBuf[i][j];
  
      //
      // For now, we make sure that the number of local dofs 
      // on a remote node are actually the number of remote dofs
      // on the remote node.  This is a consistency check.
      //
      assert(NodalDofs[*sni-1].size() == RcvBuf[i][j]);

      RemoteNodalDofs[*sni-info.nlocal-1].resize(RcvBuf[i][j++]);
      sni++;
    }
  }  
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs node sizes received." << std::endl;
  
  // That phase is finished now, we can re-use the send and receive buffers
  // In the part where we exchange the *actual* global dof numbers 
  for(unsigned int i = 0;i < nborders;i++){
    RcvBuf[i].resize(totalsizes[i]);
    borders[i].remote_dofcount = totalsizes[i];
    // Go ahead and post the receives for the actual dofs
    comm.ARecv<Mesh::IndexType>(RcvBuf[i],borders[i].rpart-1); // async recv
    if(Out && verbose)
      *Out << "Posting receive for " << RcvBuf[i].size() << " dofs from " 
	   << borders[i].rpart << std::endl;
  }
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs receive posted 2." << std::endl;

  // Now pack and send off the actual global dof numbers for the
  // dofs that I own.
  for(unsigned int i = 0;i < nborders;i++){
    SndBuf[i].resize(0);
    std::vector<Mesh::IndexType>::iterator rni = borders[i].nrecv.begin();
    while(rni != borders[i].nrecv.end()){
      std::vector<Mesh::IndexType>::iterator ldi = NodalDofs[*rni-1].begin();
      while(ldi != NodalDofs[*rni-1].end()){
	SndBuf[i].push_back(*ldi+doffset);
	ldi++;
      }
      rni++;
    }
    if(Out && verbose){
      *Out << "Posting send of " << SndBuf[i].size() << " nodal dofs to " 
	   << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,SndBuf[i]," ");
      *Out << std::endl;
    }
    borders[i].recvsize = SndBuf[i].size();
    comm.ASend<Mesh::IndexType>(SndBuf[i],borders[i].rpart-1); // asynchronous send
  }
  comm.WaitAll();
  comm.ClearRequests();
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs send posted 2." << std::endl;

  // Unpack the received global dofs into the NodalDofs array
  for(unsigned int i = 0;i < nborders;i++){
    std::vector<Mesh::IndexType>::iterator rdi = RcvBuf[i].begin();
    std::vector<Mesh::IndexType>::iterator sni = borders[i].nsend.begin();
    if(Out && verbose){
      *Out << "Unpacking recv'd dofs:" << std::endl;
      DumpContents(*Out,RcvBuf[i]," ");
      *Out << std::endl;
    }
    while(sni != borders[i].nsend.end()){
      Mesh::IndexType rnindex = *sni - info.nlocal - 1;
      std::vector<Mesh::IndexType>::iterator ndi = RemoteNodalDofs[rnindex].begin();
      while(ndi != RemoteNodalDofs[rnindex].end())
	*ndi++ = *rdi++;
      sni++;
    }
  }
  comm.Barrier();
  
}

void BuildBorderStiffness(Mesh::Connectivity &ec,std::vector<Mesh::Border> &borders,
			  Mesh::Connectivity &NodalDofs,Mesh::Connectivity &ElementDofs,
			  Mesh::Connectivity &RemoteNodalDofs,
			  Mesh::PartInfo &info,Comm::CommunicatorObject &comm,
			  std::ostream *Out,bool verbose)
{
  unsigned int nborders = borders.size();
  // For each *node* that I own on the partition boundary, there will be
  // a number of stiffness entries coming from the remote process. Which is
  // some <number> * the number of dofs on the local node. This
  // first communication is to figure out that <number>.  The size of the
  // receive buffer here is the number of local nodes on the boundary.
  std::vector<std::vector<Mesh::IndexType> > RcvBuf(nborders);
  std::vector<Mesh::IndexType> rsizes(nborders,0);
  for(unsigned int i = 0;i < nborders;i++){
    rsizes[i] = borders[i].nrecv.size();
    RcvBuf[i].resize(rsizes[i]);
    if(Out && verbose)
      *Out << "Posting receives for " << rsizes[i] << " sizes of stiffness entries "
	   << "from " << borders[i].rpart << "." << std::endl;
    comm.ARecv<Mesh::IndexType>(RcvBuf[i],borders[i].rpart-1); // async recv
  }
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs have posted recvs." << std::endl;

  // Now calculate the *number of* stiffness entries I have for the remote 
  // nodal dofs.  The number of them will be the size of the unique list of 
  // local dofs that talk to the remote node.  We'll need the dual connectivity
  // to calculate this, so we build it here on the fly.
  Mesh::Connectivity dc;
  Mesh::IndexType number_of_nodes = NodalDofs.size();
  if(Out && verbose)
    *Out << "Forming dual connectivity for " << number_of_nodes << " nodes." << std::endl;
  ec.Inverse(dc,number_of_nodes);
  if(Out && verbose)
    *Out << "dual done." << std::endl;
  std::vector<std::vector<Mesh::IndexType> > SndBuf(nborders);
  std::vector<Mesh::IndexType> ssizes(nborders,0);
  for(unsigned int i = 0;i < nborders;i++){
    borders[i].data.SendAp.resize(0);
    borders[i].data.SendAi.resize(0);
    borders[i].data.SendAp.push_back(0);
    std::vector<Mesh::IndexType>::iterator sni = borders[i].nsend.begin();
    while(sni != borders[i].nsend.end()){
      // The list is the same for each border *node*
      std::list<Mesh::IndexType> doflist;
      // Count up the dof contributions from every local
      // element that touches the remote node
      std::vector<Mesh::IndexType>::iterator ei = dc[*sni-1].begin();
      while(ei != dc[*sni-1].end()){ // for every element touching the remote node
	std::vector<Mesh::IndexType>::iterator eni = ec[*ei-1].begin();
	while(eni != ec[*ei-1].end()){ // loop thru the element's nodes
	  std::vector<Mesh::IndexType>::iterator ndi;
	  std::vector<Mesh::IndexType>::iterator ndend;
	  Mesh::IndexType offset = 0;
	  if(*eni > info.nlocal){  // then it's a remote node
	    ndi   = RemoteNodalDofs[*eni - info.nlocal - 1].begin();
	    ndend = RemoteNodalDofs[*eni - info.nlocal -1].end();
	  }
	  else{
	    offset = info.doffset;
	    ndi   = NodalDofs[*eni-1].begin();
	    ndend = NodalDofs[*eni-1].end();
	  }
	  while(ndi != ndend){ // Loop thru the node's dofs and record them into the doflist
	    doflist.push_back(*ndi + offset);
	    ndi++;
	  }
	  eni++;
	}
	eni = ElementDofs[*ei-1].begin();
	while(eni != ElementDofs[*ei-1].end()){
	  doflist.push_back(*eni+info.doffset);
	  eni++;
	}
	ei++;					
      }
      doflist.sort();
      doflist.unique();
      // The row size for one dof is doflist.size(), however, you will be
      // sending that many entries for every DOF on the remote node - so
      // the total size is rowsize * ndof.  Let the remote processor deal 
      // with that detail. 
      SndBuf[i].push_back(doflist.size());
      borders[i].data.SendAp.push_back(doflist.size() + *(borders[i].data.SendAp.rbegin()));
      std::list<Mesh::IndexType>::iterator dli = doflist.begin();
      while(dli != doflist.end())
	borders[i].data.SendAi.push_back(*dli++);
      sni++;
    }
    if(Out && verbose){
      *Out << "SendAp for color " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.SendAp," ");
      *Out << std::endl
	   << "SendAi for color " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.SendAi," ");
      *Out << std::endl;
      *Out << "Posting send of " << SndBuf[i].size() << " stiffness row sizes "
	   << "to " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,SndBuf[i]," ");
      *Out << std::endl;
    }
    comm.ASend<Mesh::IndexType>(SndBuf[i],borders[i].rpart-1); // async recv
  }
  comm.WaitAll();
  if(verbose){
    *Out << "Receive Buffer: " << std::endl;
    for(unsigned int i = 0;i < nborders;i++){
      *Out << "Partition " << i+1 << ": ";
      DumpContents(*Out,RcvBuf[i]," ");
      *Out << std::endl;
    }
  }
  comm.ClearRequests();
  comm.Barrier();
  if(Out && verbose)
    *Out << "All procs have posted sends." << std::endl;
  // By now, the remote processor has sent us an array of sizes for the nz entries
  // relevant to each local *node*.  To get the total size of the stiffness block
  // that will be sent, we need to multiply by the number of actual dofs living 
  // on the local node.  Right now, we are just interested in the dof id's that 
  // we will be receiving/node.
  for(unsigned int i = 0;i < nborders;i++){
    borders[i].data.RecvAp.resize(RcvBuf[i].size()+1);
    std::vector<Mesh::IndexType>::iterator rbi = RcvBuf[i].begin();
    std::vector<Mesh::IndexType>::iterator rai = borders[i].data.RecvAp.begin();
    *rai = 0;
    while(rbi != RcvBuf[i].end()){
      Mesh::IndexType pv = *rai++;
      *rai = pv + *rbi++;
    }
    if(Out && verbose){
      *Out << "RecvAP from partition " << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.RecvAp," ");
      *Out << std::endl;
    }
    // We can go ahead and post receives for the Ai
    // This is not the way to determine total size any more, we've changed our
    // Ap array to the traditional representation [0 ..... nentries_Ai].
    //    while(rsi != borders[i].data.RecvAp.end())
    //      totalsize += *rsi++;
    // Now we just need the last value
    std::vector<Mesh::IndexType>::reverse_iterator rsi = borders[i].data.RecvAp.rbegin();
    unsigned int totalsize = *rsi;
    borders[i].data.RecvAi.resize(totalsize);
    if(Out && verbose)
      *Out << "Posting RecvAi (" << borders[i].data.RecvAi.size()
	   << ")  from partition" << borders[i].rpart << std::endl;
    comm.ARecv<Mesh::IndexType>(borders[i].data.RecvAi,borders[i].rpart-1); // async recv
    // Note how we've already conveniently filled the SendAi, send it out now
    comm.ASend<Mesh::IndexType>(borders[i].data.SendAi,borders[i].rpart-1); // async recv
    if(Out && verbose){
      *Out << "Posting SendAi (" << borders[i].data.SendAi.size() 
	   << ") to partition" << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.SendAi," ");
      *Out << std::endl;
    }
  }
  comm.WaitAll();
  comm.ClearRequests();
  comm.Barrier();
  if(Out && verbose){
    for(unsigned int i = 0;i < nborders;i++){
      *Out << "Recieved RecvAi (" << borders[i].data.RecvAi.size() 
	   << ") from partition" << borders[i].rpart << ":" << std::endl;
      DumpContents(*Out,borders[i].data.RecvAi," ");
      *Out << std::endl;
    }
  }
}

unsigned int AllocateCommBuffers(std::vector<Mesh::Border> &borders,
				 Mesh::Connectivity &NodalDofs,std::ostream *Out,bool verbose)
{
  unsigned int nborders = borders.size();
  unsigned int nalloc = 0;
  for(unsigned int i = 0;i < nborders;i++){
    unsigned int nval_send = 0;
    unsigned int nval_recv = 0;
    std::vector<Mesh::IndexType>::iterator ni = borders[i].nrecv.begin();
    std::vector<Mesh::IndexType>::iterator Api = borders[i].data.RecvAp.begin();
    while(ni != borders[i].nrecv.end()){
      Mesh::IndexType ndof = NodalDofs[*ni - 1].size();
      Mesh::IndexType begindex = *Api++;
      nval_recv += ((*Api-begindex) * ndof);
      //      Api++;
      ni++;
    }
    if(Out && verbose)
      *Out << "Allocating buffer(" << nval_recv << ") for receiving stiffness blocks"
	   << " from " << borders[i].rpart << "." << std::endl;
    borders[i].data.RecvBuffer.resize(nval_recv,0.0);
    borders[i].data.RecvBuffer.swap(borders[i].data.RecvBuffer);
    nalloc += nval_recv;
    ni  = borders[i].nsend.begin();
    Api = borders[i].data.SendAp.begin();
    // because the first val is 0
    while(ni != borders[i].nsend.end()){
      Mesh::IndexType ndof = NodalDofs[*ni - 1].size();
      Mesh::IndexType begindex = *Api++;
      nval_send += ((*Api-begindex) * ndof);
      //      Api++;
      ni++;
    }
    if(Out && verbose)
      *Out << "Allocating buffer(" << nval_send << ") for sending stiffness blocks"
	   << " to " << borders[i].rpart << "." << std::endl;
    borders[i].data.SendBuffer.resize(nval_send,1.0);
    borders[i].data.SendBuffer.swap(borders[i].data.SendBuffer);
    nalloc+=nval_send;
  }
  return(nalloc);
}

int RegisterCommBuffers(std::vector<Mesh::Border> &borders,Comm::CommunicatorObject &comm)
{
  unsigned int nborders = borders.size();
  for(unsigned int i = 0;i < nborders;i++)
    {
      if(comm.SetSend(borders[i].data.SendBuffer,borders[i].rpart-1) < 0)
	return 1;
      if(comm.SetRecv(borders[i].data.RecvBuffer,borders[i].rpart-1) < 0)
	return 1;
    }
  return(0);
}


void CountDofs(Mesh::Connectivity &ElementDofs,Mesh::Connectivity &NodalDofs,Mesh::PartInfo &info,
	       Mesh::IndexType &nglobaldof,Mesh::IndexType &nremotedof)
{
  Mesh::IndexType number_of_elements = ElementDofs.size();
  Mesh::IndexType number_of_nodes = NodalDofs.size();
  nglobaldof = 0;
  nremotedof = 0;
  for(Mesh::IndexType pp = 0;pp < number_of_elements;pp++){
    nglobaldof += ElementDofs[pp].size();
  }
  for(Mesh::IndexType pp = 0;pp < number_of_nodes;pp++){
    if((pp+1) <= info.nlocal)
      nglobaldof += NodalDofs[pp].size();
    else
      nremotedof += NodalDofs[pp].size();
  }
}

int
main(int argc,char *argv[])
{
  std::ostream* StdOut = NULL; // program output
  std::ostream* ErrOut = NULL; // program errors
  std::ostream* LogOut = NULL; // log for each proc (if debugging on)
  std::ofstream LogFile;
  Global::GlobalObj<std::string,Mesh::IndexType,Profiler::ProfilerObj> global;
  Comm::CommunicatorObject comm(&argc,&argv);
  unsigned int rank  = comm.Rank();
  unsigned int nproc = comm.Size();
  global.Init("MeshTest",rank);

  if(rank == 0){
    StdOut = &std::cout;
    ErrOut = &std::cerr;
  }
  global.SetDebugLevel(0);
  global.FunctionEntry("main");
  MeshTestComLine comline((const char **)argv);
  comline.Initialize();
  int clerr = comline.ProcessOptions();
  std::string sverb    =  comline.GetOption("verb");
  std::string spart    =  comline.GetOption("partition");
  std::string sclone   =  comline.GetOption("clone");
  bool do_partitioning = !comline.GetOption("partition").empty();
  bool do_mesh_testing = !comline.GetOption("mesh").empty();
  bool do_assembly     = !comline.GetOption("assembly").empty();
  bool do_reorient     = !comline.GetOption("reorient").empty();
  bool debug           = !comline.GetOption("debug").empty();
  bool array_checking  = !comline.GetOption("checking").empty();
  bool do_generate     = !comline.GetOption("generate").empty();
  bool do_clone        = !comline.GetOption("clone").empty();
  bool do_dump_sparse  = !comline.GetOption("sparse").empty();
  //  bool do_metis_test   = !comline.GetOption("metis").empty();
  bool do_renumbering  = !comline.GetOption("renumber").empty();
  Comm::DataTypes IndexDT = (sizeof(Mesh::IndexType) == sizeof(size_t) ?
			     Comm::DTSIZET : Comm::DTUINT);
  if(!comline.GetOption("help").empty()){
    if(StdOut) 
      *StdOut << comline.LongUsage() << std::endl;
    comm.SetExit(1);
  }
  if(comm.Check())
    return(0);
  
  if(clerr){
    if(ErrOut) 
      *ErrOut << comline.ErrorReport() << std::endl
	      << std::endl << comline.ShortUsage() << std::endl;
    comm.SetExit(1);
  }
  if(comm.Check())
    return(1);
  int verblevel = 0;
  if(!sverb.empty()){
    verblevel = 1;
    if(sverb != ".true."){
      std::istringstream Istr(sverb);
      Istr >> verblevel;
    }
  }
  Mesh::IndexType npartitions = 0;
  if(!spart.empty() && spart != ".true."){
    std::istringstream Istr(spart);
    Istr >> npartitions;
  } 
  if(!sclone.empty() && sclone != ".true."){
    std::istringstream Istr(sclone);
    Istr >> npartitions;
  }
  if(npartitions <= 0)
    npartitions = 2;
  if(nproc > 1){
    std::ostringstream debugfilename;
    debugfilename << comline.ProgramName() << ".output." << rank;
    LogFile.open(debugfilename.str().c_str());
    global.SetDebugStream(LogFile);
    LogOut = &LogFile;
    ErrOut = &LogFile;
    if(verblevel > 1 && rank > 0)
      StdOut = &LogFile;
  }
  else{
    global.SetDebugStream(std::cerr);
    LogOut=&std::cout;
    ErrOut=&std::cout;
  }
  if(verblevel > 1 || debug){
    global.SetDebugLevel(1);
  }
  if(debug && rank > 0)
    StdOut = LogOut;
  std::vector<std::string> infiles = comline.GetArgs();
  if(infiles.size()==0) {
    if(ErrOut) 
      *ErrOut << "MeshTest::Error: No input specified." << std::endl
	      << std::endl << comline.ShortUsage() << std::endl;
    comm.SetExit(1);
  }
  if(comm.Check())
    return(1);

  if(do_generate){
    std::string MeshName;
    if(!infiles.empty())
      MeshName = infiles[0];
    unsigned int N = 10;
    if(infiles.size() > 1){
      std::istringstream Instr(infiles[1]);
      Instr >> N;
    }
    if(N < 2)
      N = 2;
    Mesh::IndexType number_of_nodes = N*N*N;
    Mesh::IndexType number_of_elements = (N-1)*(N-1)*(N-1);
    Mesh::Connectivity ec;
    ec.resize(number_of_elements);
    if(StdOut && verblevel)
      *StdOut << "Creating cube mesh with " << number_of_nodes << " nodes, "
	      << number_of_elements << " elements." << std::endl;
    Mesh::NodalCoordinates nc(number_of_nodes);
    if(LogOut && debug)
      *LogOut << "NC.size() = " << nc.Size() << std::endl;
    Mesh::IndexType node_id = 1;
    Mesh::IndexType element_id = 1;
    for(unsigned int k = 0;k < N;k++){
      double z = k;
      for(unsigned int j = 0;j < N;j++){
	double y = j;
	for(unsigned int i = 0;i < N;i++){
	  double x = i;
	  nc.x(node_id) = x;
	  nc.y(node_id) = y;
	  nc.z(node_id++) = z;
	  if(k < (N-1) && j < (N-1) && i < (N-1)){
	    Mesh::IndexType ulnode = k*N*N + j*N + i%N + 1;
	    ec[element_id-1].push_back(ulnode);
	    ec[element_id-1].push_back(ulnode+N);
	    ec[element_id-1].push_back(ulnode+N+1);
	    ec[element_id-1].push_back(ulnode+1);
	    ec[element_id-1].push_back(ulnode+N*N);
	    ec[element_id-1].push_back(ulnode+N*N+N);
	    ec[element_id-1].push_back(ulnode+N*N+N+1);
	    ec[element_id-1].push_back(ulnode+N*N+1);
	    element_id++;
	  }
	}
      }
    }
    std::ofstream Ouf;
    Ouf.open(MeshName.c_str());
    Ouf << nc << std::endl << ec << std::endl;
    Ouf.close();
  }
	  
      
  if(do_mesh_testing){
    //    Mesh::MeshObject mesh;
    std::ifstream Inf;
    std::string MeshName;
    if(!infiles.empty()){
      MeshName = infiles[0];
    }
    std::ostringstream Ostr;
    Ostr << MeshName;
    if(nproc > 1)
      Ostr << "." << rank+1 << ".pmesh";
    Inf.open(Ostr.str().c_str());
    if(!Inf){
      if(ErrOut) 
	*ErrOut << "Failed to open " << Ostr.str() << std::endl;
      comm.SetExit(1);
    }
    if(comm.Check())
      return(1);
    
    Mesh::IndexType number_of_nodes = 0;
    Mesh::NodalCoordinates nc;
    global.FunctionEntry("ReadCoord");
    Inf >> nc;
    global.FunctionExit("ReadCoord");
    number_of_nodes = nc.Size();
    if(StdOut && verblevel) 
      *StdOut << "Number of nodes: " << number_of_nodes << std::endl;
    std::vector<double> bounds(6,0.0);
    Mesh::GetCoordinateBounds(nc,bounds);
    if(StdOut && verblevel)
      *StdOut << "Mesh Bounds: (" << bounds[0] << "," << bounds[1] 
	      << ") (" << bounds[2] << "," << bounds[3] 
	      << ") (" << bounds[4] << "," << bounds[5] << ")" << std::endl;
    Mesh::Connectivity econ;
    global.FunctionEntry("ReadCon");
    Inf >> econ;
    global.FunctionExit("ReadCon");
    Inf.close();
    Mesh::IndexType number_of_elements = econ.Nelem();
    if(StdOut && verblevel) 
      *StdOut << "Number of elements: " << number_of_elements
	      << std::endl;
    Mesh::IndexType N = 0;
    Mesh::IndexType N2 = 0;
    if(infiles.size() > 1){
      std::istringstream Istr(infiles[1]);
      Istr >> N;
    }
    if(infiles.size() > 2){
      std::istringstream Istr(infiles[2]);
      Istr >> N2;
    }
    if(false){
      Mesh::Connectivity::iterator ci = econ.begin();
      Mesh::IndexType el = 1;
      while(ci != econ.end()){
	if(StdOut && verblevel) 
	  *StdOut << "Element " << el++ << ": (";
	std::vector<Mesh::IndexType>::iterator ni = ci->begin();
	while(ni != ci->end()){
	  if(StdOut && verblevel) 
	    *StdOut << *ni++;
	  if(ni != ci->end())
	    if(StdOut && verblevel) 
	      *StdOut << ",";
	}
	if(StdOut && verblevel) 
	  *StdOut << ")\n";
	ci++;
      }
    }
    Mesh::Connectivity dc;
    global.FunctionEntry("DualCon");
    econ.Inverse(dc,number_of_nodes);
    global.FunctionExit("DualCon");
    //  Mesh::NeighborHood nl;
    //  global.FunctionEntry("Neighborhood");
    //  con.GetNeighborhood(nl,dc);
    //  global.FunctionExit("Neighborhood");
    bool do_neighborhood = false;
    if(do_neighborhood){
      Mesh::Connectivity nl2;
      // Neighborhood populates an array of neighbors
      // for each element.  A neighbor is any element
      // that shares a node with a given element.
      global.FunctionEntry("Neighborhood2");
      econ.GetNeighborhood(nl2,dc);
      global.FunctionExit("Neighborhood2");
      
      
      //  std::vector<std::list<Mesh::IndexType> > element_neighbors;
      //  global.FunctionEntry("CreateAdjEl Template");
      //  Mesh::AdjEList<std::vector<std::vector<Mesh::IndexType> >,
      //    std::vector<Mesh::IndexType> >
      //    (element_neighbors,dc,number_of_elements);  
      //  global.FunctionExit("CreateAdjEl Template");
      //  Mesh::IndexType nedges = 
      //    Mesh::NumberOfEdges< std::vector<std::list<Mesh::IndexType> >,
      //    std::list<Mesh::IndexType> >(anodelist);
      
      if(false){
	//     Mesh::NeighborHood::iterator ei = nl.begin();
	//     Mesh::IndexType el = 1;
	//     while(ei != nl.end()){
	//       if(StdOut && verblevel) *StdOut << "Element " << el++ << " neighbors: ";
	//       std::set<Mesh::IndexType>::iterator ni = ei->begin();
	//       while(ni != ei->end()){
	// 	if(StdOut && verblevel) *StdOut << *ni++;
	// 	if(ni != ei->end())
	// 	  if(StdOut && verblevel) *StdOut << ",";
	//       }
	//       if(StdOut && verblevel) *StdOut << std::endl;
	//       ei++;
	//     }
	//     std::vector<std::list<Mesh::IndexType> >::iterator eni = 
	//       element_neighbors.begin();
	//     Mesh::IndexType el = 1;
	//     while(eni  != element_neighbors.end()){
	//       if(StdOut && verblevel) *StdOut << "Element " << el++ << " neighbors: ";
	//       std::list<Mesh::IndexType>::iterator ni = eni->begin();
	//       while(ni != eni->end()){
	// 	if(StdOut && verblevel) *StdOut << *ni++;
	// 	if(ni != eni->end())
	// 	  if(StdOut && verblevel) *StdOut << ",";
	//       }
	//       if(StdOut && verblevel) *StdOut << std::endl;
	//       eni++;
      //     }
	Mesh::Connectivity::iterator ei2 = nl2.begin();
	Mesh::IndexType el = 1;
	while(ei2 != nl2.end()){
	  if(StdOut && verblevel) 
	    *StdOut << "Element " << el++ << " neighbors: ";
	  std::vector<Mesh::IndexType>::iterator ni = ei2->begin();
	  while(ni != ei2->end()){
	    if(StdOut && verblevel) 
	      *StdOut << *ni++;
	    if(ni != ei2->end())
	      if(StdOut && verblevel) 
		*StdOut << ",";
	  }
	  if(StdOut && verblevel) 
	    *StdOut << std::endl;
	  ei2++;
	}
      }
    }
    // Container Orientation Test
    if(false) {
      global.FunctionEntry("List Creation");
      std::vector<Mesh::IndexType> a;
      std::vector<Mesh::IndexType> b;
      a.resize(N);
      b.resize(N);
      for(Mesh::IndexType j = 1;j <= N;j++){
	a[j-1] = j;
	b[N-j] = (j+100)%N + 1;
      }
      global.FunctionExit("List Creation");
      global.FunctionEntry("Same Orientation");
      bool same_orientation = HaveSameOrientation(a,b);
      global.FunctionExit("Same Orientation");
      global.FunctionEntry("Opp Orientation");
      bool opp_orientation = HaveOppositeOrientation(a,b);
      global.FunctionExit("Opp Orientation");
      if(StdOut && verblevel) 
	*StdOut << "Lists did " << (same_orientation ? "" : "not ")
		<< "have same orientation.\n"
		<< "Lists did " << (opp_orientation ? "" : "not ")
		<< "have opposite orientation.\n";
    }
    global.FunctionEntry("Orient elements");
    Mesh::IndexType invertedcount = 0;
    for(Mesh::IndexType element_being_processed = 1;
 	element_being_processed <= number_of_elements;
	element_being_processed++)
      {
 	Mesh::IndexType index = element_being_processed - 1;
 	Mesh::IndexType size_of_element = econ[index].size();
 	Mesh::GenericElement ge(size_of_element);
 	if(ge.Inverted(econ[index],nc)){
 	  invertedcount++;
 	  ge.ReOrient(econ[index]);
 	}
	if(!ge.ShapeOK(econ[index],nc))
	  *StdOut << "Element " << element_being_processed 
		  << " has bad shape!" << std::endl;
      }
    global.FunctionExit("Orient elements");
    nc.destroy();
    if(StdOut && verblevel) 
      *StdOut << "Number of elements reoriented: " << invertedcount 
 	      << std::endl;
    global.FunctionEntry("All Face Processing");
    Mesh::IndexType nface_estimate = static_cast<Mesh::IndexType>(2.2 * number_of_elements);
    // This next loop actually populates the data
    Mesh::IndexType number_of_faces = 0;
    Mesh::Connectivity F_N;
    Mesh::Connectivity E_F;
    Mesh::Connectivity F_E;
    std::vector<Mesh::SymbolicFace> sf;
    econ.SyncSizes();
    global.FunctionEntry("GetFaceConnectivites");
    econ.BuildFaceConnectivity(F_N,E_F,sf,dc);
    global.FunctionExit("GetFaceConnectivites");
    global.FunctionEntry("Face dual");
    E_F.Inverse(F_E,F_N.Nelem());
    global.FunctionExit("Face dual");
    global.FunctionEntry("Canned shrink");
    F_N.ShrinkWrap();
    E_F.ShrinkWrap();
    F_E.ShrinkWrap();
    std::vector<Mesh::SymbolicFace>(sf).swap(sf);
    global.FunctionExit("Canned shrink");
    global.FunctionEntry("Face syncsizes");
    F_E.SyncSizes();
    F_N.SyncSizes();
    global.FunctionExit("Face syncsizes");
  
    global.FunctionEntry("Face Stats");
    Mesh::IndexType canned_vol_tri = 0;
    Mesh::IndexType canned_vol_quad = 0;
    Mesh::IndexType canned_bndry_tri = 0;
    Mesh::IndexType canned_bndry_quad = 0;
    Mesh::Connectivity::iterator fei = F_E.begin();
    while(fei != F_E.end()){
      Mesh::IndexType faceid = (fei++ - F_E.begin()) + 1;
      if(F_E.Esize(faceid) == 1){
	if(F_N.Esize(faceid) == 3)
	  canned_bndry_tri++;
	else
	  canned_bndry_quad++;
      }
      else{
	if(F_N.Esize(faceid) == 3)
	  canned_vol_tri++;
	else
	  canned_vol_quad++;
      }	
    }
    global.FunctionExit("Face Stats");
    if(StdOut && verblevel) 
      *StdOut << "Estimated number of faces: " << nface_estimate << std::endl
	      << "Actual Number of faces: " << F_N.Nelem() << "\n"
	      << "Boundary faces: " << canned_bndry_tri + canned_bndry_quad << "\n"
	      << "Volume faces: " << canned_vol_tri + canned_vol_quad << "\n"
	      <<  canned_vol_quad +canned_bndry_quad << " Quads.\n"
	      << canned_vol_tri + canned_bndry_tri  << " Tris.\n"
	      << canned_bndry_quad+canned_bndry_tri  << " Boundary faces (" 
	      << canned_bndry_quad << "," << canned_bndry_tri << ")" << std::endl; 
    global.FunctionExit("All Face Processing");
    //    Mesh::Connectivity facecon_dual;
    //    global.FunctionEntry("Face dual");
    //    F_N.Inverse(facecon_dual,number_of_nodes);
    //    global.FunctionExit("Face dual");
    global.FunctionEntry("Edges processing");
    // For each node, determine which nodes are connected by
    // traversing the faces containing the node.
    std::vector<std::list<Mesh::IndexType> > anodelist;
    global.FunctionEntry("CreateAdjNodes Template");
    CreateAdjacentNodeList<std::vector<std::vector<Mesh::IndexType> >,
      std::vector<Mesh::IndexType> >
      (anodelist,F_N,number_of_nodes);  
    global.FunctionExit("CreateAdjNodes Template");
    Mesh::IndexType nedges = 
      NumberOfEdges< std::vector<std::list<Mesh::IndexType> >,
      std::list<Mesh::IndexType> >(anodelist);
    global.FunctionExit("Edges processing");
    if(StdOut && verblevel) 
      *StdOut << "Number of edges: " << nedges << std::endl;
    if(false){ // report the edges
      std::vector< std::list<Mesh::IndexType> >::iterator anli = anodelist.begin();
      Mesh::IndexType ii = 0;
      while(anli != anodelist.end())
	{
	  if(StdOut && array_checking) 
	    *StdOut << "Node " << ++ii << " neighbors:";
	  std::list<Mesh::IndexType>::iterator ni = anli->begin();
	  while(ni != anli->end())
	    {
	      if(StdOut && array_checking) 
		*StdOut << *ni++;
	      if(ni != anli->end())
		if(StdOut && array_checking) 
		  *StdOut << ",";
	    }
	  if(StdOut && array_checking) 
	    *StdOut << std::endl;
	  anli++;
	}
    }
  }
  
  comm.Barrier();

  // Comment
  if(do_assembly){
    global.FunctionEntry("Assembly Function");
    // Initial stab at fake assembly
    global.FunctionEntry("Initialization");
    //    Mesh::NodalCoordinates nc;
    Mesh::Connectivity econ;
    Mesh::PartInfo info;
    Mesh::IndexType datasize = sizeof(Mesh::IndexType);
    Mesh::IndexType number_of_nodes = 0;
    Mesh::IndexType number_of_elements = 0;
    std::ifstream Inf;
    if(LogOut && debug)
      *LogOut << "DataSize = " << datasize << std::endl;
    if(true) {
      global.FunctionEntry("ReadCoord");
      std::string MeshName;
      if(!infiles.empty())
	MeshName = infiles[0];
      std::ostringstream FNOstr;
      if(nproc > 1){
	FNOstr << MeshName << "." << rank+1 << ".info";
	Inf.open(FNOstr.str().c_str());
	if(!Inf){
	  if(ErrOut) 
	    *ErrOut << "Failed to open " << FNOstr.str() << std::endl;
	  comm.SetExit(1);
	}
	if(comm.Check())
	  return(1);
	Inf >> info.npart >> info.part >> info.nelem >> info.nnodes >> info.nborder >> info.nshared >> info.nown;
	info.nlocal = info.nnodes - info.nshared + info.nown;
	Inf.close();
	FNOstr.str("");
	FNOstr << MeshName << "." << rank+1 << ".pmesh";
      }
      else{
	info.nborder = 0;
	info.nshared = 0;
	FNOstr << MeshName;
      }
      Inf.open(FNOstr.str().c_str());
      if(!Inf){
	if(ErrOut) 
	  *ErrOut << "Failed to open " << FNOstr.str() << std::endl;
	comm.SetExit(1);
      }
      if(comm.Check())
	return(1);
      if(true){
	Mesh::NodalCoordinates nc;
	Inf >> nc;
	number_of_nodes = nc.Size();
      }
      global.FunctionExit("ReadCoord");
      global.FunctionEntry("ReadCon");
      Inf >> econ;
      econ.ShrinkWrap();
      global.FunctionExit("ReadCon");
      comm.Barrier();
    }
    number_of_elements = econ.Nelem();
    if(LogOut && verblevel > 1){
      Mesh::IndexType consize = GetTotalSize(econ);
      *LogOut << "Connectivity size(bytes): " << (consize+(2*number_of_elements))*datasize << std::endl;
    }
    if(nproc == 1){
      info.nnodes = number_of_nodes;
      info.nown = 0;
      info.nlocal = number_of_nodes;
      info.nelem = number_of_elements;
    }
    if(StdOut && verblevel) 
      *StdOut << "Number of elements: " << number_of_elements
	      << std::endl
	      << "Number of nodes: " << number_of_nodes << std::endl
	      << "Number of local nodes: " << info.nlocal << std::endl;
    std::vector<Mesh::Border> borders;
    std::list<Mesh::IndexType> remotecolorlist;
    Mesh::Connectivity bordernodes;
    Mesh::Connectivity gbordernodes;
    if(nproc > 1){
      Inf >> info.nborder;
      if(StdOut && verblevel)
	*StdOut << "Number of neighboring partitions: " << info.nborder << std::endl;
      borders.resize(info.nborder);
      bordernodes.resize(info.nborder);
      gbordernodes.resize(info.nborder);
      Mesh::IndexType totalrecv = 0;
      Mesh::IndexType totalsend = 0;
      for(Mesh::IndexType nn = 0;nn < info.nborder;nn++){
	Mesh::IndexType nrecv = 0;
	Mesh::IndexType nsend = 0;
	Inf >> borders[nn].rpart >> nrecv >> nsend;
	totalrecv += nrecv;
	totalsend += nsend;
	remotecolorlist.push_back(borders[nn].rpart);
	for(Mesh::IndexType ii = 0;ii < nrecv;ii++){
	  Mesh::IndexType nodeid = 0;
	  Inf >> nodeid;
	  borders[nn].nrecv.push_back(nodeid);
	  bordernodes[nn].push_back(nodeid);
	  Inf >> nodeid;
	  gbordernodes[nn].push_back(nodeid);
	}
	for(Mesh::IndexType ii = 0;ii < nsend;ii++){
	  Mesh::IndexType nodeid = 0;
	  Inf >> nodeid;
	  borders[nn].nsend.push_back(nodeid);
	  bordernodes[nn].push_back(nodeid);
	  Inf >> nodeid;
	  gbordernodes[nn].push_back(nodeid);
	}
      }
      if(LogOut && verblevel > 1){
	*LogOut << "Border description storage size(bytes): " 
		<< (totalrecv+totalsend)*datasize
		<< std::endl;
      }
    }
    Inf.close();
    if(LogOut && debug && array_checking)
      *LogOut << "BorderNodes:" << std::endl
	      << bordernodes << std::endl
	      << "GlobalBorderNodes:" << std::endl
	      << gbordernodes << std::endl;
    comm.Barrier();
    Mesh::IndexType N = 0;
    Mesh::IndexType N2 = 0;
    if(infiles.size() > 1){
      std::istringstream Istr(infiles[1]);
      Istr >> N;
      if(N < 0)
	N = 0;
    }
    if(infiles.size() > 2){
      std::istringstream Istr(infiles[2]);
      Istr >> N2;
      if(N2 < 0)
	N2 = 0;
    }
    comm.Barrier();
    if(LogOut && verblevel > 1)
      *LogOut << "All command line args processed." << std::endl;
    if(StdOut && verblevel){
      *StdOut << "Number of nodal dofs: " << N << std::endl
	      << "Number of element dofs: " << N2 << std::endl;
    }
    if(LogOut && verblevel > 1)
      *LogOut << "Building dual connectivity..." << std::endl;
    Mesh::Connectivity dc;
    global.FunctionEntry("DualCon");
    econ.Inverse(dc,number_of_nodes);
    global.FunctionExit("DualCon");
    dc.ShrinkWrap();
    Mesh::Connectivity nl2;
    comm.Barrier();
    if(LogOut && debug){
      *LogOut << "Dual Con size(bytes):" 
	      << (GetTotalSize(dc)+2*dc.size())*datasize
	      << std::endl;
      if(array_checking)
	*LogOut << "dual connectivity:" << std::endl
		<< dc << std::endl;
    }
    // Neighborhood populates an array of neighbors
    // for each element.  A neighbor is any element
    // that shares a node with a given element.
    if(LogOut && verblevel > 1)
      *LogOut << "Building element neighborhood." << std::endl;
    global.FunctionEntry("Neighborhood");
    econ.GetNeighborhood(nl2,dc);
    global.FunctionExit("Neighborhood"); 
    nl2.ShrinkWrap();
    if(LogOut && array_checking)
      *LogOut << "Element Neighborhood:" << std::endl
	      << nl2 << std::endl; 
    comm.Barrier();
    if(LogOut && debug)
      *LogOut << "Neighborhood size(bytes): " 
	      << (GetTotalSize(nl2)+(2*nl2.size()))*datasize
	      << std::endl;
    Mesh::Connectivity NodalDofs;
    NodalDofs.resize(number_of_nodes);
    NodalDofs.Sync();
    Mesh::Connectivity ElementDofs;
    ElementDofs.resize(number_of_elements);
    ElementDofs.Sync();
    std::vector<FEM::FieldData> fields;
    Mesh::IndexType nfields = 0;
    if(N > 0){
      FEM::FieldData nfield;
      nfield.Order(1);
      nfield.Components(N);
      fields.push_back(nfield);
      nfields++;
    }
    if(N2 > 0){
      FEM::FieldData efield;
      efield.Order(0);
      efield.Components(N2);
      fields.push_back(efield);
      nfields++;
    }

    // Assigns dof numbers to each node and element into the provided
    // arrays.
    //  global.FunctionEntry("Uniform DOF Assignment");
    //  Mesh::IndexType ndoftotal = Mesh::AssignUniformDofs(econ,dc,ElementFields,
    //						NodalFields,NodalDofs,
    //						ElementDofs,fields);
    //  global.FunctionExit("Uniform DOF Assignment");
    comm.Barrier();
    if(LogOut && verblevel > 1)
      *LogOut << "All processors finished prepping."
	      << std::endl
	      << "Now assigning DOFs." << std::endl;
    global.FunctionExit("Initialization");

    Mesh::IndexType ndoftotal = 0;
    Mesh::IndexType nlocal_dofs = 0;
    std::list<Mesh::IndexType> border_elements;
    //    border_elements.resize(info.nborder);
    Mesh::IndexType nvolume_nodes = number_of_nodes - info.nshared;
    //    Mesh::IndexType nlocal_nodes = nvolume_nodes + info.nown;
    std::vector<bool> is_border_element(number_of_elements,false);
    std::vector<bool> element_processed(number_of_elements,false);
    if(LogOut && verblevel)
      *LogOut << "Creating border element list for " << info.nshared
	      << " border nodes and " << nvolume_nodes << " volume nodes."
	      << std::endl;
    global.FunctionEntry("Assembly Prep");
    global.FunctionEntry("Find Border Elements");
    // At this time, we want border elements to be all of those which 
    // have entities on the partition boundary. - try something else
    for(Mesh::IndexType nn = nvolume_nodes;nn < number_of_nodes;nn++){
      //    for(Mesh::IndexType nn = nlocal_nodes;nn < number_of_nodes;nn++){
      std::vector<Mesh::IndexType>::iterator ei = dc[nn].begin();
      while(ei != dc[nn].end()){
	if(!element_processed[*ei-1]){
	  element_processed[*ei-1] = true;
	  border_elements.push_back(*ei);
	  is_border_element[*ei-1] = true;
	}
	ei++;
      }
    }
    border_elements.sort();
    global.FunctionExit("Find Border Elements");
    if(LogOut && verblevel){
      *LogOut << "Found " << border_elements.size() << " border elements."
	      << std::endl;
      if(array_checking){
	DumpContents(*LogOut,border_elements," ");
	*LogOut << std::endl;
      }
    }
    std::list<Mesh::IndexType> queue;
    if(true){
      Mesh::IndexType nprocessed = 0;
      std::vector<bool> neighbors_processed(number_of_elements,false);
      element_processed.clear();
      element_processed.resize(number_of_elements,false);
      queue.clear();
      std::list<Mesh::IndexType>::iterator qi = queue.begin();
      global.FunctionEntry("Element order queue");
      while(nprocessed < number_of_elements){
	bool from_queue = false;
	Mesh::IndexType elnum = 0;
	if(queue.empty()){
	  queue.push_back(1);
	  qi = queue.begin();
	}
	if(qi != queue.end()){
	  elnum = *qi-1;
	  from_queue = true;
	}
	else { // strangely not done, get next unprocessed element
	  for(Mesh::IndexType i = 0;i < number_of_elements;i++)
	    if(element_processed[i])
	      elnum++;
	}
	if(!element_processed[elnum]){
	  element_processed[elnum] = true;
	  if(!from_queue){
	    queue.push_back(elnum+1);
	  }
	  nprocessed++;
	}
	if(!neighbors_processed[elnum]){
	  neighbors_processed[elnum] = true;
	  std::vector<Mesh::IndexType>::iterator ni = nl2[elnum].begin();
	  while(ni != nl2[elnum].end()){
	    if(!element_processed[*ni-1]){
	      element_processed[*ni-1] = true;
	      queue.push_back(*ni);
	      nprocessed++;
	    }
	    ni++;
	  }
	}
	qi++;
      }
      global.FunctionExit("Element order queue");
      if(LogOut && debug && array_checking){
	*LogOut << "Element order queue:" << std::endl;
	DumpContents(*LogOut,queue," ");
	*LogOut << std::endl;
      }
    }
    nl2.destroy();
    Mesh::IndexType nvoldof     = 1;
    Mesh::IndexType nborderdof  = 1;
    Mesh::IndexType nshareddof  = 0;
    Mesh::IndexType nremotedof  = 0;
    Mesh::IndexType ndoftot = 0;
    global.FunctionEntry("Dof Assignment");
    if(true){
      if(nproc == 1 && false){
	std::vector<bool> node_processed(number_of_nodes,false);
	std::list<Mesh::IndexType>::iterator ei = queue.begin();
	while(ei != queue.end()){
	  std::vector<Mesh::IndexType>::iterator ni = econ[*ei-1].begin();
	  while(ni != econ[*ei-1].end()){
	    if(!node_processed[*ni-1]){
	      node_processed[*ni-1] = true;
	      NodalDofs[*ni-1].resize(N);
	      for(unsigned int i = 0;i < N;i++)
		NodalDofs[*ni-1][i] = nvoldof++;
	    }
	    ni++;
	  }
	  ElementDofs[*ei-1].resize(N2);
	  for(unsigned int i = 0;i < N2;i++)
	    ElementDofs[*ei-1][i] = nvoldof++;
	  ei++;
	}
      }
      else { // nproc > 1 - this loop works for nproc == 1 now too
	std::vector<bool> node_processed(number_of_nodes,false);
	std::list<Mesh::IndexType>::iterator ei = queue.begin();
	while(ei != queue.end()){
	  if(!is_border_element[*ei-1]){ 
	    if(array_checking)
	      *LogOut << "Assigning dofs for nonborder element " << *ei << std::endl;
	    std::vector<Mesh::IndexType>::iterator ni = econ[*ei-1].begin();
 	    while(ni != econ[*ei-1].end()){
 	      if(!node_processed[*ni-1]){
 		node_processed[*ni-1] = true;
 		NodalDofs[*ni-1].resize(N);
 		for(unsigned int i = 0;i < N;i++){
		  if(array_checking)
		    *LogOut << "node " << *ni << " getting dof " << nvoldof << std::endl;
 		  NodalDofs[*ni-1][i] = nvoldof++;
		}
 	      }
 	      ni++;
 	    }
	    ElementDofs[*ei-1].resize(N2);
 	    for(unsigned int i = 0;i < N2;i++){
	      if(array_checking)
		*LogOut << "element " << *ei << " getting dof " << nvoldof << std::endl;
 	      ElementDofs[*ei-1][i] = nvoldof++;
	    }
 	  }
	  ei++;
	}
	// 1st of three passes through border elements first picks up the 
	// elements
	ei = border_elements.begin();
	while(ei != border_elements.end()){
	  ElementDofs[*ei-1].resize(N2);
	  for(unsigned int i = 0;i < N2;i++){
	    if(array_checking)
	      *LogOut << "border element " << *ei << " getting dof " << nvoldof << std::endl;
	    ElementDofs[*ei-1][i] = nvoldof++;
	  }
	  ei++;
	}
	nvoldof--;
	// 2nd loop picks up all local nodes
	ei = border_elements.begin();
	while(ei != border_elements.end()){
	  std::vector<Mesh::IndexType>::iterator ni = econ[*ei-1].begin();
	  while(ni != econ[*ei-1].end()){
	    if(!node_processed[*ni-1] && *ni <= info.nlocal){
	      node_processed[*ni-1] = true;
	      NodalDofs[*ni-1].resize(N);
	      for(unsigned int i = 0;i < N;i++){
		if(array_checking)
		  *LogOut << "owned border node " << *ni << " getting dof " 
			  << nvoldof+nborderdof << std::endl;
		NodalDofs[*ni-1][i] = (nvoldof + nborderdof++);
	      }
	    }
	    ni++;
	  }
	  ei++;
	}
	nshareddof = nborderdof - 1;
	// 3rd loop puts tentative numbers on remote nodes
	// Note that the local doffing of remotely owned 
	// entity has potential to not match across the
	// partition boundary.  We will sync them up later
	// when we do remote dof discovery.
	ei = border_elements.begin();
	while(ei != border_elements.end()){
	  std::vector<Mesh::IndexType>::iterator ni = econ[*ei-1].begin();
	  while(ni != econ[*ei-1].end()){
	    if(!node_processed[*ni-1]){
	      node_processed[*ni-1] = true;
	      NodalDofs[*ni-1].resize(N);
	      for(unsigned int i = 0;i < N;i++){
		if(array_checking)
		  *LogOut << "remote border node " << *ni << " getting local dof " 
			  << nvoldof+nborderdof << std::endl;
		NodalDofs[*ni-1][i] = (nvoldof + nborderdof++);
	      }
	    }
	    ni++;
	  }
	  ei++;
	}
	nborderdof--;
	ndoftot = nvoldof + nborderdof;
	nremotedof = nborderdof - nshareddof;
	nlocal_dofs = nvoldof + nshareddof;
	if(StdOut && verblevel)
	  *StdOut << "NDof TOTAL (local + adjacent remote): " << ndoftot << std::endl
		  << "NDof Local  (local total) DOFS:  "  << nlocal_dofs << std::endl
		  << "NDof Remote (initial count of dofs owned remotely)  DOFS: "  
		  << nremotedof  << std::endl
		  << "NDof Border (total of dofs on partition boundaries) DOFS: "  
		  << nborderdof  << std::endl
		  << "NDof Shared (dofs on the partition boundary that I own) DOFS: " 
		  << nshareddof << std::endl;
	if(LogOut && debug){
	  *LogOut << "nborderelems = " << border_elements.size() << std::endl;
	  border_elements.sort();
	  border_elements.unique();
	  *LogOut << "after sort, size = " << border_elements.size() << std::endl;
	}
      }
    }
    global.FunctionExit("Dof Assignment");
    if(StdOut && verblevel) 
      *StdOut << "Number of local DOFs: " << nlocal_dofs << std::endl
	      << "Number of Elements on partition boundaries: " << border_elements.size()
	      << std::endl;
    if(LogOut && debug){
      *LogOut << "Total NodalDof storage(bytes): "
	      << (GetTotalSize(NodalDofs)+(number_of_nodes*2))*datasize
	      << std::endl
	      << "Total ElementDof storage(bytes): "
	      << (GetTotalSize(ElementDofs)+(number_of_elements*2))*datasize
	      << std::endl;
    }
    comm.Barrier();
    
    NodalDofs.Sync();
    NodalDofs.SyncSizes();
    ElementDofs.Sync();
    ElementDofs.SyncSizes();
    if(LogOut && debug && array_checking)
      *LogOut << "NodalDofs: " << std::endl
	      << NodalDofs << std::endl
	      << "ElementDofs: " << std::endl
	      << ElementDofs << std::endl;
    

    // Initial global dof counting
    // We need to shuffle things around here when we have
    // disparity in dofs on shared mesh entities.  Right now
    // we know there is no disparity - so we just do the 
    // straightforward thing.
    comm.Barrier();
    if(LogOut && debug)
      *LogOut << "All processors entering global dof counting." 
	      << std::endl;
    global.FunctionEntry("Global DOF counting");
    Mesh::IndexType nglobal_dofs = nlocal_dofs;
    std::vector<Mesh::IndexType> ngdofs_p(nproc,0);
    std::vector<Mesh::IndexType> nnbr_p(nproc,0);
    ngdofs_p[rank] = nglobal_dofs;
    ndoftotal = ndoftot;
    Mesh::IndexType nremote_dofs = ndoftotal - nglobal_dofs;  
    comm.Barrier();
    if(LogOut && debug)
      *LogOut << "All processors done counting local dofs."
	      << std::endl
	      << "Now gathering ..." << std::endl;
    int commstat = comm.AllGather(nglobal_dofs,ngdofs_p);
    comm.Barrier();
    if(LogOut && debug){
      *LogOut << "Global Dofcount verification = " << nglobal_dofs << std::endl
	      << "commstat = " << commstat << std::endl;
      if(array_checking){
	*LogOut << "NGDofs:" << std::endl;
	DumpContents(*LogOut,ngdofs_p);
	*LogOut << std::endl;
      }
    }
   
    info.doffset = 0;
    Mesh::IndexType ndof_global_total = 0;
    for(unsigned int nn1 = 0;nn1 < rank;nn1++)
      info.doffset += ngdofs_p[nn1];
    for(unsigned int nn2 = 0;nn2 < nproc;nn2++)
      ndof_global_total += ngdofs_p[nn2];
    if(StdOut && verblevel)
      *StdOut << "Number of global dofs on this proc: " << nglobal_dofs  << std::endl
	      << "Number of remote dofs on this proc: " << nremote_dofs  << std::endl
	      << "Number of global dofs over all procs: " << ndof_global_total << std::endl;
    if(StdOut && verblevel)
      *StdOut << "My global offset: " << info.doffset << std::endl;
    // We have enough info to form a message to all neighbors to let 
    // them know the global dof numbers of border nodes that I own.
    // Loop thru all the borders and do it
    comm.Barrier();
    if(LogOut && debug)
      *LogOut << "All procs entering global dof exchange." 
	      << std::endl;

    // This section of code does remote dof discovery.  That is, 
    // it let's all of our neighbors know the global dof id's of 
    // the mesh entities that we own on the partition boundary.
    // Likewise, we will receive the equivalent dof information from 
    // the remote partitions.
    Mesh::IndexType nremote_nodes = info.nnodes - info.nlocal;
    Mesh::Connectivity RemoteNodalDofs;
    if(LogOut && verblevel)
      *LogOut << "Number of remote nodes: " << nremote_nodes << std::endl;
    RemoteNodalDofs.resize(nremote_nodes);
    RemoteNodalDofs.Sync();
    //    TestMPI(borders,NodalDofs,RemoteNodalDofs,info,comm,LogOut,array_checking);
    //    LogFile.close();
    //    comm.Barrier();
    //    sleep(10);
    //    comm.Barrier();
    //    comm.SetExit(1);
    //    if(comm.Check())
    //      return(1);
    //    else{
    //      *StdOut << "SLEEPING" << std::endl;
    //      sleep(1000000);
    //    }
    global.FunctionEntry("Nodal DOF Exchange");
    GlobalDofExchange(borders,NodalDofs,RemoteNodalDofs,info,comm,LogOut,array_checking);
    global.FunctionExit("Nodal DOF Exchange");
    RemoteNodalDofs.ShrinkWrap();
    RemoteNodalDofs.SyncSizes();
    comm.Barrier();

    global.FunctionEntry("Count Dofs");
    Mesh::IndexType nglobdof = 0;
    Mesh::IndexType nremotedof2 = 0;
    CountDofs(ElementDofs,NodalDofs,info,nglobdof,nremotedof2);
    global.FunctionExit("Count Dofs");
    if(StdOut && verblevel)
      *StdOut << "After global exchange and recount: " << std::endl 
	      << "Number of global dofs on this proc: " << nglobdof  << std::endl
	      << "Number of remote dofs on this proc: " << nremotedof2  << std::endl;

    assert(nglobdof == nglobal_dofs);
    //    nremote_dofs = nremotedof;
    comm.Barrier();
    // Now we need to build the send and receive stiffnesses
    // The new way is to assume that the ElementDofs actually has ONLY elemental dofs, not 
    // with the nodal dofs stacked on as previously. (it only hurts performance at this point)
    global.FunctionEntry("Build Communicated Stiffness");
    BuildBorderStiffness(econ,borders,NodalDofs,ElementDofs,RemoteNodalDofs,info,comm,LogOut,array_checking);
    global.FunctionExit("Build Communicated Stiffness");
    // Allocate the communication buffers
    unsigned int nalloc = AllocateCommBuffers(borders,NodalDofs,LogOut,debug);
    if(LogOut && debug)
      *LogOut << "Allocated " << nalloc*sizeof(double) << " bytes of communication buffers." 
	      << std::endl;
    comm.Barrier();
    if(LogOut && verblevel)
      *LogOut << "Registering Comm Buffers." << std::endl;
    // Register the communication buffers with the communication object
    if(RegisterCommBuffers(borders,comm)){
      if(ErrOut)
	*ErrOut << "Error setting registering the communication buffers."
		<< std::endl;
      comm.SetExit(1);
    }
    if(comm.Check())
      return(1);    
    global.FunctionExit("Global DOF counting");
    comm.Barrier();
    if(LogOut && debug)
      *LogOut << "All processors have exchanged interpartition DOF information."
	      << std::endl
	      << "Proceeding with dof sorting." << std::endl;

    
    global.FunctionEntry("Build Sparsity");
    // Now to actually build the sparsity pattern, LD[GD]
    // We have E[LED], N[LND], Nr[GD], LD:GD
    // We'll start by building a local version of K, LD[LD],
    // and using a straightforward mapping LD:GD to obtain LD[GD].
    // Then we'll add any needed global dofs coming from remote
    // processors.
    //
    // 1. Build local K, LD[GD]
    //    - Need to somehow merge E[LED] and Nl[LND], since
    //      we always assemble element-by-element, we should
    //      rack up the LD on the elements in the order of
    //      all nodal dofs first, then element dofs. 
    //
    global.FunctionEntry("Create E[LD]");
    Mesh::Connectivity E_LD;
    E_LD.resize(number_of_elements);
    for(Mesh::IndexType en = 0;en < number_of_elements;en++){
      std::vector<Mesh::IndexType>::iterator eni = econ[en].begin();
      while(eni != econ[en].end()){ // for every node in this element
	Mesh::IndexType node_id = *eni++;
	std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[node_id-1].begin();
	while(ndi != NodalDofs[node_id-1].end())
	  E_LD[en].push_back(*ndi++); // add it's dofs to the list of dofs for this element
	
      }
      std::vector<Mesh::IndexType>::iterator edi = ElementDofs[en].begin();
      while(edi != ElementDofs[en].end())
	E_LD[en].push_back(*edi++); // add the element's dofs as well
    }
    global.FunctionEntry("E[LD] post");
    E_LD.ShrinkWrap();
    E_LD.Sync();
    E_LD.SyncSizes();
    // We need the sizes later
    std::vector<Mesh::IndexType> element_ndofs(number_of_elements,0);
    for(Mesh::IndexType ei = 0;ei < number_of_elements;ei++)
      element_ndofs[ei] = E_LD.Esize(ei+1);
    global.FunctionExit("E[LD] post");
    global.FunctionExit("Create E[LD]");
    //
    // Now E_LD is a list of all dofs for every element.  It is an important point
    // that LD are the LD (i.e. the local dofs), because they run from 1-nlocal_dofs.
    // Since both the indices of E (1-number_of_elements) and dofs both are sequential
    // lists, we can do some nifty array manipulations to obtain the information we need.
    // First, we invert E[LD] to obtain LD[E]:
    Mesh::Connectivity LD_E;
    global.FunctionEntry("Invert E[LD]");
    E_LD.Inverse(LD_E,ndoftot);
    global.FunctionExit("Invert E[LD]");
    //
    // Now, there is one slight problem with our version of LD[E], in that it includes
    // LD which are actually just local numbers for remote dofs.  We do not need the
    // stiffness rows for these.  In fact, we already have created the stiffness blocks 
    // for these rows when we created the communication buffers.  We know those dofs are 
    // all dofs with id's > nglobdof.  If we don't take care of the problem now, we pay
    // a hefty memory and processing fee in the following steps.
    LD_E.Sync();
    assert(LD_E.Nelem() == ndoftot);
    for(Mesh::IndexType di = nglobdof;di < ndoftot;di++)
      LD_E[di].resize(0);
    LD_E.ShrinkWrap();
    LD_E.SyncSizes();

    // LD[E]%E[LD] = LD[LD] := Stiffness
    Mesh::Connectivity symbstiff;
    global.FunctionEntry("Sparsity Pattern");
    LD_E.GetNeighborhood(symbstiff,E_LD,false,false);
    global.FunctionExit("Sparsity Pattern");

    comm.Barrier();
    LD_E.destroy();
    E_LD.destroy();

    if(LogOut && debug && array_checking)
      *LogOut << "All Local K:" << std::endl << symbstiff
	      << std::endl;
    
    
    // Now we need to do two things, update our "K" with the dofs which
    // are coming from remote procs and also map all the dof id's to 
    // global.  Then we're done building the parallel version of K.
    //
    // First, we need a map from local dof id's for the remote dofs 
    // on the partition boundary to the global id's that have been
    // told to us by the remote processor.  
    //
    // *Note:  This code updates the local dof rows which already
    // contain local dof numbers for remotely owned dofs.  These 
    // remote dofs are just the dofs of the remotely owned border
    // nodes.  This is just a mapping operation.
    //
    global.FunctionEntry("K update");
    if(true){
      std::vector<Mesh::IndexType> local_dof_to_global;
      local_dof_to_global.resize(ndoftot-nglobdof);
      for(Mesh::IndexType ni = info.nlocal;ni < number_of_nodes;ni++){
	std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[ni].begin();
	std::vector<Mesh::IndexType>::iterator rdi = RemoteNodalDofs[ni-info.nlocal].begin();
	assert(NodalDofs[ni].size() == RemoteNodalDofs[ni-info.nlocal].size());
	while(ndi != NodalDofs[ni].end()){
	  local_dof_to_global[*ndi-nglobdof-1] = *rdi;
	  ndi++;
	  rdi++;
	}
      }
      Mesh::Connectivity::iterator ConIt = symbstiff.begin();
      while(ConIt != symbstiff.end()){
	std::vector<Mesh::IndexType>::iterator  dofIt = ConIt->begin();
	Mesh::IndexType ind = 0;
	while(dofIt != ConIt->end()){
	  Mesh::IndexType dof_id = *dofIt++;
	  if(dof_id <= nglobdof)
	    (*ConIt)[ind++] = dof_id+info.doffset;
	  else
	    (*ConIt)[ind++] = local_dof_to_global[dof_id - nglobdof - 1];
	}
	ConIt++;
      }
    }
    //
    // Now update with the remote dof info.
    // *Note:  This part takes the list of remote dofs (some of which are not even on
    //         the border) which talk to our local dofs and updates the appropriate 
    //         local rows to contain the proper global ids of the remote dofs.  
    //         These values were not present before this section of code.  
    //         This is *not* a mapping operation.
    //
    for(unsigned int nn = 0;nn < info.nborder;nn++){
      std::vector<Mesh::IndexType>::iterator bni = borders[nn].nrecv.begin();
      std::vector<Mesh::IndexType>::iterator Api   = borders[nn].data.RecvAp.begin();
      // (1) For every node that I own on this border...
      while(bni != borders[nn].nrecv.end()){
	Mesh::IndexType begindex = *Api++;
	std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[*bni-1].begin();
	// (2) Go through it's local dof id's and..
	Mesh::IndexType nvals = *Api - begindex;
	while(ndi != NodalDofs[*bni-1].end()){
	  // (3) add the remote dof id's to the corresponding local row
	  for(unsigned int i = 0;i < nvals;i++)
	    symbstiff[*ndi-1].push_back(borders[nn].data.RecvAi[begindex+i]);
    	  ndi++;
    	}
	bni++;
      }
    }
    global.FunctionExit("K update");
    if(LogOut && debug && array_checking)
      *LogOut << "K after update:" << std::endl << symbstiff
    	      << std::endl;
    // Finally, sort and unique - oh yeah, yuck.  this will be a 
    // bottleneck. Yup, its slow alrite.
    global.FunctionEntry("K sort");
    Mesh::Connectivity::iterator ssIt = symbstiff.begin();
    while(ssIt != symbstiff.end()){
      std::sort(ssIt->begin(),ssIt->end());
      std::vector<Mesh::IndexType> tmp(ssIt->begin(),std::unique(ssIt->begin(),ssIt->end()));
      ssIt->swap(tmp);
      ssIt++;
    }
    global.FunctionExit("K sort");
    if(LogOut && debug && array_checking)
      *LogOut << "Sorted K:" << std::endl << symbstiff
    	      << std::endl;
	
    comm.Barrier();
    global.FunctionEntry("Create Stiffness");
    symbstiff.Truncate();
    symbstiff.ShrinkWrap();
    FEM::DummyStiffness<double,Mesh::IndexType,Mesh::Connectivity,std::vector<Mesh::IndexType> > k;
    k._dofs = &symbstiff;
    k.rank = rank;
    // Commented out for symbolic assembly
    k._data.resize(1000,0.0);
    k._sizes.resize(nglobdof+1,0);
    Mesh::Connectivity::iterator dci = symbstiff.begin();
    Mesh::IndexType indcnt = 1;
    while(indcnt <= nglobdof){
      k._sizes[indcnt] = dci->size()+k._sizes[indcnt-1];
      dci++;
      indcnt++;
    }
    // ngdofs_p = # global dofs on each processor

    Mesh::IndexType nnz = GetTotalSize(symbstiff);
    
    k._ndof = nglobdof;
    std::vector<Mesh::IndexType> ApLoc(1,k._sizes[nglobdof]);
    std::vector<Mesh::IndexType> ApTot(1,0);

    comm.Reduce(ApLoc,ApTot,IndexDT,Comm::SUMOP,0);
    if(LogOut && verblevel){
      std::vector<Mesh::IndexType>::iterator si = k._sizes.begin();
      std::vector<Mesh::IndexType>::iterator si2 = k._sizes.begin();
      Mesh::IndexType myav = 0;
      //	Mesh::IndexType theav = 0;
      si++;
      while(si != k._sizes.end())
	myav += (*si++ - *si2++);
      myav /= (k._sizes.size()-1);
      *LogOut << "The average NZ per stiffness row is: " << myav << std::endl;
    }

    if(LogOut && debug)
      *LogOut << "Siffness storage(bytes): " 
	      << nnz*datasize << std::endl;
    Mesh::IndexType totsize = 0;
    MPI_Reduce(&k._sizes[nglobdof],&totsize,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
    if(StdOut && verblevel) 
      *StdOut << "Total number of local Stiffness entries: " << nnz 
	      << ":" << k._sizes[nglobdof] << ":" << totsize << std::endl
	      << "Total NZ over all procs: " << ApTot[0] << std::endl;

    if(do_dump_sparse){
      std::ostringstream Ostr;
      std::ostringstream ZeroStr;
      unsigned int pwr10 = 10;
      while(nproc/pwr10 && !(rank/pwr10)){
	ZeroStr << "0";
	pwr10 *= 10;
      }
      Ostr << infiles[0] << "_sparse_" 
	   << ZeroStr.str() << rank+1;
      std::ofstream Ouf;
      Ouf.open(Ostr.str().c_str());
      Ouf << ndof_global_total << std::endl;
      Ouf << k._ndof << std::endl;
      k.Dump(Ouf);
      Ouf.close();
    }

    global.FunctionExit("Create Stiffness");
    global.FunctionExit("Build Sparsity");



    std::vector<Mesh::IndexType> Order;

    //    bool do_renumbering = true;
    if(do_renumbering){
      global.FunctionEntry("Renumbering");
      std::vector<idxtype> adj;
      std::vector<idxtype> xadj;
      std::vector<idxtype> vtxdist;
      vtxdist.resize(ngdofs_p.size()+1);
      std::vector<Mesh::IndexType>::iterator ngdi = ngdofs_p.begin();
      std::vector<idxtype>::iterator vtdi = vtxdist.begin();
      *vtdi = 0;
      while(ngdi != ngdofs_p.end()){
	idxtype nn = static_cast<idxtype>(*vtdi++ + *ngdi++);
	*vtdi = nn;
      }
      symbstiff.Truncate();
      symbstiff.graph_mode(info.doffset);
      if(array_checking && LogOut){
	*LogOut << "Graph K:" << std::endl << symbstiff
		<< std::endl;
      }
      // this function does the 1-based to 0-based mapping
      global.FunctionEntry("Container2CSR");
      MultiContainer2CSR< std::vector< std::vector<Mesh::IndexType> >, 
	std::vector<Mesh::IndexType>,
	std::vector<idxtype>,idxtype>
	(xadj,adj,symbstiff);
      symbstiff.matrix_mode(info.doffset);
      global.FunctionExit("Container2CSR");
      if(array_checking && LogOut){
	*LogOut << "Matrix K:" << std::endl << symbstiff
		<< std::endl;
      }
      std::vector<idxtype> order(ngdofs_p[rank],0);
      std::vector<idxtype> sizes(ngdofs_p.size()*2,0);
      int communicator = MPI_COMM_WORLD;
      int numflag = 0;
      int options[5];
      options[0] = 1;
      options[1] = 3;
      options[2] = 15;
      if(array_checking){
	*LogOut << "Vtxdist:" << std::endl;
	DumpContents(*LogOut,vtxdist," ");
	*LogOut << std::endl << "Xadj:" << std::endl;
	DumpContents(*LogOut,xadj," ");
	*LogOut << std::endl << "Adj:" << std::endl;
	DumpContents(*LogOut,adj," ");
	*LogOut << std::endl << "Ngdofs_p:" << std::endl;
	DumpContents(*LogOut,ngdofs_p," ");
	*LogOut << std::endl;
      }
      global.FunctionEntry("Metis-NodeID");
#ifdef PMETIS
      ParMETIS_V3_NodeND(&vtxdist[0],&xadj[0],&adj[0],&numflag,
			 options,&order[0],&sizes[0],&communicator);
#endif
      global.FunctionExit("Metis-NodeID");
      if(array_checking){
	*LogOut << "Order:" << std::endl;
	DumpContents(*LogOut,order," ");
	*LogOut << std::endl << "Sizes:" << std::endl;
	DumpContents(*LogOut,sizes," ");
	*LogOut << std::endl;
      }
      global.FunctionEntry("OrderSort");
      std::vector<idxtype> sorder(order);
      Order.resize(order.size());
      std::vector<idxtype>::iterator oi = order.begin();
      std::vector<Mesh::IndexType>::iterator Oi = Order.begin();
      while(oi != order.end())
	*Oi++ = *oi++;
      order.resize(0);
      std::sort(sorder.begin(),sorder.end());
      global.FunctionExit("OrderSort");
      std::vector<idxtype>::iterator soi = sorder.begin();
      *LogOut << "Reordered Min: " << *sorder.begin() << std::endl
	      << "Reordered Max: " << *--sorder.end() << std::endl;
      unsigned int ttcount = 0;
      idxtype max_contig = 0;
      idxtype min_contig = 0;
      //      idxtype contig_median = 0;
      idxtype contig_min = 0;
      idxtype contig_max = 0;
      int ncontig = 0;
      struct {
	int nbin;
	int rank;
      } binrank_p[nproc], binrank[nproc];
      for(unsigned int iii = 0;iii < nproc;iii++){
	binrank[iii].nbin = 0;
	binrank[iii].rank = rank;
	binrank_p[iii].rank = 0;
	binrank_p[iii].nbin = 0;
      }
      std::vector<unsigned int> nbin(nproc,0);
      int binsize = ndof_global_total/nproc;
      while(soi != sorder.end()){
	idxtype first  = *soi++;
	idxtype second = first;
	bool unbinned = true;
	unsigned int aaa = 0;
	int binmin = 0;
	while(aaa < nproc && unbinned){
	  if(first < (binmin + binsize)){
	    unbinned = false;
	    binrank[aaa].nbin++;
	  }
	  binmin += binsize;
	  aaa++;
	}
	if(soi != sorder.end()){
	  second = *soi;
	  if(second != first+1){
	    if(ncontig > 1000 && debug)
	      *LogOut << "Large contiguous block(" << ncontig << "): (" 
		      << min_contig << "," << first << ")" << std::endl;
	    if(ncontig > max_contig){
	      max_contig = ncontig;
	      contig_max = first;
	      contig_min = min_contig;
	    }
	    min_contig = second;
	    ncontig = 0;
	    if(array_checking){
	      *LogOut << "NonSequential(" << ttcount << "): " << first << ":" 
		      << second << std::endl;
	    }
	  }
	  else
	    ncontig++;
	}
	ttcount++;
      }
      if(debug){
	*LogOut << "Maximum contiguous block(" << max_contig << "): (" << contig_min
		<< "," << contig_max << ")" << std::endl;
      }
      int binmin = 0;
      int maxbin = 0;
      int solve_rank = 0;
      for(unsigned int iii = 0;iii < nproc;iii++){
	if(debug)
	  *LogOut << "Bin(" << binmin << "," << binmin+binsize << "): " << binrank[iii].nbin 
		  << std::endl;
	if(binrank[iii].nbin > maxbin){
	  maxbin = binrank[iii].nbin;
	  solve_rank = iii;
	}
	binmin+=binsize;
      }
      *LogOut << "My Preferred Solve Rank = " << solve_rank << std::endl;
      MPI_Allreduce(binrank,binrank_p,nproc,MPI_2INT,MPI_MAXLOC,MPI_COMM_WORLD);
      std::vector<unsigned int> solve_ranks(nproc,0);
      std::vector<unsigned int> srank_to_rank(nproc,0);
      std::vector<bool> solve_rank_covered(nproc,false);
      std::vector<unsigned int> rank_conflict;
      comm.Barrier();
      for(unsigned int iii = 0;iii < nproc;iii++){
	solve_ranks[binrank_p[iii].rank] = iii;
	srank_to_rank[iii] = binrank_p[iii].rank;
      }      
      //      MPI_Comm solvecomm = MPI_COMM_NULL;
      Comm::CommunicatorObject solvercomm;
      commstat = comm.Split(0,solve_ranks[rank],solvercomm);
      //      int mpistat = MPI_Comm_split(MPI_COMM_WORLD,0,binrank_p[rank].rank,&solvecomm);
      if(commstat){
	*LogOut << "MPI Error: " << commstat << std::endl;
	return(1);
      }
      //      MPI_Comm_rank(solvecomm,&solve_rank);
      //      *LogOut << "My Solve Rank(1) = " << solve_rank << std::endl;
      solve_rank = solvercomm.Rank();
      *LogOut << "My Solve Rank = " << solve_rank << std::endl;
     

      unsigned int ncovered = 0;
      for(unsigned int iii = 0;iii < nproc;iii++){
	*LogOut << "SolveRanks[" << iii << "] == " << solve_ranks[iii] 
		<< std::endl
		<< "SRank[" << iii << "] == " << srank_to_rank[iii] << std::endl;
	if(rank != iii){
	  if(solve_ranks[iii] == solve_ranks[rank])
	    rank_conflict.push_back(iii);
	  else{
	    solve_rank_covered[solve_ranks[iii]] = true;
	    ncovered++;
	  }
	}
	else{
	  solve_rank_covered[solve_ranks[iii]] = true;
	  ncovered++;
	}
      }
      if(!rank_conflict.empty()){
	comm.SetExit(1);
	*LogOut << "Conflicting for solve_rank(" << solve_rank 
		<< "): ";
	DumpContents(*LogOut,rank_conflict," ");
	*LogOut << std::endl;
      }
      if(ncovered != nproc){
	comm.SetExit(1);
	*LogOut << "Uncovered bins: ";
	unsigned int bincount = 0;
	std::vector<bool>::iterator rci = solve_rank_covered.begin();
	while(rci != solve_rank_covered.end()){
	  if(*rci++ == false)
	    *LogOut << bincount << " ";
	  bincount++;
	}
	*LogOut << std::endl;
      }
      if(comm.Check())
	return(1);
      int binxtra = ndof_global_total%nproc;
      unsigned int mynrows = binsize + (solve_rank < binxtra ? 1 : 0);
      std::vector<unsigned int> nrows_p(nproc,0);
      nrows_p[solve_rank] = mynrows;
      commstat = solvercomm.AllGather(mynrows,nrows_p);
      solvercomm.Barrier();
      int solve_doffset = 0;
      for(int iii = 0;iii < solve_rank;iii++)
	solve_doffset += nrows_p[iii];
      int my_first_row = solve_doffset;
      int my_last_row  = my_first_row + mynrows - 1;
      if(debug)
	*LogOut << "My Rows (" << mynrows << ") = (" << my_first_row 
		<< "," << my_last_row << ")" << std::endl;
      Mesh::Connectivity RemoteRows;
      RemoteRows.resize(nproc);
      for(unsigned int iii = 0;iii < nglobdof;iii++){
	unsigned int row_index = Order[iii];
	unsigned int solve_procid = 0;
	unsigned int findrow = nrows_p[0];
	while(findrow < row_index)
	  findrow += nrows_p[++solve_procid];
	RemoteRows[solve_procid].push_back(row_index);
      }
      if(array_checking)
	*LogOut << "RemoteRows:" << std::endl
		<< RemoteRows << std::endl;
      if(debug){
	*LogOut << "Row distribution: " << std::endl;
	for(unsigned int iii = 0;iii < nproc;iii++)
	  *LogOut << "Rows on Solve Rank " << iii << ": " 
		  << RemoteRows[iii].size() << std::endl;
      }
      unsigned int nlocal_rows = RemoteRows[solve_rank].size();
      if(debug)
	*LogOut << "Number of solved Rows I own: " << nlocal_rows << std::endl
		<< "Number of Rows I assemble: " << nglobdof << std::endl
		<< "Number of Rows I solve: " << mynrows << std::endl;
      int nremote_rows = 0;
      for(int iii = 0;iii < (int)nproc;iii++)
	if(iii != solve_rank)
	  nremote_rows += RemoteRows[iii].size();
      int nrecvd_rows = mynrows - nlocal_rows;
      if(debug)
	*LogOut << "Number of Remote Rows I assemble (Sent Rows): " << nremote_rows << std::endl
		<< "Number of my Rows assembled remotely (Recvd Rows): " << nrecvd_rows << std::endl;
      // Now need to communicate and discover the mapping of old remote dof numbers
      // to new.   There are two sets of dofs that we need: 
      // 1. The new numbers for the dofs on remotely owned nodes
      // 2. The new numbers for the remote dofs which we receive
      // 
      // Once we have this information, we can update our version of K so
      // that it's sparsity pattern is consistent with the renumbering.
      GlobalDofReMapExchange(econ,NodalDofs,ElementDofs,RemoteNodalDofs,
			     borders,info,Order,comm,LogOut,array_checking);
     
      if(true){
	global.FunctionEntry("Build Sparsity");
	// Now to actually build the sparsity pattern, LD[GD]
	// We have E[LED], N[LND], Nr[GD], LD:GD
	// We'll start by building a local version of K, LD[LD],
	// and using a straightforward mapping LD:GD to obtain LD[GD].
	// Then we'll add any needed global dofs coming from remote
	// processors.
	//
	// 1. Build local K, LD[GD]
	//    - Need to somehow merge E[LED] and Nl[LND], since
	//      we always assemble element-by-element, we should
	//      rack up the LD on the elements in the order of
	//      all nodal dofs first, then element dofs. 
	//
	global.FunctionEntry("Create E[LD]");
	Mesh::Connectivity E_LD;
	E_LD.resize(number_of_elements);
	for(Mesh::IndexType en = 0;en < number_of_elements;en++){
	  std::vector<Mesh::IndexType>::iterator eni = econ[en].begin();
	  while(eni != econ[en].end()){ // for every node in this element
	    Mesh::IndexType node_id = *eni++;
	    std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[node_id-1].begin();
	    while(ndi != NodalDofs[node_id-1].end())
	      E_LD[en].push_back(*ndi++); // add it's dofs to the list of dofs for this element
	    
	  }
	  std::vector<Mesh::IndexType>::iterator edi = ElementDofs[en].begin();
	  while(edi != ElementDofs[en].end())
	    E_LD[en].push_back(*edi++); // add the element's dofs as well
	}
	global.FunctionEntry("E[LD] post");
	E_LD.ShrinkWrap();
	E_LD.Sync();
	E_LD.SyncSizes();
	// We need the sizes later
	std::vector<Mesh::IndexType> element_ndofs(number_of_elements,0);
	for(Mesh::IndexType ei = 0;ei < number_of_elements;ei++)
	  element_ndofs[ei] = E_LD.Esize(ei+1);
	global.FunctionExit("E[LD] post");
	global.FunctionExit("Create E[LD]");
	//
	// Now E_LD is a list of all dofs for every element.  It is an important point
	// that LD are the LD (i.e. the local dofs), because they run from 1-nlocal_dofs.
	// Since both the indices of E (1-number_of_elements) and dofs both are sequential
	// lists, we can do some nifty array manipulations to obtain the information we need.
	// First, we invert E[LD] to obtain LD[E]:
	Mesh::Connectivity LD_E;
	global.FunctionEntry("Invert E[LD]");
	E_LD.Inverse(LD_E,ndoftot);
	global.FunctionExit("Invert E[LD]");
	//
	// Now, there is one slight problem with our version of LD[E], in that it includes
	// LD which are actually just local numbers for remote dofs.  We do not need the
	// stiffness rows for these.  In fact, we already have created the stiffness blocks 
	// for these rows when we created the communication buffers.  We know those dofs are 
	// all dofs with id's > nglobdof.  If we don't take care of the problem now, we pay
	// a hefty memory and processing fee in the following steps.
	LD_E.Sync();
	assert(LD_E.Nelem() == ndoftot);
	for(Mesh::IndexType di = nglobdof;di < ndoftot;di++)
	  LD_E[di].resize(0);
	LD_E.ShrinkWrap();
	LD_E.SyncSizes();
	
	// LD[E]%E[LD] = LD[LD] := Stiffness
	symbstiff.destroy();
	global.FunctionEntry("Sparsity Pattern");
	LD_E.GetNeighborhood(symbstiff,E_LD,false,false);
	global.FunctionExit("Sparsity Pattern");
	
	comm.Barrier();
	LD_E.destroy();
	E_LD.destroy();
	
	if(LogOut && debug && array_checking)
	  *LogOut << "All Local K:" << std::endl << symbstiff
		  << std::endl;
	
	
	// Now we need to do two things, update our "K" with the dofs which
	// are coming from remote procs and also map all the dof id's to 
	// global.  Then we're done building the parallel version of K.
	//
	// First, we need a map from local dof id's for the remote dofs 
	// on the partition boundary to the global id's that have been
	// told to us by the remote processor.  
	//
	// *Note:  This code updates the local dof rows which already
	// contain local dof numbers for remotely owned dofs.  These 
	// remote dofs are just the dofs of the remotely owned border
	// nodes.  This is just a mapping operation.
	//
	global.FunctionEntry("K update");
	if(true){
	  std::vector<Mesh::IndexType> local_dof_to_global;
	  local_dof_to_global.resize(ndoftot-nglobdof);
	  for(Mesh::IndexType ni = info.nlocal;ni < number_of_nodes;ni++){
	    std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[ni].begin();
	    std::vector<Mesh::IndexType>::iterator rdi = RemoteNodalDofs[ni-info.nlocal].begin();
	    assert(NodalDofs[ni].size() == RemoteNodalDofs[ni-info.nlocal].size());
	    while(ndi != NodalDofs[ni].end()){
	      local_dof_to_global[*ndi-nglobdof-1] = *rdi;
	      ndi++;
	      rdi++;
	    }
	  }
	  Mesh::Connectivity::iterator ConIt = symbstiff.begin();
	  while(ConIt != symbstiff.end()){
	    std::vector<Mesh::IndexType>::iterator  dofIt = ConIt->begin();
	    Mesh::IndexType ind = 0;
	    while(dofIt != ConIt->end()){
	      Mesh::IndexType dof_id = *dofIt++;
	      if(dof_id <= nglobdof)
		(*ConIt)[ind++] = Order[dof_id-1]+1;
	      else
		(*ConIt)[ind++] = local_dof_to_global[dof_id - nglobdof - 1];
	    }
	    ConIt++;
	  }
	}
	//
	// Now update with the remote dof info.
	// *Note:  This part takes the list of remote dofs (some of which are not even on
	//         the border) which talk to our local dofs and updates the appropriate 
	//         local rows to contain the proper global ids of the remote dofs.  
	//         These values were not present before this section of code.  
	//         This is *not* a mapping operation.
	//
	for(unsigned int nn = 0;nn < info.nborder;nn++){
	  std::vector<Mesh::IndexType>::iterator bni = borders[nn].nrecv.begin();
	  std::vector<Mesh::IndexType>::iterator Api   = borders[nn].data.RecvAp.begin();
	  // (1) For every node that I own on this border...
	  while(bni != borders[nn].nrecv.end()){
	    Mesh::IndexType begindex = *Api++;
	    std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[*bni-1].begin();
	    // (2) Go through it's local dof id's and..
	    Mesh::IndexType nvals = *Api - begindex;
	    while(ndi != NodalDofs[*bni-1].end()){
	      // (3) add the remote dof id's to the corresponding local row
	      for(unsigned int i = 0;i < nvals;i++)
		symbstiff[*ndi-1].push_back(borders[nn].data.RecvAi[begindex+i]);
	      ndi++;
	    }
	    bni++;
	  }
	}
	global.FunctionExit("K update");
	if(LogOut && debug && array_checking)
	  *LogOut << "K after update:" << std::endl << symbstiff
		  << std::endl;
	// Finally, sort and unique - oh yeah, yuck.  this will be a 
	// bottleneck. Yup, its slow alrite.
	global.FunctionEntry("K sort");
	Mesh::Connectivity::iterator ssIt = symbstiff.begin();
	while(ssIt != symbstiff.end()){
	  std::sort(ssIt->begin(),ssIt->end());
	  std::vector<Mesh::IndexType> tmp(ssIt->begin(),std::unique(ssIt->begin(),ssIt->end()));
	  ssIt->swap(tmp);
	  ssIt++;
	}
	global.FunctionExit("K sort");
	if(LogOut && debug && array_checking)
	  *LogOut << "Sorted K:" << std::endl << symbstiff
		  << std::endl;
	
	comm.Barrier();
	global.FunctionEntry("Create Stiffness");
	symbstiff.Truncate();
	symbstiff.ShrinkWrap();
	//	FEM::DummyStiffness<double,Mesh::IndexType,Mesh::Connectivity,std::vector<Mesh::IndexType> > k;
	k._dofs = &symbstiff;
	k.rank = rank;
	// Commented out for symbolic assembly
	k._data.resize(1000,0.0);
	k._sizes.resize(nglobdof+1,0);
	Mesh::Connectivity::iterator dci = symbstiff.begin();
	Mesh::IndexType indcnt = 1;
	while(indcnt <= nglobdof){
	  k._sizes[indcnt] = dci->size()+k._sizes[indcnt-1];
	  dci++;
	  indcnt++;
	}
      } 
      Mesh::IndexType nnz = GetTotalSize(symbstiff);
      
      k._ndof = nglobdof;
      std::vector<Mesh::IndexType> ApLoc(1,k._sizes[nglobdof]);
      std::vector<Mesh::IndexType> ApTot(1,0);
      // comm.Reduce does not work - since we convert datatypes :P      
      //      comm.Reduce(ApLoc,ApTot,MPI_SUM,0);
      if(LogOut && verblevel){
	std::vector<Mesh::IndexType>::iterator si = k._sizes.begin();
	std::vector<Mesh::IndexType>::iterator si2 = k._sizes.begin();
	Mesh::IndexType myav = 0;
	//	Mesh::IndexType theav = 0;
	si++;
	while(si != k._sizes.end())
	  myav += (*si++ - *si2++);
	myav /= (k._sizes.size()-1);
	*LogOut << "The average NZ per stiffness row is: " << myav << std::endl;
      }
      
      if(LogOut && debug)
	*LogOut << "Stiffness storage(bytes): " 
		<< nnz*datasize << std::endl;
      Mesh::IndexType totsize = 0;
      MPI_Reduce(&k._sizes[nglobdof],&totsize,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
      if(debug) 
	*LogOut << "Total number of local Stiffness entries: " 
		<< k._sizes[nglobdof] << std::endl
		<< "Total NZ over all procs: " << totsize << std::endl;
      if(do_dump_sparse){
	std::ostringstream Ostr;
	std::ostringstream ZeroStr;
	unsigned int pwr10 = 10;
	while(nproc/pwr10 && !(rank/pwr10)){
	  ZeroStr << "0";
	  pwr10 *= 10;
	}
	Ostr << infiles[0] << "_sparse_rnm_" 
	     << ZeroStr.str() << rank+1;
	std::ofstream Ouf;
	Ouf.open(Ostr.str().c_str());
	Ouf << ndof_global_total << std::endl;
	Ouf << k._ndof << std::endl;
	DumpContents(Ouf,Order," ");
	Ouf << std::endl;
	k.Dump(Ouf);
	Ouf.close();
      }
      global.FunctionExit("Create Stiffness");
      global.FunctionExit("Build Sparsity");
      std::vector<bool> element_processed(number_of_elements,false);
      std::vector<bool> is_border_element(number_of_elements,false);
      for(unsigned int iii = 0;iii < info.nlocal;iii++){
	std::vector<Mesh::IndexType>::iterator ndi = NodalDofs[iii].begin();
	bool node_processed = false;
	while(ndi != NodalDofs[iii].end() && !node_processed){
	  if(Order[*ndi] <= (unsigned int)my_first_row ||
	     Order[*ndi] >= (unsigned int)my_last_row){
	    std::vector<Mesh::IndexType>::iterator bei = dc[iii].begin();
	    while(bei != dc[iii].end())
	      is_border_element[*bei++] = true;
	    node_processed = true;
	  }
	  ndi++;
	}
      }
      for(unsigned int iii = info.nlocal;iii < number_of_nodes;iii++){
	std::vector<Mesh::IndexType>::iterator bei = dc[iii].begin();
	while(bei != dc[iii].end())
	  is_border_element[*bei++] = true;
      }
      for(unsigned int iii = 0;iii < number_of_elements;iii++){
	std::vector<Mesh::IndexType>::iterator edi = ElementDofs[iii].begin();
	bool element_processed = false;
	while(edi != ElementDofs[iii].end() && !element_processed){
	  if(Order[*edi] <= (unsigned int)my_first_row || 
	     Order[*edi] >= (unsigned int)my_last_row){
	    is_border_element[iii] = true;
	    element_processed = true;
	  }
	  edi++;
	}
      }
      std::vector<bool>::iterator bei = is_border_element.begin();
      unsigned int nborder_elements = 0;
      while(bei != is_border_element.end()){
	if(*bei++)
	  nborder_elements++;
      }
      if(debug)
	*LogOut << "Border Elements: " << nborder_elements 
		<< "/" << number_of_elements << std::endl;
      global.FunctionExit("Renumbering");      
    }
    
    // Border Node Mapping
    // In order to to border assemblies, I need to have this map:
    // BN:<C,id>, (i.e. For each Boundary Node, it's color and id
    // where id is the index of it's entry in the nsend array.
    global.FunctionEntry("Border Node Mapping");
    std::vector<std::pair<unsigned int,unsigned int> > boundary_node_map(info.nnodes-info.nlocal);
    for(unsigned int i = 0;i < info.nborder;i++){
      std::vector<Mesh::IndexType>::iterator sni = borders[i].nsend.begin();
      unsigned int index = 0;
      while(sni != borders[i].nsend.end()){
	boundary_node_map[*sni - info.nlocal - 1] = std::make_pair(i+1,index+1);
	sni++;
	index++;
      }
    }
    global.FunctionExit("Border Node Mapping");
    std::list<Mesh::IndexType> nonborder_queue;
    if(true){
      // Re-define border elements to be only those touching remote nodes
      global.FunctionEntry("Find Border Elements");
      is_border_element.clear();
      is_border_element.resize(0);
      is_border_element.swap(is_border_element);
      element_processed.clear();
      element_processed.resize(0);
      element_processed.resize(number_of_elements,false);
      border_elements.clear();
      border_elements.resize(0);
      for(Mesh::IndexType nn = info.nlocal;nn < number_of_nodes;nn++){
	std::vector<Mesh::IndexType>::iterator ei = dc[nn].begin();
	while(ei != dc[nn].end()){
	  if(!element_processed[*ei-1]){
	    element_processed[*ei-1] = true;
	    border_elements.push_back(*ei);
	  }
	  ei++;
	}
      }
      border_elements.sort();
      std::list<Mesh::IndexType>::iterator ei = queue.begin();
      while(ei != queue.end()){
	if(!element_processed[*ei-1])
	  nonborder_queue.push_back(*ei);
	ei++;
      }
      if(LogOut && verblevel){
	*LogOut << "Found " << border_elements.size() << " border elements." << std::endl
		<< "Found " << nonborder_queue.size() << " non border elements."  
		<< std::endl;
	if(array_checking){
	  *LogOut << "Border Elements:" << std::endl;
	  DumpContents(*LogOut,border_elements," ");
	  *LogOut << std::endl << "NonBorder Elements:" << std::endl;
	  DumpContents(*LogOut,nonborder_queue," ");
	  *LogOut << std::endl;
	}
      }
      global.FunctionExit("Find Border Elements");
    }    
    element_processed.clear();
    element_processed.resize(0);
    element_processed.resize(number_of_elements,false);
    global.FunctionExit("Assembly Prep");
    if(false){
      global.FunctionEntry("Search Test");
      Mesh::Connectivity::iterator ci = symbstiff.begin();
      while(ci != symbstiff.end()){
	std::vector<Mesh::IndexType>::iterator ssrbegin = ci->begin();
	std::vector<Mesh::IndexType>::iterator ssrend   = ci->end();
	std::vector<Mesh::IndexType>::iterator ssi      = ci->begin();
	std::vector<Mesh::IndexType>::iterator curpos   = ci->begin();
	std::vector<Mesh::IndexType>::iterator last_one = ci->begin();
	while(ssi != ci->end()){
	  last_one = curpos;
	  curpos = std::lower_bound(ssrbegin,ssrend,*ssi++);
	}
	ci++;
      }
      global.FunctionExit("Search Test");
    }
    //    return(0);
    comm.Barrier();
    // Now go for a true assembly
    // First, work up some fake data, howabout all 1.0's.
    std::vector<double> dofdat(100000,1.0);

    global.FunctionEntry("Assembly Processing");
    Mesh::IndexType nsearches = 0;
    nsearches = 0;
    
    // Post receives
    global.FunctionEntry("Posting receives");
    comm.RecvAll();
    global.FunctionExit("Posting receives");
    comm.Barrier();
    
    Mesh::IndexType bordersearch = 0;
    if(true){
      global.FunctionEntry("Send Buffer Assembly");
      std::list<Mesh::IndexType>::iterator bei = border_elements.begin();
      while(bei != border_elements.end()){
	Mesh::IndexType eindex = *bei++ - 1;
	element_processed[eindex]=true;
	Mesh::IndexType endof  = element_ndofs[eindex]; 
	if(!((eindex+1) > 0 && (eindex+1) <= ElementDofs.Nelem())){
	  *ErrOut << "ERROR: ElementDofs[" << eindex << "].size() = " << ElementDofs[eindex].size()
		  << " and ElementDofs.size() = " << ElementDofs.Nelem() << std::endl;
	  return(1);
	}
	Mesh::IndexType nedof  = ElementDofs.Esize(eindex+1); 
	Mesh::IndexType search = FEM::AssembleBorderElementII(econ[eindex],ElementDofs[eindex],endof,nedof,
							      NodalDofs,RemoteNodalDofs,k,borders,
							      boundary_node_map,dofdat,info,LogOut,debug);
	bordersearch+=search;
	
      }
      global.FunctionExit("Send Buffer Assembly");
    }
    comm.Barrier();
    global.FunctionEntry("Posting Sends");
    comm.SendAll();
    global.FunctionExit("Posting Sends");
    comm.Barrier();
    if(false){
      global.FunctionEntry("Assembly Loop");
      // loop over elements
      //      for(Mesh::IndexType nnn = 0;nnn < number_of_elements;nnn++){
      std::list<Mesh::IndexType>::iterator ei = nonborder_queue.begin();
      while(ei != nonborder_queue.end()){
	//      std::list<Mesh::IndexType>::iterator ei = queue.begin();
	//      while(ei != queue.end()){
	Mesh::IndexType nnn = *ei++ -1;
	Mesh::IndexType search = 0;
	//	if(!element_processed[nnn]){
	  Mesh::IndexType endof = element_ndofs[nnn];
	  Mesh::IndexType nedof = ElementDofs.Esize(nnn+1);
	  search = FEM::AssembleLocalElementII(econ[nnn],ElementDofs[nnn],endof,nedof,
					       NodalDofs,k,dofdat,info);
	  nsearches+=search;
	  //	}
      }
      global.FunctionExit("Assembly Loop");
    }
    Mesh::IndexType fast_searches = 0;
    if(true){
      global.FunctionEntry("Fast Assembly");
      fast_searches = FEM::FastAssembleLocalElements(nonborder_queue,econ,k,
						     NodalDofs,ElementDofs,element_ndofs,info);
      global.FunctionExit("Fast Assembly");
    }
    comm.Barrier();
    global.FunctionEntry("Wait for messages");
    comm.WaitAll();
    global.FunctionExit("Wait for messages");
    comm.Barrier();
    // Need to add code to assemble the recv'd dofs to K
    global.FunctionEntry("Recv Buffer Assembly");
    Mesh::IndexType postsearch = 0;
    if(true){
      for(unsigned int i = 0;i < info.nborder;i++){
	int search = RecvBufAssembly(borders[i],NodalDofs,k,dofdat);
	postsearch+=search;
      }
    }
    global.FunctionExit("Recv Buffer Assembly");
    global.FunctionExit("Assembly Processing");
    if(StdOut && verblevel) 
      *StdOut << "Slow Ass searches = " << nsearches    << std::endl
	      << "Border searches   = " << bordersearch << std::endl
	      << "Recv Buf searches = " << postsearch   << std::endl
	      << "Fast Ass searches = " << fast_searches << std::endl;
    comm.ClearRequests();
    comm.Barrier();
    global.FunctionExit("Assembly Function");
  }
  comm.Barrier();




  //  bool do_partitioning = true;
  if(do_partitioning){
    if(nproc > 1){
      if(ErrOut)
	*ErrOut << "Error: partitioning is for serial runs." << std::endl;
      comm.SetExit(1);
    }
    if(comm.Check())
      return(1);
    global.FunctionEntry("Partitioning");
    std::ifstream Inf;
    std::string MeshName(infiles[0]);
    std::string MeshBaseName(MeshName.substr(0,MeshName.find_last_of(".")));
    Inf.open(MeshName.c_str());
    if(!Inf){
      if(ErrOut) 
	*ErrOut << "Failed to open " << MeshName << std::endl;
      return(1);
    }
    global.FunctionEntry("Read Mesh");
    Mesh::IndexType number_of_nodes = 0;
    if(true){
      Mesh::NodalCoordinates nc;
      Inf >> nc;
      number_of_nodes = nc.Size();
      nc.destroy();
    }
    if(StdOut && verblevel) 
      *StdOut << "Number of nodes: " << number_of_nodes << std::endl;
    Mesh::Connectivity econ;
    Inf >> econ;
    global.FunctionExit("Read Mesh");
    Inf.close();
    econ.ShrinkWrap();
    Mesh::IndexType number_of_elements = econ.Nelem();
    if(StdOut && verblevel) 
      *StdOut << "Number of elements: " << number_of_elements
	      << std::endl;
    Mesh::IndexType N = npartitions;
    econ.SyncSizes();
    // re-orient the elements if needed
    if(do_reorient){ // don't do this for 2d geometry (yet)
      Mesh::NodalCoordinates nc;
      Inf.open(MeshName.c_str());
      if(!Inf){
	if(ErrOut) 
	  *ErrOut << "Failed to open " << MeshName << std::endl;
	return(1);
      }
      global.FunctionEntry("ReadCoords");
      Inf >> nc;
      Inf.close();
      global.FunctionExit("ReadCoords");
      global.FunctionEntry("Element Reorientation");
      for(Mesh::IndexType element_being_processed = 1;
	  element_being_processed <= number_of_elements;
	  element_being_processed++)
	{
	  Mesh::IndexType index = element_being_processed - 1;
	  Mesh::IndexType size_of_element = econ.Esize(element_being_processed);
	  Mesh::GenericElement ge(size_of_element);
	  if(ge.Inverted(econ[index],nc))
	    ge.ReOrient(econ[index]);
	}
      global.FunctionExit("Element Reorientation");
    }
    //    Mesh::IndexType N2 = 0;
    //    if(infiles.size() > 1){
    //      std::istringstream Istr(infiles[1]);
    //      Istr >> N;
    //    }
    //    if(infiles.size() > 2){
    //      std::istringstream Istr(infiles[2]);
    //      Istr >> N2;
    //    }
    Mesh::Connectivity dc;
    global.FunctionEntry("DualCon");
    econ.Inverse(dc,number_of_nodes);
    global.FunctionExit("DualCon");
    dc.ShrinkWrap();
    Mesh::IndexType number_of_faces = 0;
    Mesh::Connectivity F_N;
    Mesh::Connectivity E_F;
    std::vector<Mesh::SymbolicFace> sf;
    econ.SyncSizes();
    global.FunctionEntry("GetFaceConnectivites");
    econ.BuildFaceConnectivity(F_N,E_F,sf,dc);
    global.FunctionExit("GetFaceConnectivites");
    global.FunctionEntry("Shrinking");    
    std::vector<Mesh::SymbolicFace>(sf).swap(sf);
    global.FunctionExit("Shrinking");
    global.FunctionEntry("Creating Graph");
    Mesh::Connectivity nl2(number_of_elements);
    std::vector<Mesh::SymbolicFace>::iterator sfi = sf.begin();
    while(sfi != sf.end()){
      if(sfi->second.first > 0){
	nl2[sfi->first.first-1].push_back(sfi->second.first);
	nl2[sfi->second.first-1].push_back(sfi->first.first);
      }
      sfi++;
    }
    global.FunctionExit("Creating Graph");
   // -------------- Partitioning --------------------
    std::vector<idxtype> adj;
    std::vector<idxtype> xadj;
    global.FunctionEntry("Container2CSR");
    MultiContainer2CSR< std::vector< std::vector<Mesh::IndexType> >, 
      std::vector<Mesh::IndexType>,
      std::vector<idxtype>,idxtype>
      (xadj,adj,nl2);
    global.FunctionExit("Container2CSR");
    nl2.destroy();
    std::vector<idxtype> order(number_of_elements,0);
    std::vector<idxtype> sizes(2,0);
    idxtype wgtflag = 0;
    idxtype numflag = 0;
    idxtype nparts = N;
    idxtype options[5];
    idxtype edgecuts = 0;
    options[0] = 1;
    options[1] = 3;
    options[2] = 15;
    std::vector<idxtype> vtxdist(2,0);
    vtxdist[1] = number_of_elements;
    std::vector<idxtype> part(number_of_elements,0);
    global.FunctionEntry("Metis");
    int comm = MPI_COMM_WORLD;
    //    int ncon = 0;
    //    float stupid = 0;
    //    float ubvecstupid = 0;
    if(array_checking){
      *LogOut << "Vtxdist:" << std::endl;
      DumpContents(*LogOut,vtxdist," ");
      *LogOut << std::endl << "Xadj:" << std::endl;
      DumpContents(*LogOut,xadj," ");
      *LogOut << std::endl << "Adj:" << std::endl;
      DumpContents(*LogOut,adj," ");
      *LogOut << std::endl;
    }
    if(false){
      global.FunctionEntry("Renumbering");    
#ifdef PMETIS  
      ParMETIS_V3_NodeND(&vtxdist[0],&xadj[0],&adj[0],&numflag,
			 options,&order[0],&sizes[0],&comm);
#endif
      global.FunctionExit("Renumbering");      
      if(array_checking){
	*LogOut << "Order:" << std::endl;
	DumpContents(*LogOut,order," ");
	*LogOut << std::endl << "Sizes:" << std::endl;
	DumpContents(*LogOut,sizes," ");
	*LogOut << std::endl;
      }
    }
#ifdef PMETIS
    ParMETIS_V3_PartKway(&vtxdist[0],&xadj[0],&adj[0],NULL,NULL,&wgtflag,
			 &numflag,NULL,&nparts,NULL,NULL,
			 options,&edgecuts,&part[0],&comm);
#else
    METIS_PartGraphKway(&vtxdist[1],&xadj[0],&adj[0],NULL,NULL,&wgtflag,
			&numflag,&nparts,options,&edgecuts,&part[0]);
#endif
    global.FunctionExit("Metis");
    // ------------------------------------------------
    // Take part and create the 1-shifted version into 
    // a multicontainer
    Mesh::Connectivity pcon;
    pcon.resize(number_of_elements);
    for(Mesh::IndexType i = 0;i < number_of_elements;i++)
      pcon[i].push_back(part[i]+1);
    pcon.Sync();
    pcon.SyncSizes();
    //    std::ofstream Ouf;
    //    Ouf.open("part");
    //    Ouf << pcon << std::endl;
    //    Ouf.close();
    if(StdOut && verblevel) 
      *StdOut << "Number of edgecuts: " << edgecuts << std::endl;
    if(false){
      if(StdOut && verblevel) 
	*StdOut << "Partition Table: ";
      DumpContents(*StdOut,part," ");
      if(StdOut && verblevel) 
	*StdOut << std::endl;
    }
    part.resize(0);
    part.swap(part);
    Mesh::IndexType Nparts = nparts;

   // Invert the partition table
    global.FunctionEntry("Dual Pcon");
    Mesh::Connectivity dual_pcon;  // E[C] -> C[E]
    pcon.Inverse(dual_pcon,nparts);
    global.FunctionExit("Dual Pcon");
    dual_pcon.SyncSizes();

    global.FunctionEntry("GetAdj nodecolors");
    Mesh::Connectivity nodecolors; // N[E]%E[C] = N[C]
    dc.GetAdjacent(nodecolors,pcon,nparts,true);
    nodecolors.SyncSizes();
    global.FunctionExit("GetAdj nodecolors");
    
    global.FunctionEntry("Inverse colornodes");
    Mesh::Connectivity colorsofnodes; // N[C] -> C[N]
    nodecolors.Inverse(colorsofnodes,nparts);
    global.FunctionExit("Inverse colornodes");
    colorsofnodes.SyncSizes();

    global.FunctionEntry("ColorNeighbors");
    Mesh::Connectivity colorneighbors; // C[N]%N[C] = C[C]
    colorsofnodes.GetNeighborhood(colorneighbors,nodecolors,true,true);
    global.FunctionExit("ColorNeighbors");
    colorneighbors.SyncSizes();

    // This is enough information to fully describe the partitioning
    // to each partition and form the communication lists.  We can
    // write this to file, or just pass it along in messages to 
    // each partition.

    // For each color, show me the taco
    // produces C[BN] (border nodes for each color)
    Mesh::Connectivity bordercon;
    nodecolors.InverseDegenerate(bordercon,nparts);
    // Lets determine and report the communication imbalance
    bool report_each = false;
    if(report_each)
      if(StdOut && verblevel) 
	*StdOut << "Nshared nodes/partition:" << std::endl;
    double mean_nbnodes = 0;
    Mesh::IndexType max_nbnodes = 0;
    Mesh::IndexType min_nbnodes = number_of_nodes;
    for(Mesh::IndexType iii = 0;iii < Nparts;iii++){
      Mesh::IndexType nbsize = bordercon[iii].size();
      if(report_each && StdOut)
	*StdOut << "Partition " << iii+1 << ": "
		<< nbsize << std::endl;
      mean_nbnodes += static_cast<double>(nbsize);
      if(nbsize > max_nbnodes) max_nbnodes = nbsize;
      if(nbsize < min_nbnodes) min_nbnodes = nbsize;
    }
    mean_nbnodes /= static_cast<double>(nparts);
    if(StdOut && verblevel) 
      *StdOut << "(Max/Min/Mean) shared nodes/partition: (" << max_nbnodes << "/"
	      << min_nbnodes << "/" << mean_nbnodes << ")" << std::endl;

    Mesh::IndexVec bn_index_map;
    std::map<Mesh::IndexType,Mesh::IndexType> bnodemap;
    std::list<Mesh::IndexType> bnodes;
    global.FunctionEntry("Flatten");
    Flatten<Mesh::Connectivity,std::vector<Mesh::IndexType>,std::list<Mesh::IndexType> >
      (bordercon,bnodes);
    global.FunctionExit("Flatten");
    global.FunctionEntry("SortUnique");
    bnodes.sort();
    bnodes.unique();
    global.FunctionExit("SortUnique");
    Mesh::IndexType number_of_border_nodes = bnodes.size();
    bn_index_map.resize(number_of_border_nodes);
    Mesh::IndexType i = 0;
    std::list<Mesh::IndexType>::iterator eli = bnodes.begin();
    while(eli != bnodes.end()){
      bn_index_map[i] = *eli;
      bnodemap.insert(std::make_pair(*eli,i+1));
      eli++;
      i++;
    }
    
    global.FunctionEntry("MapBorderNodes");
    // C[BN] -> BN[C], but need BN to be mapped in 1-NBN range
    Mesh::Connectivity bordercon_mapped;    
    MapElements<Mesh::Connectivity,Mesh::Connectivity,
      std::vector<Mesh::IndexType>,std::map<Mesh::IndexType,Mesh::IndexType> >
      (bordercon,bordercon_mapped,bnodemap);
    global.FunctionExit("MapBorderNodes");
    Mesh::Connectivity bcon_dual;
    bordercon_mapped.Sync();
    bordercon_mapped.Inverse(bcon_dual,number_of_border_nodes);
    // Now we have BN[C], which is a list of colors for each border
    // node.  We need to step through and assign each border node 
    // an actual owner.  This will produce a list
    // of which color (only 1 of them) is the owner of a particular BN.
    // Or, bn[C], where bn are the subset of BN that C owns.
    // We do this in a way so as to equalize the number of communicated
    // nodes on each processors.
    std::srand((unsigned)std::time(0));
    Mesh::Connectivity owner;
    std::vector<Mesh::IndexType> ncomnodes(nparts,0);
    owner.resize(number_of_border_nodes);
    Mesh::Connectivity::iterator bni = bcon_dual.begin();
    Mesh::IndexType ownind = 0;
    while(bni != bcon_dual.end()){
      Mesh::IndexType ncolors = bni->size();
      Mesh::IndexType color = (std::rand()%ncolors)+1;
      for(Mesh::IndexType ci = 0;ci < ncolors;ci++)
	if(ncomnodes[(*bni)[color-1] - 1] > ncomnodes[(*bni)[ci]-1])
	  color = ci+1;
      std::vector<Mesh::IndexType>::iterator bnici = bni->begin();
      for(Mesh::IndexType ci = 1;ci < color;ci++)
	bnici++;
      //      if(bordercon[*bnici-1].size() > static_cast<Mesh::IndexType>(mean_nbnodes))
      //	color = (std::rand()%ncolors)+1;
      owner[ownind++].push_back(*bnici);
      ncomnodes[*bnici-1]++;
      bni++;
    }
    // bn[C] -> C[bn]
    Mesh::Connectivity dual_owner;
    owner.Sync();
    owner.Inverse(dual_owner,nparts);
    mean_nbnodes = 0;
    min_nbnodes = number_of_nodes;
    max_nbnodes = 0;
    for(Mesh::IndexType iii = 0;iii < Nparts;iii++){
      Mesh::IndexType nbsize = dual_owner[iii].size();
      if(report_each)
	if(StdOut && verblevel) 
	  *StdOut << "Partition " << iii+1 << ": "
		  << nbsize << std::endl;
      mean_nbnodes += static_cast<double>(nbsize);
      if(nbsize > max_nbnodes) max_nbnodes = nbsize;
      if(nbsize < min_nbnodes) min_nbnodes = nbsize;
    }
    mean_nbnodes /= static_cast<double>(nparts);
    if(StdOut && verblevel) 
      *StdOut << "(Max/Min/Mean) recvd nodes/partition: (" << max_nbnodes << "/"
	      << min_nbnodes << "/" << mean_nbnodes << ")" << std::endl;
    
    
    // Populate some information about each color
    std::vector<Mesh::PartInfo> ColorInfo(nparts);
    for(Mesh::IndexType colind = 0;colind < Nparts;colind++){
      ColorInfo[colind].part    = colind+1;                      // partition id
      ColorInfo[colind].nelem   = dual_pcon[colind].size();      // number of elems
      ColorInfo[colind].nborder = colorneighbors[colind].size(); // number of borders
      ColorInfo[colind].nnodes  = colorsofnodes[colind].size();  // total num nodes
      ColorInfo[colind].nshared = bordercon[colind].size();      // num shared nodes
      ColorInfo[colind].nown    = dual_owner[colind].size();     // num shared nodes owned
    }

    global.FunctionEntry("Colorization");
    global.FunctionEntry("Colorization 1");
    // First, lets loop through the nodes and populate each color's
    // non-shared node list.
    Mesh::Connectivity PartitionedNodes;
    PartitionedNodes.resize(Nparts);
    Mesh::IndexVec psize;
    Primitive::MultiContainer<Mesh::IndexType,Mesh::IndexType>::VecMap pglob2loc;
    pglob2loc.resize(Nparts);
    psize.resize(Nparts,0);
    for(Mesh::IndexType jjj = 0;jjj < number_of_nodes;jjj++){
      Mesh::IndexType node_id = jjj+1;
      if(nodecolors.Esize(node_id) == 1){ // then it's not shared
	Mesh::IndexType thecolor = nodecolors[node_id-1][0];
	PartitionedNodes[thecolor-1].push_back(node_id);
	psize[thecolor-1]++;
	pglob2loc[thecolor-1].insert(std::make_pair(node_id,psize[thecolor-1]));
      }
    }
    // Now lets add each color's shared nodes that it actually owns
    // into the node list
    for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
      std::vector<Mesh::IndexType>::iterator doi = dual_owner[jjj].begin();
      while(doi != dual_owner[jjj].end()){
	Mesh::IndexType border_node_id = *doi++;
	Mesh::IndexType global_node_id = bn_index_map[border_node_id-1];
	PartitionedNodes[jjj].push_back(global_node_id);
	psize[jjj]++;
	pglob2loc[jjj].insert(std::make_pair(global_node_id,psize[jjj]));
      }
    } 
    //   
    // That's all we can conveniently do at this time.  We need
    // to do the following sections of code before we can populate
    // the rest of the nodes, which are nodes on each color that
    // are actually owned by another color
    //
    global.FunctionExit("Colorization 1");


    global.FunctionEntry("Colorization 2");
    // Initial prep of border description for each color
    Primitive::MultiContainer<Mesh::Border>::VecVec ColorBorders;
    Mesh::VecMap   colormap;
    ColorBorders.resize(nparts);
    colormap.resize(nparts);
    for(Mesh::IndexType iii = 0;iii<Nparts;iii++){
      ColorBorders[iii].resize(colorneighbors[iii].size());
      Mesh::IndexVec::iterator cni = colorneighbors[iii].begin();
      Mesh::IndexType colorind = 0;
      while(cni != colorneighbors[iii].end()){
	colormap[iii].insert(std::make_pair(*cni,colorind+1));
	ColorBorders[iii][colorind].rpart = *cni;
	colorind++;
	cni++;
      }
    }
    global.FunctionExit("Colorization 2");

    global.FunctionEntry("Colorization 3");
    // Populate each color's border structure
    // Loop through all colors
    for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
      Mesh::IndexType this_color = jjj+1;
      // Loop thru all border nodes this color owns
      Mesh::IndexVec::iterator nrcvIt = dual_owner[jjj].begin();
      while(nrcvIt != dual_owner[jjj].end()){
	Mesh::IndexType border_node_id = *nrcvIt++;
	Mesh::IndexType global_bn_id = bn_index_map[border_node_id-1];
	// Loop through the colors that touch this border node
	Mesh::IndexVec::iterator ncIt = nodecolors[global_bn_id-1].begin();
	while(ncIt != nodecolors[global_bn_id-1].end()){
	  Mesh::IndexType remote_color_id = *ncIt++;
	  if(remote_color_id != this_color){
	    // Add this node to the border data structure for the local color
	    Mesh::IndexType local_border_id = colormap[jjj][remote_color_id];
	    assert(ColorBorders[jjj][local_border_id-1].rpart == remote_color_id);
	    ColorBorders[jjj][local_border_id-1].nrecv.push_back(global_bn_id);
	    // Add this node to the border data structure for remote procs
	    Mesh::IndexType remote_border_id = colormap[remote_color_id-1][this_color];
	    assert(ColorBorders[remote_color_id-1][remote_border_id-1].rpart == this_color);
	    ColorBorders[remote_color_id-1][remote_border_id-1].nsend.push_back(global_bn_id);
	  }
	}
	
      }
    }
    global.FunctionExit("Colorization 3");

    // Now to continue populating the nodes for each color - we now have 
    // naturally processed the nodes on each color that are remotely 
    // owned, so we can conveniently add them to each color's nodelist
    // Loop over colors
    global.FunctionEntry("Colorization 4");
    for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
      std::list<Mesh::IndexType> sendnodes;
      // Loop over each color's border
      for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nborder;iii++){
	std::vector<Mesh::IndexType>::iterator bni = ColorBorders[jjj][iii].nsend.begin();
	while(bni != ColorBorders[jjj][iii].nsend.end())
	  sendnodes.push_back(*bni++);
      }    
      sendnodes.sort();
      sendnodes.unique();
      std::list<Mesh::IndexType>::iterator sni = sendnodes.begin();
      while(sni != sendnodes.end()){
	Mesh::IndexType send_node_id = *sni++;
	PartitionedNodes[jjj].push_back(send_node_id);
	psize[jjj]++;
	pglob2loc[jjj].insert(std::make_pair(send_node_id,psize[jjj]));
      }
    }
    global.FunctionExit("Colorization 4");
    // 
    // Now do some color validation
    global.FunctionEntry("Color validation");
    for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
      assert(psize[jjj] == ColorInfo[jjj].nnodes);
    }   
    global.FunctionExit("Color validation");
    global.FunctionExit("Colorization");
    global.FunctionEntry("Writing Partitions");
    // Loop through all colors
    global.FunctionEntry("Read Coords");
    Mesh::NodalCoordinates nc;
    Inf.open(MeshName.c_str());
    if(!Inf){
      if(ErrOut) 
	*ErrOut << "Failed to open " << MeshName << std::endl;
      return(1);
    }
    Inf >> nc;
    Inf.close();
    global.FunctionExit("Read Coords");
    std::vector<std::string> nodal_props(number_of_nodes);
    std::vector<std::string> elemental_props(number_of_elements);
    bool mesh_has_properties = false;
    std::string propname = MeshBaseName + ".prop";
    Inf.open(propname.c_str());
    if(Inf){
      mesh_has_properties = true;
      global.FunctionEntry("Read Properties");
      for(unsigned int propcount = 0;propcount < number_of_nodes;propcount++)
	std::getline(Inf,nodal_props[propcount]);
      std::getline(Inf,elemental_props[0]);
      for(unsigned int propcount = 0;propcount < number_of_elements;propcount++)
	std::getline(Inf,elemental_props[propcount]);
      Inf.close();
      global.FunctionExit("Read Properties");
    }
    for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
      std::ostringstream Ostr;
      Ostr << MeshBaseName << "." << jjj+1 << ".info";
      std::ofstream OutFile;
      std::ofstream OutProp;
      OutFile.open(Ostr.str().c_str());
      if(mesh_has_properties){
	Ostr.str("");
	Ostr << MeshBaseName << "." << jjj+1 << ".prop";
	OutProp.open(Ostr.str().c_str());
      }
      OutFile << Nparts << std::endl << jjj+1 << std::endl 
	      << ColorInfo[jjj].nelem << std::endl
	      << ColorInfo[jjj].nnodes << std::endl
	      << ColorInfo[jjj].nborder << std::endl
	      << ColorInfo[jjj].nshared << std::endl
	      << ColorInfo[jjj].nown << std::endl;
      OutFile.close();
      Ostr.str("");
      Ostr << MeshBaseName << "." << jjj+1 << ".pmesh";
      OutFile.open(Ostr.str().c_str());
      // Write the nodal coordinates
      OutFile << ColorInfo[jjj].nnodes << std::endl;
      for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nnodes;iii++){
	Mesh::IndexType node_id = PartitionedNodes[jjj][iii];
	OutFile << nc.x(node_id) << " " << nc.y(node_id) << " "
		<< nc.z(node_id) << std::endl;
	if(OutProp)
	  OutProp << nodal_props[node_id-1] << std::endl;
      }
      if(OutProp)
	OutProp << std::endl;
      // Write the element connectivity
      OutFile << ColorInfo[jjj].nelem << std::endl;
      for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nelem;iii++){
	Mesh::IndexType eid = dual_pcon[jjj][iii];
	std::vector<Mesh::IndexType>::iterator eni = econ[eid-1].begin();
	while(eni != econ[eid-1].end()){
	  OutFile << pglob2loc[jjj][*eni];
	  if(eni != econ[eid-1].end())
	    OutFile << " ";
	  eni++;
	}
	OutFile << std::endl;
	if(OutProp)
	  OutProp << elemental_props[eid-1] << std::endl;
      }
      if(OutProp)
	OutProp.close();
      // Write the borders
      OutFile << ColorInfo[jjj].nborder << std::endl;
      for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nborder;iii++){
	OutFile << ColorBorders[jjj][iii].rpart << " "
		<< ColorBorders[jjj][iii].nrecv.size() << " "
		<< ColorBorders[jjj][iii].nsend.size() << std::endl;
	std::vector<Mesh::IndexType>::iterator rni = ColorBorders[jjj][iii].nrecv.begin();
	while(rni != ColorBorders[jjj][iii].nrecv.end()){
	  OutFile << pglob2loc[jjj][*rni] << " " << *rni << std::endl;
	  rni++;
	}
	rni = ColorBorders[jjj][iii].nsend.begin();
	while(rni != ColorBorders[jjj][iii].nsend.end()){
	  OutFile << pglob2loc[jjj][*rni] << " " << *rni << std::endl;
	  rni++;
	}
      }
      OutFile.close();
    }
    global.FunctionExit("Writing Partitions");
    global.FunctionExit("Partitioning");
  }






  //
  // Renumbering test for elements 
  //
//   bool do_partition_elements = false;
//   if(do_partition_elements){
//     if(nproc > 1){
//       if(ErrOut)
// 	*ErrOut << "Error: partitioning is for serial runs." << std::endl;
//       comm.SetExit(1);
//     }
//     if(comm.Check())
//       return(1);
//     global.FunctionEntry("Partitioning");
//     std::ifstream Inf;
//     std::string MeshName(infiles[0]);
//     Inf.open(MeshName.c_str());
//     if(!Inf){
//       if(ErrOut) 
// 	*ErrOut << "Failed to open " << MeshName << std::endl;
//       return(1);
//     }
//     global.FunctionEntry("Read Mesh");
//     Mesh::IndexType number_of_nodes = 0;
//     if(true){
//       Mesh::NodalCoordinates nc;
//       Inf >> nc;
//       number_of_nodes = nc.Size();
//       nc.destroy();
//     }
//     if(StdOut && verblevel) 
//       *StdOut << "Number of nodes: " << number_of_nodes << std::endl;
//     Mesh::Connectivity econ;
//     Inf >> econ;
//     global.FunctionExit("Read Mesh");
//     Inf.close();
//     econ.ShrinkWrap();
//     Mesh::IndexType number_of_elements = econ.Nelem();
//     if(StdOut && verblevel) 
//       *StdOut << "Number of elements: " << number_of_elements
// 	      << std::endl;
//     Mesh::IndexType N = npartitions;
//     econ.SyncSizes();
//     // re-orient the elements if needed
//     if(do_reorient){ // don't do this for 2d geometry (yet)
//       Mesh::NodalCoordinates nc;
//       Inf.open(MeshName.c_str());
//       if(!Inf){
// 	if(ErrOut) 
// 	  *ErrOut << "Failed to open " << MeshName << std::endl;
// 	return(1);
//       }
//       global.FunctionEntry("ReadCoords");
//       Inf >> nc;
//       Inf.close();
//       global.FunctionExit("ReadCoords");
//       global.FunctionEntry("Element Reorientation");
//       for(Mesh::IndexType element_being_processed = 1;
// 	  element_being_processed <= number_of_elements;
// 	  element_being_processed++)
// 	{
// 	  Mesh::IndexType index = element_being_processed - 1;
// 	  Mesh::IndexType size_of_element = econ.Esize(element_being_processed);
// 	  Mesh::GenericElement ge(size_of_element);
// 	  if(ge.Inverted(econ[index],nc))
// 	    ge.ReOrient(econ[index]);
// 	}
//       global.FunctionExit("Element Reorientation");
//     }
//     //    Mesh::IndexType N2 = 0;
//     //    if(infiles.size() > 1){
//     //      std::istringstream Istr(infiles[1]);
//     //      Istr >> N;
//     //    }
//     //    if(infiles.size() > 2){
//     //      std::istringstream Istr(infiles[2]);
//     //      Istr >> N2;
//     //    }
//     Mesh::Connectivity dc;
//     global.FunctionEntry("DualCon");
//     econ.Inverse(dc,number_of_nodes);
//     global.FunctionExit("DualCon");
//     dc.ShrinkWrap();
//     Mesh::Connectivity nl2;
//     // Neighborhood populates an array of neighbors
//     // for each element.  A neighbor is any element
//     // that shares a node with a given element.
//     global.FunctionEntry("Neighborhood2");
//     econ.GetNeighborhood(nl2,dc);
//     global.FunctionExit("Neighborhood2"); 
//     nl2.ShrinkWrap();
//     // -------------- Partitioning --------------------
//     std::vector<idxtype> adj;
//     std::vector<idxtype> xadj;
//     global.FunctionEntry("Container2CSR");
//     Mesh::Container2CSR< std::vector< std::vector<Mesh::IndexType> >, 
//       std::vector<Mesh::IndexType>,
//       std::vector<idxtype>,idxtype>
//       (xadj,adj,nl2);
//     global.FunctionExit("Container2CSR");
//     //    nl2.destroy();
//     //    dc.destroy();
//     //    econ.destroy();
//     std::vector<idxtype> order(number_of_elements,0);
//     std::vector<idxtype> sizes(2,0);
//     idxtype nelem = static_cast<idxtype>(number_of_elements);
//     idxtype wgtflag = 0;
//     idxtype numflag = 0;
//     idxtype nparts = N;
//     idxtype options[5];
//     idxtype edgecuts = 0;
//     options[0] = 1;
//     options[1] = 3;
//     options[2] = 15;
//     std::vector<idxtype> vtxdist(2*Npartitions,0);
//     vtxdist[1] = number_of_elements;
//     std::vector<idxtype> part(number_of_elements,0);
//     global.FunctionEntry("Metis");
//     int comm = MPI_COMM_WORLD;
//     int ncon = 0;
//     float stupid = 0;
//     float ubvecstupid = 0;
//     if(array_checking){
//       *LogOut << "Vtxdist:" << std::endl;
//       DumpContents(*LogOut,vtxdist," ");
//       *LogOut << std::endl << "Xadj:" << std::endl;
//       DumpContents(*LogOut,xadj," ");
//       *LogOut << std::endl << "Adj:" << std::endl;
//       DumpContents(*LogOut,adj," ");
//       *LogOut << std::endl;
//     }
//     global.FunctionEntry("Renumbering");      
//     ParMETIS_V3_NodeND(&vtxdist[0],&xadj[0],&adj[0],&numflag,
// 		       options,&order[0],&sizes[0],&comm);
//     global.FunctionExit("Renumbering");      
//     if(array_checking){
//       *LogOut << "Order:" << std::endl;
//       DumpContents(*LogOut,order," ");
//       *LogOut << std::endl << "Sizes:" << std::endl;
//       DumpContents(*LogOut,sizes," ");
//       *LogOut << std::endl;
//     }
//     ParMETIS_V3_PartKway(&vtxdist[0],&xadj[0],&adj[0],NULL,NULL,&wgtflag,
// 			 &numflag,NULL,&nparts,NULL,NULL,
// 			 options,&edgecuts,&part[0],&comm);
//     //    METIS_PartGraphKway(&vtxdist[0],&xadj[0],&adj[0],NULL,NULL,&wgtflag,
//     //			&numflag,&nparts,options,&edgecuts,&part[0],&comm);
//     global.FunctionExit("Metis");
//     // ------------------------------------------------
//     // Take part and create the 1-shifted version into 
//     // a multicontainer
//     Mesh::Connectivity pcon;
//     pcon.resize(number_of_elements);
//     for(Mesh::IndexType i = 0;i < number_of_elements;i++)
//       pcon[i].push_back(part[i]+1);
//     pcon.Sync();
//     pcon.SyncSizes();
//     //    std::ofstream Ouf;
//     //    Ouf.open("part");
//     //    Ouf << pcon << std::endl;
//     //    Ouf.close();
//     if(StdOut && verblevel) 
//       *StdOut << "Number of edgecuts: " << edgecuts << std::endl;
//     if(false){
//       if(StdOut && verblevel) 
// 	*StdOut << "Partition Table: ";
//       DumpContents(*StdOut,part," ");
//       if(StdOut && verblevel) 
// 	*StdOut << std::endl;
//     }
//     part.resize(0);
//     part.swap(part);
//     Mesh::IndexType Nparts = nparts;

//    // Invert the partition table
//     global.FunctionEntry("Dual Pcon");
//     Mesh::Connectivity dual_pcon;  // E[C] -> C[E]
//     pcon.Inverse(dual_pcon,nparts);
//     global.FunctionExit("Dual Pcon");
//     dual_pcon.SyncSizes();

//     global.FunctionEntry("GetAdj nodecolors");
//     Mesh::Connectivity nodecolors; // N[E]%E[C] = N[C]
//     dc.GetAdjacent(nodecolors,pcon,nparts,true);
//     nodecolors.SyncSizes();
//     global.FunctionExit("GetAdj nodecolors");
    
//     global.FunctionEntry("Inverse colornodes");
//     Mesh::Connectivity colorsofnodes; // N[C] -> C[N]
//     nodecolors.Inverse(colorsofnodes,nparts);
//     global.FunctionExit("Inverse colornodes");
//     colorsofnodes.SyncSizes();

//     global.FunctionEntry("ColorNeighbors");
//     Mesh::Connectivity colorneighbors; // C[N]%N[C] = C[C]
//     colorsofnodes.GetNeighborhood(colorneighbors,nodecolors,true,true);
//     global.FunctionExit("ColorNeighbors");
//     colorneighbors.SyncSizes();

//     // This is enough information to fully describe the partitioning
//     // to each partition and form the communication lists.  We can
//     // write this to file, or just pass it along in messages to 
//     // each partition.

//     // For each color, show me the taco
//     // produces C[BN] (border nodes for each color)
//     Mesh::Connectivity bordercon;
//     nodecolors.InverseDegenerate(bordercon,nparts);
//     // Lets determine and report the communication imbalance
//     bool report_each = false;
//     if(report_each)
//       if(StdOut && verblevel) 
// 	*StdOut << "Nshared nodes/partition:" << std::endl;
//     double mean_nbnodes = 0;
//     Mesh::IndexType max_nbnodes = 0;
//     Mesh::IndexType min_nbnodes = number_of_nodes;
//     for(Mesh::IndexType iii = 0;iii < Nparts;iii++){
//       Mesh::IndexType nbsize = bordercon[iii].size();
//       if(report_each && StdOut)
// 	*StdOut << "Partition " << iii+1 << ": "
// 		<< nbsize << std::endl;
//       mean_nbnodes += static_cast<double>(nbsize);
//       if(nbsize > max_nbnodes) max_nbnodes = nbsize;
//       if(nbsize < min_nbnodes) min_nbnodes = nbsize;
//     }
//     mean_nbnodes /= static_cast<double>(nparts);
//     if(StdOut && verblevel) 
//       *StdOut << "(Max/Min/Mean) shared nodes/partition: (" << max_nbnodes << "/"
// 	      << min_nbnodes << "/" << mean_nbnodes << ")" << std::endl;

//     Mesh::IndexVec bn_index_map;
//     std::map<Mesh::IndexType,Mesh::IndexType> bnodemap;
//     std::list<Mesh::IndexType> bnodes;
//     global.FunctionEntry("Flatten");
//     Flatten<Mesh::Connectivity,std::vector<Mesh::IndexType>,std::list<Mesh::IndexType> >
//       (bordercon,bnodes);
//     global.FunctionExit("Flatten");
//     global.FunctionEntry("SortUnique");
//     bnodes.sort();
//     bnodes.unique();
//     global.FunctionExit("SortUnique");
//     Mesh::IndexType number_of_border_nodes = bnodes.size();
//     bn_index_map.resize(number_of_border_nodes);
//     Mesh::IndexType i = 0;
//     std::list<Mesh::IndexType>::iterator eli = bnodes.begin();
//     while(eli != bnodes.end()){
//       bn_index_map[i] = *eli;
//       bnodemap.insert(std::make_pair(*eli,i+1));
//       eli++;
//       i++;
//     }
    
//     global.FunctionEntry("MapBorderNodes");
//     // C[BN] -> BN[C], but need BN to be mapped in 1-NBN range
//     Mesh::Connectivity bordercon_mapped;    
//     MapElements<Mesh::Connectivity,Mesh::Connectivity,
//       std::vector<Mesh::IndexType>,std::map<Mesh::IndexType,Mesh::IndexType> >
//       (bordercon,bordercon_mapped,bnodemap);
//     global.FunctionExit("MapBorderNodes");
//     Mesh::Connectivity bcon_dual;
//     bordercon_mapped.Sync();
//     bordercon_mapped.Inverse(bcon_dual,number_of_border_nodes);
//     // Now we have BN[C], which is a list of colors for each border
//     // node.  We need to step through and assign each border node 
//     // an actual owner.  This will produce a list
//     // of which color (only 1 of them) is the owner of a particular BN.
//     // Or, bn[C], where bn are the subset of BN that C owns.
//     // We do this in a way so as to equalize the number of communicated
//     // nodes on each processors.
//     std::srand((unsigned)std::time(0));
//     Mesh::Connectivity owner;
//     std::vector<Mesh::IndexType> ncomnodes(nparts,0);
//     owner.resize(number_of_border_nodes);
//     Mesh::Connectivity::iterator bni = bcon_dual.begin();
//     Mesh::IndexType ownind = 0;
//     while(bni != bcon_dual.end()){
//       Mesh::IndexType ncolors = bni->size();
//       Mesh::IndexType color = (std::rand()%ncolors)+1;
//       for(Mesh::IndexType ci = 0;ci < ncolors;ci++)
// 	if(ncomnodes[(*bni)[color-1] - 1] > ncomnodes[(*bni)[ci]-1])
// 	  color = ci+1;
//       std::vector<Mesh::IndexType>::iterator bnici = bni->begin();
//       for(Mesh::IndexType ci = 1;ci < color;ci++)
// 	bnici++;
//       //      if(bordercon[*bnici-1].size() > static_cast<Mesh::IndexType>(mean_nbnodes))
//       //	color = (std::rand()%ncolors)+1;
//       owner[ownind++].push_back(*bnici);
//       ncomnodes[*bnici-1]++;
//       bni++;
//     }
//     // bn[C] -> C[bn]
//     Mesh::Connectivity dual_owner;
//     owner.Sync();
//     owner.Inverse(dual_owner,nparts);
//     mean_nbnodes = 0;
//     min_nbnodes = number_of_nodes;
//     max_nbnodes = 0;
//     for(Mesh::IndexType iii = 0;iii < Nparts;iii++){
//       Mesh::IndexType nbsize = dual_owner[iii].size();
//       if(report_each)
// 	if(StdOut && verblevel) 
// 	  *StdOut << "Partition " << iii+1 << ": "
// 		  << nbsize << std::endl;
//       mean_nbnodes += static_cast<double>(nbsize);
//       if(nbsize > max_nbnodes) max_nbnodes = nbsize;
//       if(nbsize < min_nbnodes) min_nbnodes = nbsize;
//     }
//     mean_nbnodes /= static_cast<double>(nparts);
//     if(StdOut && verblevel) 
//       *StdOut << "(Max/Min/Mean) recvd nodes/partition: (" << max_nbnodes << "/"
// 	      << min_nbnodes << "/" << mean_nbnodes << ")" << std::endl;
    
    
//     // Populate some information about each color
//     std::vector<Mesh::PartInfo> ColorInfo(nparts);
//     for(Mesh::IndexType colind = 0;colind < Nparts;colind++){
//       ColorInfo[colind].part    = colind+1;                      // partition id
//       ColorInfo[colind].nelem   = dual_pcon[colind].size();      // number of elems
//       ColorInfo[colind].nborder = colorneighbors[colind].size(); // number of borders
//       ColorInfo[colind].nnodes  = colorsofnodes[colind].size();  // total num nodes
//       ColorInfo[colind].nshared = bordercon[colind].size();      // num shared nodes
//       ColorInfo[colind].nown    = dual_owner[colind].size();     // num shared nodes owned
//     }

//     global.FunctionEntry("Colorization");
//     global.FunctionEntry("Colorization 1");
//     // First, lets loop through the nodes and populate each color's
//     // non-shared node list.
//     Mesh::Connectivity PartitionedNodes;
//     PartitionedNodes.resize(Nparts);
//     Mesh::IndexVec psize;
//     Primitive::MultiContainer<Mesh::IndexType,Mesh::IndexType>::VecMap pglob2loc;
//     pglob2loc.resize(Nparts);
//     psize.resize(Nparts,0);
//     for(Mesh::IndexType jjj = 0;jjj < number_of_nodes;jjj++){
//       Mesh::IndexType node_id = jjj+1;
//       if(nodecolors.Esize(node_id) == 1){ // then it's not shared
// 	Mesh::IndexType thecolor = nodecolors[node_id-1][0];
// 	PartitionedNodes[thecolor-1].push_back(node_id);
// 	psize[thecolor-1]++;
// 	pglob2loc[thecolor-1].insert(std::make_pair(node_id,psize[thecolor-1]));
//       }
//     }
//     // Now lets add each color's shared nodes that it actually owns
//     // into the node list
//     for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
//       std::vector<Mesh::IndexType>::iterator doi = dual_owner[jjj].begin();
//       while(doi != dual_owner[jjj].end()){
// 	Mesh::IndexType border_node_id = *doi++;
// 	Mesh::IndexType global_node_id = bn_index_map[border_node_id-1];
// 	PartitionedNodes[jjj].push_back(global_node_id);
// 	psize[jjj]++;
// 	pglob2loc[jjj].insert(std::make_pair(global_node_id,psize[jjj]));
//       }
//     } 
//     //   
//     // That's all we can conveniently do at this time.  We need
//     // to do the following sections of code before we can populate
//     // the rest of the nodes, which are nodes on each color that
//     // are actually owned by another color
//     //
//     global.FunctionExit("Colorization 1");


//     global.FunctionEntry("Colorization 2");
//     // Initial prep of border description for each color
//     Primitive::MultiContainer<Mesh::Border>::VecVec ColorBorders;
//     Mesh::VecMap   colormap;
//     ColorBorders.resize(nparts);
//     colormap.resize(nparts);
//     for(Mesh::IndexType iii = 0;iii<Nparts;iii++){
//       ColorBorders[iii].resize(colorneighbors[iii].size());
//       Mesh::IndexVec::iterator cni = colorneighbors[iii].begin();
//       Mesh::IndexType colorind = 0;
//       while(cni != colorneighbors[iii].end()){
// 	colormap[iii].insert(std::make_pair(*cni,colorind+1));
// 	ColorBorders[iii][colorind].rpart = *cni;
// 	colorind++;
// 	cni++;
//       }
//     }
//     global.FunctionExit("Colorization 2");

//     global.FunctionEntry("Colorization 3");
//     // Populate each color's border structure
//     // Loop through all colors
//     for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
//       Mesh::IndexType this_color = jjj+1;
//       // Loop thru all border nodes this color owns
//       Mesh::IndexVec::iterator nrcvIt = dual_owner[jjj].begin();
//       while(nrcvIt != dual_owner[jjj].end()){
// 	Mesh::IndexType border_node_id = *nrcvIt++;
// 	Mesh::IndexType global_bn_id = bn_index_map[border_node_id-1];
// 	// Loop through the colors that touch this border node
// 	Mesh::IndexVec::iterator ncIt = nodecolors[global_bn_id-1].begin();
// 	while(ncIt != nodecolors[global_bn_id-1].end()){
// 	  Mesh::IndexType remote_color_id = *ncIt++;
// 	  if(remote_color_id != this_color){
// 	    // Add this node to the border data structure for the local color
// 	    Mesh::IndexType local_border_id = colormap[jjj][remote_color_id];
// 	    assert(ColorBorders[jjj][local_border_id-1].rpart == remote_color_id);
// 	    ColorBorders[jjj][local_border_id-1].nrecv.push_back(global_bn_id);
// 	    // Add this node to the border data structure for remote procs
// 	    Mesh::IndexType remote_border_id = colormap[remote_color_id-1][this_color];
// 	    assert(ColorBorders[remote_color_id-1][remote_border_id-1].rpart == this_color);
// 	    ColorBorders[remote_color_id-1][remote_border_id-1].nsend.push_back(global_bn_id);
// 	  }
// 	}
	
//       }
//     }
//     global.FunctionExit("Colorization 3");

//     // Now to continue populating the nodes for each color - we now have 
//     // naturally processed the nodes on each color that are remotely 
//     // owned, so we can conveniently add them to each color's nodelist
//     // Loop over colors
//     global.FunctionEntry("Colorization 4");
//     for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
//       std::list<Mesh::IndexType> sendnodes;
//       // Loop over each color's border
//       for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nborder;iii++){
// 	std::vector<Mesh::IndexType>::iterator bni = ColorBorders[jjj][iii].nsend.begin();
// 	while(bni != ColorBorders[jjj][iii].nsend.end())
// 	  sendnodes.push_back(*bni++);
//       }    
//       sendnodes.sort();
//       sendnodes.unique();
//       std::list<Mesh::IndexType>::iterator sni = sendnodes.begin();
//       while(sni != sendnodes.end()){
// 	Mesh::IndexType send_node_id = *sni++;
// 	PartitionedNodes[jjj].push_back(send_node_id);
// 	psize[jjj]++;
// 	pglob2loc[jjj].insert(std::make_pair(send_node_id,psize[jjj]));
//       }
//     }
//     global.FunctionExit("Colorization 4");
//     // 
//     // Now do some color validation
//     global.FunctionEntry("Color validation");
//     for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
//       assert(psize[jjj] == ColorInfo[jjj].nnodes);
//     }   
//     global.FunctionExit("Color validation");
//     global.FunctionExit("Colorization");
//     global.FunctionEntry("Writing Partitions");
//     // Loop through all colors
//     global.FunctionEntry("Read Coords");
//     Mesh::NodalCoordinates nc;
//     Inf.open(MeshName.c_str());
//     if(!Inf){
//       if(ErrOut) 
// 	*ErrOut << "Failed to open " << MeshName << std::endl;
//       return(1);
//     }
//     Inf >> nc;
//     Inf.close();
//     global.FunctionExit("Read Coords");
//     for(Mesh::IndexType jjj = 0;jjj < Nparts;jjj++){
//       std::ostringstream Ostr;
//       Ostr << MeshName << "." << jjj+1 << ".info";
//       std::ofstream OutFile;
//       OutFile.open(Ostr.str().c_str());
//       OutFile << ColorInfo[jjj].nelem << std::endl
// 	      << ColorInfo[jjj].nnodes << std::endl
// 	      << ColorInfo[jjj].nborder << std::endl
// 	      << ColorInfo[jjj].nshared << std::endl
// 	      << ColorInfo[jjj].nown << std::endl;
//       OutFile.close();
//       Ostr.str("");
//       Ostr << MeshName << "." << jjj+1 << ".pmesh";
//       OutFile.open(Ostr.str().c_str());
//       // Write the nodal coordinates
//       OutFile << ColorInfo[jjj].nnodes << std::endl;
//       for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nnodes;iii++){
// 	Mesh::IndexType node_id = PartitionedNodes[jjj][iii];
// 	OutFile << nc.x(node_id) << " " << nc.y(node_id) << " "
// 		<< nc.z(node_id) << std::endl;
//       }
//       // Write the element connectivity
//       OutFile << ColorInfo[jjj].nelem << std::endl;
//       for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nelem;iii++){
// 	Mesh::IndexType eid = dual_pcon[jjj][iii];
// 	std::vector<Mesh::IndexType>::iterator eni = econ[eid-1].begin();
// 	while(eni != econ[eid-1].end()){
// 	  OutFile << pglob2loc[jjj][*eni];
// 	  if(eni != econ[eid-1].end())
// 	    OutFile << " ";
// 	  eni++;
// 	}
// 	OutFile << std::endl;
//       }
//       // Write the borders
//       OutFile << ColorInfo[jjj].nborder << std::endl;
//       for(Mesh::IndexType iii = 0;iii < ColorInfo[jjj].nborder;iii++){
// 	OutFile << ColorBorders[jjj][iii].rpart << " "
// 		<< ColorBorders[jjj][iii].nrecv.size() << " "
// 		<< ColorBorders[jjj][iii].nsend.size() << std::endl;
// 	std::vector<Mesh::IndexType>::iterator rni = ColorBorders[jjj][iii].nrecv.begin();
// 	while(rni != ColorBorders[jjj][iii].nrecv.end()){
// 	  OutFile << pglob2loc[jjj][*rni] << " " << *rni << std::endl;
// 	  rni++;
// 	}
// 	rni = ColorBorders[jjj][iii].nsend.begin();
// 	while(rni != ColorBorders[jjj][iii].nsend.end()){
// 	  OutFile << pglob2loc[jjj][*rni] << " " << *rni << std::endl;
// 	  rni++;
// 	}
//       }
//       OutFile.close();
//     }
//     global.FunctionExit("Writing Partitions");
//     global.FunctionExit("Partitioning");
//  }
  //  bool do_partitioning = true;
  if(do_clone){
    if(nproc > 1){
      if(ErrOut)
	*ErrOut << "Error: Cloning is a serial operation." << std::endl;
      comm.SetExit(1);
    }
    if(infiles.empty()){
      if(ErrOut)
	*ErrOut << "Error: Nothing to clone." << std::endl;
      comm.SetExit(1);
    }
    if(comm.Check())
      return(1);
    global.FunctionEntry("Cloning");
    std::ifstream Inf;
    std::string MeshName(infiles[0]);
    unsigned int N = 0;
    if(infiles.size() > 1){
      std::istringstream Instr(infiles[1]);
      Instr >> N;
    }
    if(N > 0)
      npartitions = N;
    if(StdOut && verblevel)
      *StdOut << "Cloning " << MeshName << " into " << npartitions << " partitions."
	      << std::endl;
    Inf.open(MeshName.c_str());
    if(!Inf){
      if(ErrOut) 
	*ErrOut << "Failed to open " << MeshName << std::endl;
      return(1);
    }
    global.FunctionEntry("Read Mesh");
    Mesh::IndexType number_of_nodes = 0;
    Mesh::NodalCoordinates nc;
    Inf >> nc;
    number_of_nodes = nc.Size();
    std::vector<double> bounds(6,0.0);
    Mesh::GetCoordinateBounds(nc,bounds);
    if(StdOut && verblevel)
      *StdOut << "Number of nodes: " << number_of_nodes << std::endl
	      << "Mesh Bounds: (" << bounds[0] << "," << bounds[1] 
	      << ") (" << bounds[2] << "," << bounds[3] 
	      << ") (" << bounds[4] << "," << bounds[5] << ")" << std::endl;
    double Xmin = bounds[0];
    double Xmax = bounds[1];
    Mesh::Connectivity econ;
    Inf >> econ;
    Inf.close();
    econ.ShrinkWrap();
    econ.SyncSizes();
    global.FunctionExit("Read Mesh");
    Mesh::IndexType number_of_elements = econ.Nelem();
    if(StdOut && verblevel) 
      *StdOut << "Number of elements: " << number_of_elements
	      << std::endl;
    global.FunctionEntry("Forming Clones");
    // Identify the boundary nodes
    std::vector<Mesh::IndexType> left_owned;
    std::vector<Mesh::IndexType> left_remote;
    std::vector<Mesh::IndexType> right_owned;
    std::vector<Mesh::IndexType> right_remote;
    std::vector<Mesh::IndexType> non_boundary;
    bool left_polarity = true;
    bool right_polarity = false;
    for(Mesh::IndexType nn = 1;nn <= number_of_nodes;nn++){
      if(nc.x(nn) == Xmin){
	if(left_polarity){
	  left_owned.push_back(nn);
	  left_polarity = false;
	}
	else{
	  left_remote.push_back(nn);
	  left_polarity = true;
	}
      }
      else if(nc.x(nn) == Xmax){
	if(right_polarity){
	  right_owned.push_back(nn);
	  right_polarity = false;
	}
	else{
	  right_remote.push_back(nn);
	  right_polarity = true;
	}
      }
      else
	non_boundary.push_back(nn);
    }
    Mesh::IndexType nnbn   = non_boundary.size();
    Mesh::IndexType nleft  = left_owned.size() + left_remote.size();
    Mesh::IndexType nright = right_owned.size() + right_remote.size();
    Mesh::IndexType nshared_nodes = nleft + nright;
    N = npartitions;
    Mesh::IndexType nown_left  = left_owned.size();
    Mesh::IndexType nown_right = right_owned.size();
    if(StdOut && verblevel)
      *StdOut << "Number of internal nodes/partition: " << nnbn << std::endl
	      << "Number of shared nodes/partition: " << nshared_nodes << std::endl;
    // Map the nodes to new id's based on the partitioning
    global.FunctionEntry("Mapping Clone Genes");
    std::vector<Mesh::IndexType> NodeMap(number_of_nodes,0);
    std::vector<Mesh::IndexType>::iterator ni = non_boundary.begin();
    Mesh::IndexType node_id = 1;
    while(ni != non_boundary.end())
      NodeMap[*ni++ - 1] = node_id++;
    ni = left_owned.begin();
    while(ni != left_owned.end())
      NodeMap[*ni++ - 1] = node_id++;
    ni = right_owned.begin();
    while(ni != right_owned.end())
      NodeMap[*ni++ - 1] = node_id++;
    ni = left_remote.begin();
    while(ni != left_remote.end())
      NodeMap[*ni++ - 1] = node_id++;
    ni = right_remote.begin();
    while(ni != right_remote.end())
      NodeMap[*ni++ - 1] = node_id++;
    node_id--;
    assert(node_id == number_of_nodes);
    global.FunctionExit("Mapping Clone Genes");
    global.FunctionExit("Forming Clones");
    global.FunctionEntry("Writing Clones");
    for(Mesh::IndexType jjj = 0;jjj < N;jjj++){
      std::ostringstream Ostr;
      Ostr << MeshName << "." << jjj+1 << ".info";
      std::ofstream OutFile;
      OutFile.open(Ostr.str().c_str());
      OutFile << N << std::endl << jjj+1 << std::endl
	      << number_of_elements << std::endl
	      << number_of_nodes << std::endl
	      << 2 << std::endl
	      << nshared_nodes << std::endl
	      << nown_left+nown_right << std::endl;
      OutFile.close();
      Ostr.str("");
      Ostr << MeshName << "." << jjj+1 << ".pmesh";
      OutFile.open(Ostr.str().c_str());
      // Write the nodal coordinates
      OutFile << number_of_nodes << std::endl;
      std::vector<Mesh::IndexType>::iterator nbni = non_boundary.begin();
      while(nbni != non_boundary.end()){
	OutFile << nc.x(*nbni) << " " << nc.y(*nbni) << " "
		<< nc.z(*nbni) << std::endl;
	nbni++;
      }
      nbni = left_owned.begin();
      while(nbni != left_owned.end()){
	OutFile << nc.x(*nbni) << " " << nc.y(*nbni) << " "
		<< nc.z(*nbni) << std::endl;
	nbni++;
      }
      nbni = right_owned.begin();
      while(nbni != right_owned.end()){
	OutFile << nc.x(*nbni) << " " << nc.y(*nbni) << " "
		<< nc.z(*nbni) << std::endl;
	nbni++;
      }
      nbni = left_remote.begin();
      while(nbni != left_remote.end()){
	OutFile << nc.x(*nbni) << " " << nc.y(*nbni) << " "
		<< nc.z(*nbni) << std::endl;
	nbni++;
      }
      nbni = right_remote.begin();
      while(nbni != right_remote.end()){
	OutFile << nc.x(*nbni) << " " << nc.y(*nbni) << " "
		<< nc.z(*nbni) << std::endl;
	nbni++;
      }
      // Write the element connectivity
      OutFile << number_of_elements << std::endl;
      for(Mesh::IndexType iii = 0;iii < number_of_elements;iii++){
	std::vector<Mesh::IndexType>::iterator eni = econ[iii].begin();
	while(eni != econ[iii].end()){
	  OutFile << NodeMap[*eni-1];
	  if(eni != econ[iii].end())
	    OutFile << " ";
	  eni++;
	}
	OutFile << std::endl;
      }
      // Write the borders
      OutFile << 2 << std::endl;
      OutFile << (jjj == 0 ? N : jjj) << " "
	      << nown_left << " "
	      << nleft - nown_left << std::endl;
      std::vector<Mesh::IndexType>::iterator rni = left_owned.begin();
      while(rni != left_owned.end()){
	  OutFile << NodeMap[*rni-1] << " " << *rni << std::endl;
	  rni++;
      }
      rni = left_remote.begin();
      while(rni != left_remote.end()){
	OutFile << NodeMap[*rni-1] << " " << *rni << std::endl;
	rni++;
      }
      OutFile << (jjj == (N-1) ? 1 : jjj+2) << " "
	      << nown_right << " "
	      << nright - nown_right << std::endl;
      rni = right_owned.begin();
      while(rni != right_owned.end()){
	  OutFile << NodeMap[*rni-1] << " " << *rni << std::endl;
	  rni++;
      }
      rni = right_remote.begin();
      while(rni != right_remote.end()){
	OutFile << NodeMap[*rni-1] << " " << *rni << std::endl;
	rni++;
      }
      OutFile.close();
    }
    global.FunctionExit("Writing Clones");
    global.FunctionExit("Cloning");
  }
  comm.Barrier();
  if(LogOut && debug)
    *LogOut << "All processors made it to the end." << std::endl;
  global.FunctionExit("main");
  comm.Barrier();
  if(LogOut && debug)
    *LogOut << "All processors made it to the end." << std::endl;
  global.Finalize();
  //  global.Profiler.Finalize();
  if(StdOut && verblevel)
    //    global.Report(Profiler.SummarizeSerialExecution(*StdOut);
    global.Report(*StdOut);
  return(0);
}

// AssembleOrderedRow(DOFType &row,DOFContType &j,DOFDataContType &dofdat,KType &k) {
//     typename DOFContType::iterator jIt = j.begin();
//     typename DOFDataContType::iterator ddIt = dofdat.begin();
//     if(*jIt == 0)
//       jIt++;
//     Mesh::IndexType kindex = k.get_index(row,j);
//     Mesh::IndexType kend   = k.get_index(row+1);
//     while(kindex < kend){
//       Mesh::IndexType jind = *jIt;
//       if(k.get_dof(kindex) == jind){
//          k[kindex++] = *ddIt++;
//          jIt++;
//       }
//       kindex++;
//     }
// }

  class TestBase {
  private:
    int private_data;
  protected:
    int protected_data;
  public:
    virtual ~TestBase(){};
    virtual int f1(){ return 0; };
    virtual int f2(){ return 0; };
  };
  
  
  template <typename DataType>
    class TestDerived1 : public TestBase
    {
    protected:
      std::vector<DataType> DerivedNative;
    };
  template <typename DataType>
    class TestDerived2 : public TestBase
    {
    protected:
      std::vector<DataType> DerivedNative;
    public:
      int f3() { return 0; };
    };
  class TestNativeType1 {
  private:
  int private_data;
public:
  int public_data;
};
class TestNativeType2 {
private:
  int private_data;
public:
  int public_data;
};
int TestFunction()
{
  typedef TestBase* TestBasePtr;
  std::vector<TestBasePtr> TestVector;
  TestDerived1<TestNativeType1> t11;
  TestDerived1<TestNativeType2> t12;
  TestDerived2<TestNativeType1> t21;
  TestVector.push_back(&t11);
  TestVector.push_back(&t12);
  TestVector.push_back(&t21);
  return(0);
}
