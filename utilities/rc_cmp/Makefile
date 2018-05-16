OBJS = \
          main.o file.o partition.o points.o fileselect.o util.o metric.o interpolation.o tecplotwriter.o

EXENAME = rc_cmp

CC = g++
#CC = /opt/mpich-mx-gcc-32/bin/mpicxx
CCOPTS = -g -Wno-deprecated -Wno-write-strings
LINK = g++
#LINK = /opt/mpich-mx-gcc-32/bin/mpicxx
LINKOPTS = -o $(EXENAME)

$(EXENAME):  $(OBJS)
	$(LINK) $(LINKOPTS) $(OBJS)

clean:
	-rm *.o $(EXENAME)

util.o : util.cpp util.h
	$(CC) -c $(CCOPTS) util.cpp

tecplotwriter.o : tecplotwriter.cpp tecplotwriter.h file.h points.h datatypedef.h
	$(CC) -c $(CCOPTS) tecplotwriter.cpp

metric.o : metric.cpp metric.h file.cpp file.h CImg.h colormaps.h datatypedef.h
	$(CC) -c $(CCOPTS) metric.cpp

fileselect.o : fileselect.cpp fileselect.h file.cpp file.h
	$(CC) -c $(CCOPTS) fileselect.cpp

points.o : points.cpp points.h datatypedef.h
	$(CC) -c $(CCOPTS) points.cpp

partition.o : partition.cpp partition.h points.cpp points.h util.cpp util.h datatypedef.h
	$(CC) -c $(CCOPTS) partition.cpp

interpolation.o : interpolation.cpp interpolation.h points.cpp points.h util.cpp util.h datatypedef.h
	$(CC) -c $(CCOPTS) interpolation.cpp

file.o : file.cpp file.h points.cpp points.h partition.cpp partition.h util.cpp util.h datatypedef.h
	$(CC) -c $(CCOPTS) file.cpp

main.o : main.cpp main.h file.cpp file.h fileselect.cpp fileselect.h metric.cpp metric.h interpolation.cpp interpolation.h
	$(CC) -c $(CCOPTS) main.cpp

