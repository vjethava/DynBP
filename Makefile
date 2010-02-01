CXX=g++
# CXXFLAGS=-O2 -g -I./inc -DDEBUG -I/home/vjethava/share/include
# LDFLAGS=-L/home/vjethava/share/lib -lm -lgslcblas -lgsl -lcxcore -lcv -lhighgui -lcvaux -lml
CXXFLAGS=-O2 -g -I./inc -DDEBUG 
LDFLAGS=-L/usr/local/lib -lm -lgslcblas -lgsl -lcxcore -lcv -lhighgui -lcvaux -lml
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
TOBJDIR=./obj/test
BINDIR=./bin
DATADIR=./data
SOBJDIR=./obj/src
TSDIR=./test
SSDIR=./src
MAKE=make

vpath %.cpp $(TSDIR)
vpath %.cpp $(SSDIR)
vpath %.o $(TOBJDIR)
vpath %.o $(SOBJDIR)

all:
	cd $(SSDIR); $(MAKE) all; 
	cd $(TSDIR); $(MAKE) clean; make ptrue;
	@echo "Finished make all"

clean:
	rm -rf $(SOBJDIR)/*.o $(TOBJDIR)/*.o $(BINDIR)/video1
	@echo "Finished clean" 

potts:
	@rm -rf $(BINDIR)/potts
	cd $(SSDIR); $(MAKE) all; 
	cd $(TSDIR); $(MAKE) clean; make potts;
	@echo "Finished make all"

tracker: 
	@rm -rf $(BINDIR)/potts
	cd $(SSDIR); $(MAKE) all; 
	cd $(TSDIR); $(MAKE) clean; make tracker;
	@echo "Finished make all"

testAMPL:
	@rm -rf $(BINDIR)/testAMPL
	cd $(SSDIR); $(MAKE) all;
	cd $(TSDIR); $(MAKE) testAMPL;
	@echo "Finished make testAMPL"


ptrue:
	@rm -rf $(BINDIR)/ptrue
	cd $(SSDIR); $(MAKE) all;
	cd $(TSDIR); $(MAKE) ptrue;
	@echo "Finished make ptrue"

MVI:
	@rm -rf $(BINDIR)/MVI
	cd $(SSDIR); $(MAKE) all;
	cd $(TSDIR); $(MAKE) clean; $(MAKE) MVI;
	@echo "Finished make MVI"

testGraphBP:
	@rm -rf $(BINDIR)/testGraphBP
	cd $(SSDIR); $(MAKE) all;
	cd $(TSDIR); $(MAKE) clean; $(MAKE) testGraphBP;
	@echo "Finished make testGraphBP"
	
testGraphGBP: 
	@rm -rf $(BINDIR)/testGraphGBP
	cd $(SSDIR); $(MAKE) all;
	cd $(TSDIR); $(MAKE) clean; $(MAKE) testGraphGBP;
	@echo "Finished make testGraphGBP"	

GraphGBP: testGraphGBP 
	@echo "Running testGraphGBP"
	$(BINDIR)/testGraphGBP 