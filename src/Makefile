CXX=g++
# LDFLAGS=-L/usr/local/lib -lm -lgslcblas -lgsl -lcxcore -lcv -lhighgui -lcvaux -lml
# CXXFLAGS=-O2 -g -I../inc -DDEBUG


CXXFLAGS=-O2 -g -I../inc -I/home/vjethava/share/include -DDEBUG 
LDFLAGS=-L/home/vjethava/share/lib -lm -lgslcblas -lgsl -lcxcore -lcv -lhighgui -lcvaux -lml

SRCS=state.cpp belief.cpp node.cpp bpgraph.cpp myedge.cpp mygraph.cpp main.cpp ising.cpp VideoReader.cpp VideoGraph.cpp PottsGraph.cpp PottsModel.cpp PottsCB.cpp VR.cpp MVI.cpp GraphBP.cpp 
OBJS=$(SRCS:.cpp=.o) 
OBJDIR=../obj/src
BINDIR=../bin
DATADIR=../data
INCDIR=../inc
vpath %.h $(INCDIR)
vpath %.o $(OBJDIR)

all: $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	@echo "finished compilation"
	
%.o: %.cpp
	@echo "compiling files $< to $(OBJDIR)/$@"
	@$(CXX) -c $(CXXFLAGS) $< -o $(OBJDIR)/$@
	
%.cpp: %.h

clean: 
	@echo "cleaning the $(OBJDIR) directory"
	rm -rf $(OBJDIR)/*.o

