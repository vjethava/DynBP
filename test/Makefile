CXX=g++
CXXFLAGS=-O2 -g -I../inc -I/home/vjethava/share/include -DDEBUG 
LDFLAGS=-L/home/vjethava/share/lib -lm -lgslcblas -lgsl -lcxcore -lcv -lhighgui -lcvaux -lml
#CXXFLAGS=-O2 -g -I../inc -DDEBUG 
#LDFLAGS=-L/usr/local/lib -lm -lgslcblas -lgsl -lcxcore -lcv -lhighgui -lcvaux -lml
SRCS=$(wildcard *.cpp)
#boost_regex_example.cpp family_tree.cpp gsl_bessel_example.cpp BGL_quickTour.cpp test_belief.cpp property_map_example.cpp vertex_name_property.cpp testMyGraph.cpp
OBJS=$(SRCS:.cpp=.o)
OBJDIR=../obj/test
BINDIR=../bin
DATADIR=../data
SOBJDIR=../obj/src
TSDIR=.
SSDIR=../src
MAKE=make

vpath %.cpp $(TSDIR)
vpath %.cpp $(SSDIR)
vpath %.o $(OBJDIR)
vpath %.o $(SOBJDIR)

all: boost_regex family_tree gsl_bessel BGL_quickTour
	
%.o: %.cpp
	@echo "compiling file $<"
	$(CXX) -c $(CXXFLAGS) $< -o $(OBJDIR)/$@

%.cpp: %.h

BGL_quickTour: BGL_quickTour.o
	$(CXX) $(LDFLAGS) $< -o $(BINDIR)/$@
	@$(BINDIR)/$@

gsl_bessel: gsl_bessel_example.o
	@echo "gsl_bessel test example"
	$(CXX) $(LDFLAGS) $< -o $(BINDIR)/$@
	@$(BINDIR)/$@

boost_regex: boost_regex_example.o
	@echo "boost_regex test example"
	$(CXX) $(LDFLAGS) -lboost_regex $< -o $(BINDIR)/$@
	@$(BINDIR)/$@ < $(DATADIR)/regex_example.txt
	
PM_example: property_map_example.o 
	$(CXX) $(LDFLAGS) $< -o $(BINDIR)/$@
	@$(BINDIR)/$@

family_tree: family_tree.o
	@echo "boostgraph test example"
	$(CXX) $(LDFLAGS) $< -o $(BINDIR)/$@
	@$(BINDIR)/$@
	
test_belief: clean test_belief.o 
	$(CXX) $(LDFLAGS) $(SOBJDIR)/belief.o $(SOBJDIR)/state.o $(SOBJDIR)/bpgraph.o $(SOBJDIR)/node.o $(OBJDIR)/test_belief.o -o $(BINDIR)/$@
	@$(BINDIR)/$@

bptest: bptest.o
	$(CXX) $(LDFLAGS) $(SOBJDIR)/bpgraph.o $(SOBJDIR)/node.o  $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(OBJDIR)/bptest.o -o $(BINDIR)/$@
	@$(BINDIR)/$@
	
vertex_name_property: vertex_name_property.o
	$(CXX) $(LDFLAGS) $(SOBJDIR)/bpgraph.o $(SOBJDIR)/node.o  $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(OBJDIR)/vertex_name_property.o -o $(BINDIR)/$@
	$(BINDIR)/$@
	
myGraph: clean testMyGraph.o
	@$(CXX) $(LDFLAGS) $(SOBJDIR)/myedge.o $(SOBJDIR)/node.o  $(SOBJDIR)/ising.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o $(OBJDIR)/testMyGraph.o -o $(BINDIR)/tester
	@echo "Linking over" 
runMyGraph:	myGraph
	$(BINDIR)/tester
	
video1: video1.o
	@$(CXX) $(LDFLAGS) $(SOBJDIR)/VideoReader.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(OBJDIR)/video1.o -o $(BINDIR)/video1
#	$(BINDIR)/video1 ../data/jumping_girl_input.avi

potts:	pottsMain.o
	@$(CXX) $(LDFLAGS) $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VideoReader.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(OBJDIR)/pottsMain.o $(SOBJDIR)/PottsModel.o -o $(BINDIR)/potts
	
tracker: trackerMain.o 
	$(CXX) $(LDFLAGS) $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/VideoReader.o $(SOBJDIR)/PottsCB.o $(SOBJDIR)/PottsModel.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/Tracker.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(OBJDIR)/trackerMain.o -o $(BINDIR)/tracker
	
clean: 
	@echo "cleaning the $(OBJDIR) directory"
	@rm -rf $(OBJDIR)/*.o
	
optical_flow: optical_flow_demo_linux.full_api.o
	$(CXX) $(LDFLAGS) $(OBJDIR)/$< -o $(BINDIR)/$@
	$(BINDIR)/$@ ../data/jumping_girl_input.avi
		
ptrue:	Porig.o
	$(CXX) $(LDFLAGS) $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VideoReader.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(SOBJDIR)/PottsModel.o $(OBJDIR)/Porig.o -o $(BINDIR)/ptrue	
	
MVI: TestMVI.o
	$(CXX) $(LDFLAGS) $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VR.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/MVI.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(SOBJDIR)/PottsModel.o $(OBJDIR)/TestMVI.o -o $(BINDIR)/MVI	
	@echo "Finished compiling $@"

image_gen: imgGen.o 
	$(CXX) $(LDFLAGS) $< -o $(BINDIR)/image_gen
	@echo "Finished compiling $@"
	
testGraphBP: testGraphBP.o
	$(CXX) $(LDFLAGS)  $(SOBJDIR)/GraphBP.o $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VideoReader.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(SOBJDIR)/PottsModel.o $(OBJDIR)/testGraphBP.o -o $(BINDIR)/testGraphBP
	
testGraphGBP: testGraphGBP.o
	$(CXX) $(LDFLAGS)  $(SOBJDIR)/GraphGBP.o $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VideoReader.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(SOBJDIR)/PottsModel.o $(OBJDIR)/testGraphGBP.o -o $(BINDIR)/testGraphGBP

testAMPL: testAMPL.o
	$(CXX) $(LDFLAGS) $(SOBJDIR)/PottsGraph.o  $(SOBJDIR)/VideoReader.o $(SOBJDIR)/node.o  $(SOBJDIR)/VideoNode.o $(SOBJDIR)/VideoGraph.o $(SOBJDIR)/state.o $(SOBJDIR)/belief.o $(SOBJDIR)/mygraph.o $(SOBJDIR)/main.o  $(SOBJDIR)/GraphPPM.o $(SOBJDIR)/PottsModel.o $(OBJDIR)/testAMPL.o -o $(BINDIR)/testAMPL
	@echo "Finished compiling $@"

