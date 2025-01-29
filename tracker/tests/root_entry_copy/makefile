CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra
ROOTFLAGS = $(shell root-config --cflags --libs)

CopyTreeBranches: test_copy.cc
	$(CXX) $(CXXFLAGS) -g -o CopyTreeBranches test_copy.cc $(ROOTFLAGS)

clean:
	rm -f CopyTreeBranches