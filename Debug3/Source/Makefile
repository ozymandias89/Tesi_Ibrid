CXX = g++
CXXFLAGS = -g -Wall -O
LDADD =

UPDIR_BIN = ../bin/
CPX_INCDIR  = /opt/CPLEX_Studio/cplex/cplex/include
CPX_LDPATH  = /opt/CPLEX_Studio/cplex/cplex/lib/x86-64_linux/static_pic
CPX_LDADD = -lcplex -lm -pthread

OBJ = main.o
EXECS = main

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(CPX_INCDIR) -c $^ -o $@

main: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(EXECS) -L$(CPX_LDPATH) $(CPX_LDADD)
	mv $(OBJ) $(UPDIR_BIN)
	mv $(EXECS) $(UPDIR_BIN)

clean:
	rm -rf $(UPDIR_BIN)$(OBJ) $(UPDIR_BIN)$(EXECS)

.PHONY: clean