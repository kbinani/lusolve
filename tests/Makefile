CXXFLAGS=-Igtest/include -Igtest
CXX=g++
OBJS=lusolve_test.o gtest_main.o gtest-all.o
PROGRAMS=lusolve_test.out

all: $(PROGRAMS)

lusolve_test.out: $(OBJS)

lusolve_test.o: lusolve_test.cpp ../lusolve.hpp

gtest_main.o: gtest/src/gtest_main.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

gtest-all.o: gtest/src/gtest-all.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

.o.out:
	$(CXX) $^ -o $@

clean:
	$(RM) $(PROGRAMS) $(OBJS)
