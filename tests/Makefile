TARGET=mc-test
CXX=g++ -std=c++17
CXXFLAGS= -Wall -g -MMD -I../include/
LD=g++
OBJS=test.o

.PHONY: all run clean purge

all: $(TARGET)

$(TARGET): $(OBJS)
	$(LD) -o $(TARGET) $(OBJS)

test.o: test.cpp
	$(CXX) -c $(CXXFLAGS) test.cpp -o test.o

run: ${TARGET}
	@./${TARGET}

purge:
	@$(MAKE) clean
	@rm -f ./${TARGET} *~

clean:
	@rm -f *.o *.d


