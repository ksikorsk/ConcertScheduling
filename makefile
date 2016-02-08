CXX := g++
CXXFLAGS := -g # -Wall

BIN=./bin/
SOURCE=./source/

all: clean build

clean:
	rm -f -r $(BIN)

bin:
	mkdir $(BIN)

HEADERS=$(notdir $(shell find $(SOURCE) -maxdepth 1 -name "*.hh"))
SRCS=$(patsubst %.hh, %, $(HEADERS))
OBJS=$(addprefix $(BIN), $(addsuffix .o, $(SRCS)))
build: bin $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(BIN)ConcertScheduling $(OBJS) $(SOURCE)ConcertScheduling.cc

$(SRCS): $(SOURCE)$(addsuffix .hh, $@) $(SOURCE)$(addsuffix .cc, $@)
	$(CXX) $(CXXFLAGS) -o $(BIN)$@.o -c $(SOURCE)$@.cc	

run:
	$(BIN)ConcertScheduling