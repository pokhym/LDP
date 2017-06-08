EXENAME = test

COMPILER = g++
WARNINGS = -Wall
COMPILER_OPTS = -c -g $(WARNINGS) -std=c++11
LINKER = g++

MAIN_OBJS = main.o util.o data.o
MAIN_DEPS = main.cpp util.h data.h

DATA_DEPS = data.h
UTIL_DEPS = util.h data.h

CODE_CLN = *.o $(EXENAME)

all: $(EXENAME)

$(EXENAME): main.o data.o util.o
	$(LINKER) -g -Wall -std=c++11 -o $(EXENAME) $(MAIN_OBJS)

main.o: $(MAIN_DEPS)
	$(COMPILER) $(COMPILER_OPTS) main.cpp

util.o: $(UTIL_DEPS)
	$(COMPILER) $(COMPILER_OPTS) util.cpp

data.o: $(DATA_DEPS)
	$(COMPILER) $(COMPILER_OPTS) data.cpp

clean:
	-rm -f $(CODE_CLN)
