EXENAME = test

COMPILER = g++
WARNINGS = -Wall
COMPILER_OPTS = -c -g -O0 -Werror $(WARNINGS)
LINKER = g++

MAIN_OBJS = main.o data.o util.o
MAIN_DEPS = data.h data.cpp data.o util.h util.cpp util.o main.cpp

DATA_DEPS = data.h
UTIL_DEPS = util.h

CODE_CLN = *.o $(EXENAME)

all: $(EXENAME)

$(EXENAME): main.o
	$(LINKER) $(MAIN_OBJS) -o $(EXENAME)

main.o: $(MAIN_DEPS)
	$(COMPILER) $(COMPILER_OPTS) main.cpp

util.o: $(UTIL_DEPS)
	$(COMPILER) $(COMPILER_OPTS) util.cpp

data.o: $(DATA_DEPS)
	$(COMPILER) $(COMPILER_OPTS) data.cpp

clean:
	-rm -f $(CODE_CLN)
