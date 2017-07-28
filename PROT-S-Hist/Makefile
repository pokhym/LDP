EXENAME = test

COMPILER = g++
CC = cc
WARNINGS = -Wall
CXXFLAGS = -c -g $(WARNINGS) -std=c++11
CXXFLAGS       += -I./c++ -std=c++11								\
		    -Wall -Wextra -pedantic -Wno-missing-braces -Wwrite-strings			\
		    -Wpointer-arith -Wnon-virtual-dtor -Woverloaded-virtual			\
		    -Wsign-promo -Wswitch -Wreturn-type
CXXFLAGS       += -O3
CFLAGS		= # -O3 already defined
LINKER = g++

MAIN_OBJS = main.o util.o data.o djelic_bch.o
MAIN_DEPS = main.cpp util.h data.h djelic_bch.h

DATA_DEPS = data.h
UTIL_DEPS = util.h data.h

CODE_CLN = *.o $(EXENAME)

all: $(EXENAME)

$(EXENAME): main.o data.o util.o  djelic_bch.o
	$(LINKER) -g -Wall -std=c++11 -o $(EXENAME) $(MAIN_OBJS)

main.o: $(MAIN_DEPS)
	$(COMPILER) $(CXXFLAGS) main.cpp

util.o: $(UTIL_DEPS)
	$(COMPILER) $(CXXFLAGS) util.cpp

data.o: $(DATA_DEPS)
	$(COMPILER) $(CXXFLAGS) data.cpp

djelic:
	git clone https://github.com/pjkundert/bch.git $@

c++/ezpwd/bch \
djelic/include \
djelic/lib/bch.c: djelic

.PHONY: djelictest
djelictest:	djelic/Documentation/bch/nat_tu_tool

djelic/Documentation/bch/nat_tu_tool: djelic
	cd djelic/Documentation/bch && make && ./nat_tu_short.sh

djelic_bch.c:	CFLAGS += -I standalone -I djelic/Documentation/bch/standalone -I djelic/include
djelic_bch.o:	CFLAGS += -I standalone -I djelic/Documentation/bch/standalone -I djelic/include
djelic_bch.o:	djelic_bch.c		djelic/lib/bch.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	-rm -f $(CODE_CLN)

clean:
	-rm -f $(CODE_CLN)
