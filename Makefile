
CFLAGS	:= $(shell root-config --cflags)
LIBS    := $(shell root-config --glibs) -lMinuit
INCLUDE := -I$(ROOTSYS)/include -I./

TARGET	= test_hypo

SRCS=hypotheses.cc mclimit_csm.cc
OBJS=$(subst .cc,.o,$(SRCS))

all: $(TARGET)

%.o: %.cc
	g++ -O2 $(CFLAGS) $(INCLUDE) -c -o $@ $<

$(TARGET): $(OBJS)
	g++ -O2 $(CFLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

clean:
	rm -f $(TARGET) *.o *~ *.eps *.pdf *.histogram

