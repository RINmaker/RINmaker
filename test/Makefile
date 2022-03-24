CC = g++
INCLUDES = ./libraries/ ./libraries/gtest/include/ ./libraries/gtest/
CPPFLAGS = --std=c++17 -Wall -Wno-reorder -DNDEBUG -O3 -funroll-loops -finline-functions -frename-registers -pthread $(foreach dir, $(INCLUDES), -I$(dir))

LDFLAGS =
LDLIBS  = -lpthread

SRCS = $(wildcard *.cpp) ./libraries/pugixml/pugixml.cpp 

OBJS      = $(subst .cpp,.o,$(SRCS))
OBJS_TEST = $(subst .cpp,.test.o,$(SRCS))

EXE 	 = RING
EXE_TEST = RING_Test

all: $(EXE) $(EXE_TEST)

$(EXE): $(OBJS) gtest.a
	$(CC) $(LDFLAGS) -o ./out/$(EXE) $(OBJS) gtest.a $(LDLIBS)

$(EXE_TEST): $(OBJS_TEST) gtest.a
	$(CC) $(LDFLAGS) -o ./out/$(EXE_TEST) $(OBJS_TEST) gtest.a $(LDLIBS)

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $< -o $@

%.test.o: %.cpp
	$(CC) $(CPPFLAGS) -D_TEST -c $< -o $@

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CC) $(CPPFLAGS) -MM $^ >>./.depend

clean:
	$(RM) $(OBJS) $(OBJS_TEST)

distclean: clean
	$(RM) $(EXE)
	$(RM) *~ .depend

.PHONY: $(EXE) ring

include .depend

### GTEST ###
GTEST_DIR     = ./libraries/gtest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRCS_   = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

GTEST_CPPFLAGS += -isystem $(GTEST_DIR)/include
GTEST_CXXFLAGS += -g -Wall -Wextra -pthread

gtest-all.o : $(GTEST_SRCS_)
	$(CC) $(GTEST_CPPFLAGS) -I$(GTEST_DIR) $(GTEST_CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^
### END GTEST ###