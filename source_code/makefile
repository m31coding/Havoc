CC = gcc
CXX = g++

DIRS = $(shell find . -type d)

CFLAGS = -Wall -g -DGIT_COMMIT='$(GIT_COMMIT)'
CXXFLAGS = $(CFLAGS)

SRCS_C = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.c))
SRCS_CXX = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.cpp))
SRCS = $(SRCS_C) $(SRCS_CXX)
OBJS_C = $(SRCS_C:%.c=%.o)
OBJS_CXX = $(SRCS_CXX:%.cpp=%.o)
OBJS = $(OBJS_C) $(OBJS_CXX)
HDRS = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.h))
TEMP = $(SRCS:%=%~) $(HDRS:%=%~)

LIBS = -lm

EXENAME = havoc 

$(EXENAME): $(OBJS)
	$(CXX) $(OBJS) $(CXXFLAGS) -o $(EXENAME) $(LIBS)

$(OBJS): $(HDRS) makefile

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXENAME) $(TEMP) constants~

%.o: %.c
	$(CC) $*.c -c $(CFLAGS) -o $*.o

%.o: %.cpp
	$(CXX) $*.cpp -c $(CXXFLAGS) -o $*.o

#current git commit
GIT_COMMIT = "$(shell git rev-parse HEAD)"