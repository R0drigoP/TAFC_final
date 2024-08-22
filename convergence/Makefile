CC = g++
CFLAGS = -c -Wall
LDFLAGS = 

SOURCES = main.cpp molecule.cpp 
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main
TILDE := $(wildcard */*~) $(wildcard *~)
TXT := $(wildcard */*.txt) $(wildcard *.txt)

all: 
	$(CC) main.cpp molecule.cpp `root-config --cflags --glibs` -o main

parallel:
	$(CC) -fopenmp main.cpp molecule.cpp `root-config --cflags --glibs` -o main

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(TILDE) $(TXT)
