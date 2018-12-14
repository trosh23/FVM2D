### Begin MACROS
#Define compiler
CC=g++
#Define linker
LD=g++
#Compilation flags
CFLAGS=−c
#Define sources files to compile
SOURCES=main.cpp
#Define generated object files from SOURCES
OBJECTS=$(SOURCES:.cpp=.o)
#Define name executable
EXEC = FVM2D.exe
### End MACROS


### Begin targets
build: $(SOURCES) $(EXEC)
$(EXEC) : $(OBJECTS)
	$(LD) $(OBJECTS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) $< −o $@

#Define additional targets
clean:
	rm -rf $(OBJECTS)
all : clean build

### End targets