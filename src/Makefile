PROG = CoNET
CC = g++
CPPFLAGS = -O3 -std=c++14 -Wall -Wextra -lpthread -g
LDFLAGS = -lpthread
OBJS = CoNET.o logger.o parameters.o breakpoints.o csv_reader.o move_type.o

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS)  -o $(PROG) $(OBJS) -lpthread
CoNETo:
	$(CC) $(CPPFLAGS) -c CoNET.cpp
parameters.o: src/parameters/parameters.h
	$(CC) $(CPPFLAGS) -c src/parameters/parameters.cpp 
logger.o: src/utils/logger/logger.h
	$(CC) $(CPPFLAGS) -c src/utils/logger/logger.cpp 
breakpoints.o: src/utils/breakpoints/breakpoints.h
	$(CC) $(CPPFLAGS) -c src/utils/breakpoints/breakpoints.cpp
csv_reader.o: src/cell_provider/csv_reader/csv_reader.h
	$(CC) $(CPPFLAGS) -c src/cell_provider/csv_reader/csv_reader.cpp
move_type.o: src/moves/move_type.h
	$(CC) $(CPPFLAGS) -c src/moves/move_type.cpp
clean:
	rm -f core $(PROG) $(OBJS)
