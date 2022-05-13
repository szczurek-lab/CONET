PROG = CONET
CC = g++
CPPFLAGS = -O3 -std=c++14 -Wall -Wextra -lpthread -g
LDFLAGS = -lpthread
OBJS = CONET.o logger.o parameters.o types.o csv_reader.o move_type.o

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS)  -o $(PROG) $(OBJS) -lpthread
CONET.o:
	$(CC) $(CPPFLAGS) -c CONET.cpp
parameters.o: src/parameters/parameters.h
	$(CC) $(CPPFLAGS) -c src/parameters/parameters.cpp 
logger.o: src/utils/logger/logger.h
	$(CC) $(CPPFLAGS) -c src/utils/logger/logger.cpp 
types.o: src/types.h
	$(CC) $(CPPFLAGS) -c src/types.cpp
csv_reader.o: src/cell_provider/csv_reader/csv_reader.h
	$(CC) $(CPPFLAGS) -c src/cell_provider/csv_reader/csv_reader.cpp
move_type.o: src/moves/move_type.h
	$(CC) $(CPPFLAGS) -c src/moves/move_type.cpp
clean:
	rm -f core $(PROG) $(OBJS)
