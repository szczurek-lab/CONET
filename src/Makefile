BOOST_LIB = /usr/local/lib
BOOST_INC = /usr/local/include

PROG = CONET
CC = g++
CPPFLAGS = -I${BOOST_INC} -O3 -std=c++17 -lpthread -g -pg
LDFLAGS = -lpthread -pg -g -lboost_program_options
OBJS = CONET.o logger.o parameters.o types.o csv_reader.o move_type.o

$(PROG) : $(OBJS)
	$(CC) -I${BOOST_INC} -L${BOOST_LIB} -o $(PROG) $(OBJS) $(LDFLAGS)
CONET.o:
	$(CC) $(CPPFLAGS) -c CONET.cpp
parameters.o: src/parameters/parameters.h
	$(CC) $(CPPFLAGS) -c src/parameters/parameters.cpp
logger.o: src/utils/logger/logger.h
	$(CC) $(CPPFLAGS) -c src/utils/logger/logger.cpp 
types.o: src/types.h
	$(CC) $(CPPFLAGS) -c src/types.cpp
csv_reader.o: src/input_data/csv_reader.h
	$(CC) $(CPPFLAGS) -c src/input_data/csv_reader.cpp
move_type.o: src/moves/move_type.h
	$(CC) $(CPPFLAGS) -c src/moves/move_type.cpp
clean:
	rm -f core $(PROG) $(OBJS)
