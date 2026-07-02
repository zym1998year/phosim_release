CC = ${CCP}
CFLAGS = -g -Wall -std=c++14 -c -I${CFIO_INC_DIR} -I${FFTW3_INC_DIR} -I ../
LFLAGS = -g -Wall -I${CFIO_INC_DIR} -L${CFIODIR} -Wl,-rpath,${CFIODIR} -lcfitsio -I${FFTW3_INC_DIR} -L${FFTW3_DIR} -lfftw3 -lm -lz
OBJECTS =   measure.o

measure: $(OBJECTS)
	 $(CC) $(OBJECTS) $(LFLAGS) -o ./measure

measure.o: measure.cpp
	$(CC) $(CFLAGS) measure.cpp -o measure.o

all: measure

# g++ -I${CFIO_INC_DIR} -L${CFIODIR} -lcfitsio -lm measure.cpp -o measure"
