CC 	= g++

ERASE = rm -f
CFLAGS = -g  -O3 -ffast-math -Wall -c -I${CFIO_INC_DIR} -I ../
LFLAGS = -g -O3 -ffast-math -Wall  -lm -I${CFIO_INC_DIR} -Wl,-rpath,${CFIODIR} -L${CFIODIR} -lcfitsio -I ../ -L ../



OBJECTS = random.o readtext.o galaxycomplex.o cgparams.o sphericalcell.o catalogobject.o starobject.o galaxyobject.o

all: phosimcatgen

phosimcatgen: phosimcatgen.o $(OBJECTS)
	 $(CC) -o phosimcatgen  phosimcatgen.o $(OBJECTS) $(LFLAGS)

phosimcatgen.o: phosimcatgen.cpp
	$(CC) $(CFLAGS)  phosimcatgen.cpp -o  phosimcatgen.o

cgparams.o: cgparams.cpp
	$(CC) $(CFLAGS) cgparams.cpp  -o cgparams.o

catalogobject.o: catalogobject.cpp
	$(CC) $(CFLAGS) catalogobject.cpp  -o catalogobject.o

starobject.o: starobject.cpp
	$(CC) $(CFLAGS) starobject.cpp  -o starobject.o

galaxyobject.o: galaxyobject.cpp
	$(CC) $(CFLAGS) galaxyobject.cpp  -o galaxyobject.o

random.o: ../ancillary/random.cpp
	$(CC) $(CFLAGS) ../ancillary/random.cpp -o random.o

readtext.o: ../ancillary/readtext.h ../ancillary/readtext.cpp
	$(CC) $(CFLAGS) ../ancillary/readtext.cpp -o readtext.o

galaxycomplex.o: galaxycomplex.cpp
	$(CC) $(CFLAGS) galaxycomplex.cpp  -o galaxycomplex.o

sphericalcell.o: sphericalcell.cpp
	$(CC) $(CFLAGS) sphericalcell.cpp  -o sphericalcell.o

clean:
	$(ERASE) $(OBJECTS) phosimcatgen.o
