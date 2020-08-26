CC = mpic++
CFLAGS = -Wall -O3 --std=c++11
BFLAGS = -I/home/alessandro/miniconda3/envs/MPI/include -L/home/alessandro/miniconda3/envs/MPI/lib -lboost_system -pthread -lboost_mpi -lboost_serialization

es10_2.exe : es10_2.o random.o
	$(CC) random.o es10_2.o -o es10_2.exe $(BFLAGS)
es10_2.o : es10_2.cpp classi.h
	$(CC) -c es10_2.cpp -o es10_2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es10_2.exe
