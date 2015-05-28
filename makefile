harmonic_cosine: cosine.o
	gcc -pthread -o cosine cosine.o -lgmp

harmonic_cosine.o: cosine.c
	gcc -c cosine.c -Wall -pedantic -ansi -g

clean:
	rm -rf *.o
	rm -rf *~
	rm cosine