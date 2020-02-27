PHONY.:all

all: apr czebyszew aprox intrp prosta

apr: main.o splines.o points.o aproksymator.o gaus/libge.a
	$(CC) -o apr main.o splines.o points.o aproksymator.o -L gaus -l ge -lm

czebyszew: main.o splines.o points.o aproksymator_baza_czebysz.o gaus/libge.a
	$(CC) -o czebyszew main.o splines.o points.o aproksymator_baza_czebysz.o -L gaus -l ge -lm

aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

aproksymator_baza_czebysz.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_baza_czebysz.c
aproksymator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator.c

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta

