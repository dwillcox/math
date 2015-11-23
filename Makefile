all:
	gfortran -Og -c polynomial.f90
	gfortran -Og -c quadtester.f90
	gfortran -o tester polynomial.o quadtester.o

clean:
	rm *.o *.mod tester
