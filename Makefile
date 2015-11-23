all:
	gfortran -Og -c polynomial.f90 -ffpe-trap=invalid,zero,overflow,underflow
	gfortran -Og -c quadtester.f90 -ffpe-trap=invalid,zero,overflow,underflow
	gfortran -o tester polynomial.o quadtester.o

clean:
	rm *.o *.mod tester
