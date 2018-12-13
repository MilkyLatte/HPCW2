stencil: finalStencil.c
	mpicc -xHost -std=c99 -Ofast -Wall -qopt-report=1 -qopt-report-phase=vec $^ -o $@

