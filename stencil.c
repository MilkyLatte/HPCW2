#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>

#define MASTER 0

int main(int argc, char* argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  MPI_Init( &argc, &argv );

  int size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  int nx = atoi(argv[1]);
  int ny = atoi(argv[2])/size;
  int niters = atoi(argv[3]);


  // Allocate the image
  float *image = malloc(sizeof(float)*nx*ny);
  float *tmp_image = malloc(sizeof(float)*nx*ny);

  int rank;
  int flag;
  enum bool {FALSE, TRUE};
  char hostname[MPI_MAX_PROCESSOR_NAME];

  int dest;
  int source;
  int tag = 0;
  MPI_Status status;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  MPI_Finalize();

  return EXIT_SUCCESS;
}

void init_image(const int nx, const int ny, float * restrict image, float * restrict  tmp_image) {
  // Zero everything
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      image[j+i*ny] = 0.0f;
      tmp_image[j+i*ny] = 0.0f;
    }
  }

  // Checkerboard
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int jj = j*ny*0.125f; jj < (j+1)*ny*0.125f; ++jj) {
        for (int ii = i*nx*0.125f; ii < (i+1)*nx*0.125f; ++ii) {
          if ((i+j)%2) image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}
