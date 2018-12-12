#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>
#include <sys/time.h>

#define MASTER 0
#define OUTPUT_FILE "stencil.pgm"


void init_image(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image, int rank);
void stencil(const int nx, const int ny, float * restrict image, float * restrict  tmp_image, int rank, int size);
void topStencil(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image, int rank, int size);
void botStencil(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image, int rank, int size);

void output_image(const char * file_name, const int nx, const int ny, float * restrict image);
double wtime(void);

int calc_nrows_from_rank(int rank, int size, int ny);
void haloExchange(const int ncols, const int nrows, int up, int down, float * image, float * sendbuf, float * recvbuf, MPI_Status status);
void sendInitial(const int local_ncols, const int local_nrows, float * restrict image, int size);
void joinImage(const int local_nrows, const int local_ncols, float * restrict fullImage, float * restrict image, int size, int rank, MPI_Status status);

int main(int argc, char* argv[]){
  if (argc != 4) {
    fprintf(stderr, "Usage: %s nx ny niters\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int size;
  int rank;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );


  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int niters = atoi(argv[3]);

  int up;
  int down;
  int flag;
  MPI_Status status;
  int local_nrows;
  int local_ncols;
  float *fullImage;
  float *full_tmp_image;
  float *image;
  float *tmp_image;
  float *sbuf;
  float *rbuf;


  up = (rank == MASTER) ? (rank + size - 1) : (rank - 1);
  down = (rank + 1) % size;

  local_nrows = calc_nrows_from_rank(rank, size, nx);
  local_ncols = ny;

  sbuf = (float*)malloc(sizeof(float)*local_ncols);
  rbuf = (float*)malloc(sizeof(float)*local_ncols);




  if (rank == MASTER) {
      fullImage = malloc(sizeof(float)*nx*ny);
      full_tmp_image = malloc(sizeof(float)*nx*ny);
      image = malloc(sizeof(float) * (local_nrows+2) * local_ncols);
      tmp_image = malloc(sizeof(float) * (local_nrows+2) * local_ncols);
      init_image(nx, ny, fullImage, full_tmp_image, rank);
      sendInitial(local_ncols, local_nrows, fullImage, size);
      for (int i = 0; i < local_ncols*(local_nrows+2); i++){
        if (i < local_ncols){
          image[i] = rank;
        } else if (i > (local_ncols * (local_nrows+1)-1)){
          image[i] = rank;
        } else {
          image[i] = fullImage[i - local_ncols];
        }
      }
  } else {
    float *recvbuf = malloc(sizeof(float) * local_nrows * local_ncols);
    image = malloc(sizeof(float) * (local_nrows+2) * local_ncols);
    tmp_image = malloc(sizeof(float) * (local_nrows+2) * local_ncols);

    MPI_Recv(recvbuf, local_nrows*local_ncols, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, &status);
    for (int i = 0; i < local_ncols*(local_nrows+2); i++){
      if (i < local_ncols){
        image[i] = rank;
      } else if (i > (local_ncols * (local_nrows+1)-1)){
        image[i] = rank;
      } else {
        image[i] = recvbuf[i - local_ncols];
      }
    }
  }
  for (int i = 0; i < local_nrows + 2; ++i) {
    for (int j = 0; j < local_ncols; ++j) {
      tmp_image[j+i*ny] = 0.0f;
    }
  }
 haloExchange(local_ncols, local_nrows, up, down, image, sbuf, rbuf, status);
//  haloExchange(local_ncols, local_nrows, up, down, image, sbuf, rbuf, status);

 double tic = wtime();

 for (int t = 0; t < niters; ++t) {
   if (rank == 0) {
     topStencil(local_nrows+2, local_ncols, image, tmp_image, rank, size);
   } else if (rank == size - 1){
     botStencil(local_nrows+2, local_ncols, image, tmp_image, rank, size);
   } else {
     stencil(local_nrows+2, local_ncols, image, tmp_image, rank, size);
   }
   haloExchange(local_ncols, local_nrows, up, down, tmp_image, sbuf, rbuf, status);
   if (rank == 0) {
     topStencil(local_nrows+2, local_ncols, tmp_image, image, rank, size);
   } else if (rank == size - 1){
     botStencil(local_nrows+2, local_ncols, tmp_image, image, rank, size);
   } else {
     stencil(local_nrows+2, local_ncols, tmp_image, image, rank, size);
   }
   haloExchange(local_ncols, local_nrows, up, down, image, sbuf, rbuf, status);
 }
 double toc = wtime();

//  if (rank == 0) {
//    for (int j = 0; j < local_nrows+2; j++) {
//      for (int i = 0; i < local_ncols; i++) {
//        printf(" %f |", tmp_image[j*(local_ncols)+i]);
//        //printf(" %d |", j*(local_ncols+2)+i);
//      }
//      printf("\n");
//    }
// }

// Output
 printf("------------------------------------\n");
 printf(" runtime: %lf s\n", toc-tic);
 printf("------------------------------------\n");
 joinImage(local_nrows, local_ncols, fullImage, image, size, rank, status);

 if (rank == 0) {
   // for (int j = 0; j < nx; j++) {
   //   for (int i = 0; i < ny; i++) {
   //     printf(" %f |", fullImage[j*(ny)+i]);
   //   }
   //   printf("\n");
   // }
   output_image(OUTPUT_FILE, nx, ny, fullImage);
 }




  MPI_Finalize();
}

void sendInitial(const int local_ncols, const int local_nrows, float * restrict image, int size){
  float *sendbuf;
  sendbuf = (float*)malloc(sizeof(float)*local_ncols*local_nrows);
  for (int k = 1; k < size; k++) {
    for (int i = 0; i < local_ncols*local_nrows; i++){
      sendbuf[i] = image[i + k*local_ncols*local_nrows];
    }
    MPI_Send(sendbuf, local_ncols*local_nrows, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
  }
  free(sendbuf);
}

void init_image(const int nx, const int ny, float * restrict image, float * restrict  tmp_image, int rank) {
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
          if ((i+j)%2)

	image[jj+ii*ny] = 100.0f;
        }
      }
    }
  }
}

void haloExchange(const int ncols, const int nrows, int up, int down, float * image, float * sendbuf, float * recvbuf, MPI_Status status) {
  int tag = 0;

  for (int i = 0; i < ncols; i++){
    sendbuf[i] = image[1 * (ncols) + i];
  }
  MPI_Sendrecv(sendbuf, ncols, MPI_FLOAT, up, tag, recvbuf, ncols, MPI_FLOAT, down, tag, MPI_COMM_WORLD, &status);
  for (int i = 0; i < ncols; i++){
    image[(nrows+1) * ncols + i] = recvbuf[i];
  }

  for(int i = 0; i < ncols; i++){
    sendbuf[i] = image[(nrows) * ncols + i];
  }
  MPI_Sendrecv(sendbuf, ncols, MPI_FLOAT, down, tag, recvbuf, ncols, MPI_FLOAT, up, tag, MPI_COMM_WORLD, &status);
  for(int i = 0; i < ncols; i++){
    image[i] = recvbuf[i];
  }
}

void topStencil(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image, int rank, int size){
  float x = 0.6f;
  float y = 0.1f;
  for (int i = 1; i < nx-1; ++i) {
    for (int j = 1; j < ny-1; ++j) {
      float variable = image[j+i*ny] * x;
      if (i != 1) variable += image[j  +(i-1)*ny] * y;
      variable += image[j  +(i+1)*ny] * y;
      variable += image[j-1+i*ny] * y;
      variable += image[j+1+i*ny] * y;
      tmp_image[j+i*ny] = variable;
    }
  }

  for (int j = 1; j < ny-1; ++j) {
    float variable = image[j+(nx-1)*ny] * x;
    variable += image[j +(nx-2)*ny] * y;
    variable += image[j-1+(nx-1)*ny] * y;
    variable += image[j+1+(nx-1)*ny] * y;
    tmp_image[j+(nx-1)*ny] = variable;
  }

  for (int i = 1; i < nx-1; ++i){
   float variable = image[i*ny] * x;
   if (i != 1) variable += image[(i-1)*ny] * y;
   variable += image[(i+1)*ny] * y;
   variable += image[1+i*ny] * y;
   tmp_image[i*ny] = variable;
  }

  for (int i = 1; i < nx-1; ++i) {
    float variable = image[(ny-1)+i*ny]*x;
    if (i != 1) variable += image[(ny-1)+(i-1)*ny] * y;
    variable += image[(ny-1)+(i+1)*ny] * y;
    variable += image[(ny-1)-1+i*ny] * y;
    tmp_image[(ny-1)+i*ny] = variable;
  }

  float leftTop = image[(nx-1) * ny] * x;
  leftTop += image[((nx-1)-1) * ny] * y;
  leftTop += image[1+(nx-1)*ny] * y;
  tmp_image[(nx-1)*ny] = leftTop;
  float rightTop = image[(ny-1)+(nx-1) * ny] * x;
  rightTop += image[(ny-1) + (nx-1-1) * ny] * y;
  rightTop += image[((ny-1) - 1) + (nx-1) * ny] * y;
  tmp_image[(ny-1)+(nx-1) * ny] = rightTop;
}
void botStencil(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image, int rank, int size){
  float x = 0.6f;
  float y = 0.1f;
  for (int i = 1; i < nx-1; ++i) {
    for (int j = 1; j < ny-1; ++j) {
      float variable = image[j+i*ny] * x;
      variable += image[j  +(i-1)*ny] * y;
      if (i != nx-2) variable += image[j  +(i+1)*ny] * y;
      variable += image[j-1+i*ny] * y;
      variable += image[j+1+i*ny] * y;
      tmp_image[j+i*ny] = variable;
    }
  }

  for (int j = 1; j < ny-1; ++j) {
    float variable = image[j] * x;
    variable += image[j + nx] * y;
    variable += image[j-1] * y;
    variable += image[j+1] * y;
    tmp_image[j] = variable;
  }

  for (int i = 1; i < nx-1; ++i){
   float variable = image[i*ny] * x;
   variable += image[(i-1)*ny] * y;
   if (i != nx-2) variable += image[(i+1)*ny] * y;
   variable += image[1+i*ny] * y;
   tmp_image[i*ny] = variable;
  }

  for (int i = 1; i < nx-1; ++i) {
    float variable = image[(ny-1)+i*ny]*x;
    variable += image[(ny-1)+(i-1)*ny] * y;
    if (i != nx-2) variable += image[(ny-1)+(i+1)*ny] * y;
    variable += image[(ny-1)-1+i*ny] * y;
    tmp_image[(ny-1)+i*ny] = variable;
  }

  float leftBot = image[0] * x;
  leftBot += image[ny] * x;
  leftBot += image[1] * x;
  tmp_image[0] = leftBot;
  float rightBot = image[ny-1] * x;
  rightBot += image[(ny-1)+ny] * y;
  rightBot += image[(ny-1)-1] * y;
  tmp_image[ny-1] = rightBot;

}
void stencil(const int nx, const int ny, float * restrict  image, float * restrict  tmp_image, int rank, int size) {
  float x = 0.6f;
  float y = 0.1f;
  for (int i = 1; i < nx-1; ++i) {
    for (int j = 1; j < ny-1; ++j) {
      float variable = image[j+i*ny] * x;
      variable += image[j  +(i-1)*ny] * y;
      variable += image[j  +(i+1)*ny] * y;
      variable += image[j-1+i*ny] * y;
      variable += image[j+1+i*ny] * y;
      tmp_image[j+i*ny] = variable;
    }
  }
  // when i = 0;
  for (int j = 1; j < ny-1; ++j) {
    float variable = image[j] * x;
    variable += image[j + nx] * y;
    variable += image[j-1] * y;
    variable += image[j+1] * y;
    tmp_image[j] = variable;
  }

  // when i == nx - 1
  for (int j = 1; j < ny-1; ++j) {
    float variable = image[j+(nx-1)*ny] * x;
    variable += image[j +(nx-2)*ny] * y;
    variable += image[j-1+(nx-1)*ny] * y;
    variable += image[j+1+(nx-1)*ny] * y;
    tmp_image[j+(nx-1)*ny] = variable;
  }



  // when j == 0
  for (int i = 1; i < nx-1; ++i){
   float variable = image[i*ny] * x;
   variable += image[(i-1)*ny] * y;
   variable += image[(i+1)*ny] * y;
   variable += image[1+i*ny] * y;
   tmp_image[i*ny] = variable;
  }
  // when j == nx -1
  for (int i = 1; i < nx-1; ++i) {
    float variable = image[(ny-1)+i*ny]*x;
    variable += image[(ny-1)+(i-1)*ny] * y;
    variable += image[(ny-1)+(i+1)*ny] * y;
    variable += image[(ny-1)-1+i*ny] * y;
    tmp_image[(ny-1)+i*ny] = variable;
  }

  if (rank != 10){
    float leftBot = image[0] * x;
    leftBot += image[ny] * x;
    leftBot += image[1] * x;
    tmp_image[0] = leftBot;
    float rightBot = image[ny-1] * x;
    rightBot += image[(ny-1)+ny] * y;
    rightBot += image[(ny-1)-1] * y;
    tmp_image[ny-1] = rightBot;
  }

  if (rank != 10){
    float leftTop = image[(nx-1) * ny] * x;
    leftTop += image[((nx-1)-1) * ny] * y;
    leftTop += image[1+(nx-1)*ny] * y;
    tmp_image[(nx-1)*ny] = leftTop;
    float rightTop = image[(ny-1)+(nx-1) * ny] * x;
    rightTop += image[(ny-1) + (nx-1-1) * ny] * y;
    rightTop += image[((ny-1) - 1) + (nx-1) * ny] * y;
    tmp_image[(ny-1)+(nx-1) * ny] = rightTop;
  }

}

void joinImage(const int local_nrows, const int local_ncols, float * restrict fullImage, float * restrict image, int size, int rank, MPI_Status status){
  float * buf = malloc(sizeof(float) * local_nrows * local_ncols);
  if (rank == 0) {
    for(int i = 0; i < local_ncols * local_nrows; i++) {
      fullImage[i] = image[i + local_ncols];
    }
    for (int k = 1; k < size; k++) {
      MPI_Recv(buf, local_ncols*local_nrows, MPI_FLOAT, k, 0, MPI_COMM_WORLD, &status);
      for(int i = 0; i < local_ncols * local_nrows; i++) {
        fullImage[i + local_ncols*local_nrows*k] = buf[i];
      }
    }
  } else {
    for(int i = 0; i < local_ncols * local_nrows; i++) {
      buf[i] = image[i + local_ncols];
    }
    MPI_Send(buf, local_ncols*local_nrows, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
  }
}

void output_image(const char * file_name, const int nx, const int ny, float * restrict image) {

  // Open output file
  FILE *fp = fopen(file_name, "w");
  if (!fp) {
    fprintf(stderr, "Error: Could not open %s\n", OUTPUT_FILE);
    exit(EXIT_FAILURE);
  }

  // Ouptut image header
  fprintf(fp, "P5 %d %d 255\n", nx, ny);

  // Calculate maximum value of image
  // This is used to rescale the values
  // to a range of 0-255 for output
  float maximum = 0.0f;
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      if (image[j+i*ny] > maximum)
        maximum = image[j+i*ny];
    }
  }

  // Output image, converting to numbers 0-255
  for (int i = 0; i < ny; ++i) {
    for (int j = 0; j < nx; ++j) {
      fputc((char)(255*image[j+i*ny]/maximum), fp);
    }
  }

  // Close the file
  fclose(fp);

}

double wtime(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec*1e-6;
}

int calc_nrows_from_rank(int rank, int size, const int nx)
{
  int nrows;

  nrows = nx / size;       /* integer division */
  if ((nx % size) != 0) {  /* if there is a remainder */
    if (rank == size - 1)
      nrows += nx % size;  /* add remainder to last rank */
  }

  return nrows;
}
