/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = icc
OPT = -O3 -avx2
CFLAGS = -Wall -std=gnu99 $(OPT) -funroll-loops -march=core2 -ftree-vectorize -ffast-math
#MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
#LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
LDLIBS = -lrt  -I$(MKLROOT)/include -Wl,-L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

*/
#include <immintrin.h>


const char* dgemm_desc = "Simple blocked dgemm.";

#if !defined(BLOCK_SIZE)
#define BLOCK_SIZE 64
#endif

#define min(a,b) (((a)<(b))?(a):(b))
/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* restrict C)
{
  /* For each row i of A */
  for (int i = 0; i < M; ++i)
    /* For each column j of B */ 
    for (int j = 0; j < N; ++j) 
    {
      /* Compute C(i,j) */
      double cij = C[i+j*lda];
      __m256d csum[3][4] = {{0}};
      for (int k = 0; k < K; ++k)
        for (int bi = 0; bi < 4; bi++){
        __m256d bb = _mm256_loadu_pd(B+k+bi*8);
          for (int ai = 0; ai < 3; ai++){
        __m256d aa = _mm256_broadcast_sd(A+ai+k);
            csum[ai][bi] += aa * bb;
          }
        }
            // Accumulate the results into C.
        for (int ai = 0; ai < 3; ai++) {
          for (int bi = 0; bi < 4; bi++) {
          cij += _mm256_hadd_pd(C,ai,bi*8,csum[ai][bi]);
           }
          }
       }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format. 
 * On exit, A and B maintain their input values. */  
void square_dgemm (int lda, double* A, double* B, double* restrict C)
{
  /* For each block-row of A */ 
  for (int i = 0; i < lda; i += BLOCK_SIZE)
    /* For each block-column of B */
    for (int j = 0; j < lda; j += BLOCK_SIZE)
      /* Accumulate block dgemms into block of C */
      for (int k = 0; k < lda; k += BLOCK_SIZE)
      {
	/* Correct block dimensions if block "goes off edge of" the matrix */
	int M = min (BLOCK_SIZE, lda-i);
	int N = min (BLOCK_SIZE, lda-j);
	int K = min (BLOCK_SIZE, lda-k);

	/* Perform individual block dgemm */
	do_block(lda, M, N, K, A + i + k*lda, B + k + j*lda, C + i + j*lda);
      }
}
