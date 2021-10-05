/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 * */

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "likwid-stuff.h"

const char* dgemm_desc = "Blocked dgemm, OpenMP-enabled";

void copy_matrix_block(double **S, double **D, int brl, int bcl, int bs)
{
  //std::cout << "inside copy - row location is: " << brl << "column location is: " << bcl << "\n";
  for (int row = 0; row < bs; row++)
  {
     for (int col = 0; col < bs; col++)
     {
        //std::cout << "before copy - row is: " << row+brl << "; column is: " << col+bcl << "\n";
        D[row][col] = S[row + brl][col + bcl];
        //std::cout << "after copy - D value is: " << D[row][col] << "; S value is: " << S[row+brl][col+bcl] << "\n";
     }
  }
}

void matrix_multiply(double **AS, double **BS, double **PROD, int num_rows, int num_cols)
{
   for (int row = 0; row < num_rows; row++)
   {
      for (int col = 0; col < num_cols; col++)
      {
         for (int k = 0; k < num_cols; k++)
         {
            //std::cout << "multiplying AAA to BBB - for: [" << row << "][" << col << "]\n";
            PROD[row][col] += AS[row][k] * BS[k][col];
         }
      }
   }
}

void copy_block_to_matrix(double **S, double **D, int brl, int bcl, int bs)
{
  for (int row = 0; row < bs; row++)
  {
     for (int col = 0; col < bs; col++)
     {
        std::cout << "cpoying block at: [" << brl+row << "][" << bcl+col << "]\n";
        D[brl+row][bcl+col] = S[row][col];
     }
  }
}

void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   // insert your code here: implementation of blocked matrix multiply with copy optimization and OpenMP parallelism enabled

   // be sure to include LIKWID_MARKER_START(MY_MARKER_REGION_NAME) inside the block of parallel code,
   // but before your matrix multiply code, and then include LIKWID_MARKER_STOP(MY_MARKER_REGION_NAME)
   // after the matrix multiply code but before the end of the parallel code block.

  // declare and dynamically allocate 2D arrays
  double **AA, **BB, **CC;

  AA = new double *[n];
  BB = new double *[n];
  CC = new double *[n];
  for (int i = 0; i < n; i++)
  {
     AA[i] = new double [n];
     BB[i] = new double [n];
     CC[i] = new double [n];
  }
  
  // copy column major vector A, B, and C into 2D arrays AA, BB, and CC respectively 
  for (int k = 0, i = 0; i < n*n; k++, i+=n)
  {
     for (int j = 0; j < n; j++)
     {
        AA[j][k] = A[i+j];
        BB[j][k] = B[i+j];
        CC[j][k] = C[i+j];
     }
  }
  
  // block matrix multiplication logic

  #pragma omp parallel
  {
      double **AAA, **BBB, **CCC;  // matrix block arrays
      // allocate memory for block matrix copy
      AAA = new double *[block_size];
      BBB = new double *[block_size];
      CCC = new double *[block_size];
      for (int i = 0; i < block_size; i++)
      {
         AAA[i] = new double [block_size];
         BBB[i] = new double [block_size];
         CCC[i] = new double [block_size];
      }
//#ifdef LIKWID_PERFMON
      LIKWID_MARKER_START(MY_MARKER_REGION_NAME);
//#endif
      #pragma omp for
      for (int ii = 0; ii < n; ii += block_size)  // partition rows by block size; iterate for n/block_size blocks
      {
        for (int jj = 0; jj < n; jj += block_size) // partition columns by block size; iterate for n/block_size blocks
        {
          // std::cout << "copy of product matrix block\n";
          copy_matrix_block(CC, CCC, ii, jj, block_size);
          for (int kk = 0; kk < n; kk += block_size)  // for each row and column of blocks
          {
            //std::cout << "copying AA matrix block to AAA - stating at: [" << ii << "][" << kk << "]\n";
            copy_matrix_block(AA, AAA, ii, kk, block_size);
            //std::cout << "copying BB matrix block to BBB - stating at: [" << kk << "][" << jj << "]\n";
            copy_matrix_block(BB, BBB, kk, jj, block_size);
            // basic matrix multiple applied to matrix blocks
            //std::cout << "multiplying AAA and BBB - product in CCC \n";
            matrix_multiply(AAA, BBB, CCC, block_size, block_size);
            //std::cout << "at this point BMMCO is done for a block\n";
            //std::cout << "ii is: " << ii << "; jj is: " << jj << "; kk is: " << kk << "\n";
          }
          // copy block product to produc matrix
          std::cout << "copy product block CCC to product matric CC for: " << ii << jj << "\n";
          copy_block_to_matrix(CCC, CC, ii, jj, block_size);
        }
      } //end #pragma omp for
      for (int i = 0; i < block_size; i++)
      {
         delete [] AAA[i];
         delete [] BBB[i];
         delete [] CCC[i];
      }
      delete [] AAA;
      delete [] BBB;
      delete [] CCC;
//#ifdef LIKWID_PERFMON
      LIKWID_MARKER_STOP(MY_MARKER_REGION_NAME);
//#endif
  } // end #pragma omp parallel
  
  // copy 2d array CC to column major vector C
  for (int i = 0; i < n; i++)
  {
     for (int j = 0; j < n; j++)
     {
        C[i*n+j] = CC[j][i];
     }
  }

  // release allocated memory
  for (int i = 0; i < n; i++)
  {
     delete [] AA[i];
     delete [] BB[i];
     delete [] CC[i];
  }
  delete [] AA;
  delete [] BB;
  delete [] CC;
}
