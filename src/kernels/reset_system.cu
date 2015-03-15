#ifndef FIBERS_REEST_MATRIX_KERNEL_
#define FIBERS_RESET_MATRIX_KERNEL_

#include "../common.h"
#include "constants.cu"

__global__ void reset_system(
    #ifdef VALIDATE
    int *validation,
    #endif //VALIDATE
    float *a_matrix,
    float *b_vector
)
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= NUMBER_OF_FIBERS) return;

  int x_row_index;
  int y_row_index;
  int z_row_index;

  int x_column_index;
  int y_column_index;
  int z_column_index;

  for (int force_index_i = 0; force_index_i < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_i)
  {
    x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 0;
    y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 1;
    z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 2;

    b_vector[x_row_index] = 0;
    b_vector[y_row_index] = 0;
    b_vector[z_row_index] = 0;

    for (int j = 0; j < NUMBER_OF_FIBERS; ++j)
    {
      for (int force_index_j = 0; force_index_j < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_j)
      {
        x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 0;
        y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 1;
        z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 2;

        if(i==j && force_index_i == force_index_j) {
          a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
          a_matrix[x_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[x_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[y_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
          a_matrix[y_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[z_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[z_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
        } else {
          a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[x_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[x_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[y_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[y_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[z_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[z_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
          a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 0;
        }
#ifdef VALIDATE
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 0;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 0;

        validation[x_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[x_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[x_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[x_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[x_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 0;
        validation[x_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 1;

        validation[x_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[x_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[x_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[x_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[x_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 0;
        validation[x_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 2;

        validation[y_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[y_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[y_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[y_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[y_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 1;
        validation[y_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 0;

        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 1;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 1;

        validation[y_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[y_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[y_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[y_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[y_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 1;
        validation[y_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 2;

        validation[z_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[z_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[z_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[z_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[z_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 2;
        validation[z_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 0;

        validation[z_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[z_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[z_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[z_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[z_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 2;
        validation[z_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 1;

        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_i;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 2;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 2;
#endif //VALIDATE
      }
    }
  }
}

#endif //FIBERS_RESET_MATRIX_KERNEL_
