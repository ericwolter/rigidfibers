#ifndef FIBERS_ASSEMBLE_SYSTEM_1D_KERNEL_
#define FIBERS_ASSEMBLE_SYSTEM_1D_KERNEL_

#include "constants.cu"
#include "compute_inner_integral_analytically.cu"
#include "compute_inner_integral_numerically.cu"

__global__ void
assemble_system_1D(
        #ifdef VALIDATE
        int *validation,
        #endif //VALIDATE
    const float4 *positions,
    const float4 *orientations,
    float *a_matrix,
    float *b_vector
    )
{
  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= NUMBER_OF_FIBERS) return;

  const float c  = logf(SLENDERNESS * SLENDERNESS * M_E);
  const float d  = -c;
  const float e  = 2.0f;
  const float cc = 1.0f;
  const float D1 = 0.75f / (d - 2.0f * cc);

  const float4 position_i = positions[i];
  const float4 orientation_i = orientations[i];

  float4 external_force;
  external_force.x = 0.5f * 0.0f;
  external_force.y = 0.5f * 0.0f;
  external_force.z = 0.5f * -1.0f;

  int x_row_index;
  int y_row_index;
  int z_row_index;

  int x_column_index;
  int y_column_index;
  int z_column_index;

  float G[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 6];
  float GF[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 3];

  for (int j = 0; j < NUMBER_OF_FIBERS; ++j)
  {
    const float4 position_j = positions[j];
    const float4 orientation_j = orientations[j];

    for (int force_index_j = 0; force_index_j < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_j)
    {
      if (i == j)
      {
        x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 0;
        y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 1;
        z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_j + 2;

        x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 0;
        y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 1;
        z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 2;

        a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
        a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = 1;
        a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = 1;

#ifdef VALIDATE
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_j;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 0;
        validation[x_row_index * 6 + x_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 0;

        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_j;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 1;
        validation[y_row_index * 6 + y_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 1;

        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0] = i;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1] = j;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2] = force_index_j;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3] = force_index_j;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4] = 2;
        validation[z_row_index * 6 + z_column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5] = 2;
#endif //VALIDATE
        continue;
      }

      float Q1;
      float Q2;
      float Q3;

      int force_index_i = 0;

      // theta in equation 23
      float T11 = 0.0f;
      float T22 = 0.0f;
      float T33 = 0.0f;
      float T12 = 0.0f;
      float T13 = 0.0f;
      float T23 = 0.0f;

      float TF1 = 0.0f;
      float TF2 = 0.0f;
      float TF3 = 0.0f;
      float QF;

#ifdef USE_ANALYTICAL_INTEGRATION
          compute_G_analytic(position_i, orientation_i, position_j, orientation_j, force_index_j, external_force, G, GF, i == 89 && j == 21);
#else
          compute_G_numeric(position_i, orientation_i, position_j, orientation_j, force_index_j, external_force, G, GF, i == 9 && j == 88 && force_index_j == 1);
#endif //USE_ANALYTICAL_INTEGRATION

      for (int quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
      {
        const float quadrature_weight = quadrature_weights[quadrature_index_i];
        const float legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
        T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;

        if (force_index_j == 0)
        {
          TF1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          TF2 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          TF3 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
        }
      }

      Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
      Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
      Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

      x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 0;
      y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 1;
      z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 2;

      x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 0;
      y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 1;
      z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 2;

      a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.x * Q1;
      a_matrix[x_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.x * Q2;
      a_matrix[x_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.x * Q3;
      a_matrix[y_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.y * Q1;
      a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.y * Q2;
      a_matrix[y_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.y * Q3;
      a_matrix[z_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.z * Q1;
      a_matrix[z_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.z * Q2;
      a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = D1 * orientation_i.z * Q3;

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

      if (force_index_j == 0)
      {
        QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

        b_vector[x_row_index] -= D1 * orientation_i.x * QF;
        b_vector[y_row_index] -= D1 * orientation_i.y * QF;
        b_vector[z_row_index] -= D1 * orientation_i.z * QF;
      }

      for (force_index_i = 1; force_index_i < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index_i)
      {
        const float gamma = 0.5f * (2.0f * (force_index_i + 1) + 1.0f) / (d + e - cc * lambda[force_index_i]);

        T11 = 0.0f;
        T22 = 0.0f;
        T33 = 0.0f;
        T12 = 0.0f;
        T13 = 0.0f;
        T23 = 0.0f;

        TF1 = 0.0f;
        TF2 = 0.0f;
        TF3 = 0.0f;

        for (int quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
        {
          const float quadrature_weight = quadrature_weights[quadrature_index_i];
          const float legendre_polynomial = legendre_polynomials[quadrature_index_i + force_index_i * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
          T11 += quadrature_weight * G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          T22 += quadrature_weight * G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          T33 += quadrature_weight * G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          T12 += quadrature_weight * G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          T13 += quadrature_weight * G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          T23 += quadrature_weight * G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;

          if (force_index_j == 0)
          {
            TF1 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            TF2 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
            TF3 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] * legendre_polynomial;
          }
        }

        Q1 = T11 * orientation_i.x + T12 * orientation_i.y + T13 * orientation_i.z;
        Q2 = T12 * orientation_i.x + T22 * orientation_i.y + T23 * orientation_i.z;
        Q3 = T13 * orientation_i.x + T23 * orientation_i.y + T33 * orientation_i.z;

        x_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 0;
        y_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 1;
        z_row_index = i * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index_i + 2;

        x_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 0;
        y_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 1;
        z_column_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + force_index_j * DIMENSIONS + 2;

        a_matrix[x_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T11 - eigen[force_index_i] * orientation_i.x * Q1);
        a_matrix[x_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T12 - eigen[force_index_i] * orientation_i.x * Q2);
        a_matrix[x_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T13 - eigen[force_index_i] * orientation_i.x * Q3);
        a_matrix[y_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T12 - eigen[force_index_i] * orientation_i.y * Q1);
        a_matrix[y_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T22 - eigen[force_index_i] * orientation_i.y * Q2);
        a_matrix[y_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T23 - eigen[force_index_i] * orientation_i.y * Q3);
        a_matrix[z_row_index + x_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T13 - eigen[force_index_i] * orientation_i.z * Q1);
        a_matrix[z_row_index + y_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T23 - eigen[force_index_i] * orientation_i.z * Q2);
        a_matrix[z_row_index + z_column_index * TOTAL_NUMBER_OF_ROWS] = gamma * (T33 - eigen[force_index_i] * orientation_i.z * Q3);

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

        if (force_index_j == 0)
        {
          QF = TF1 * orientation_i.x + TF2 * orientation_i.y + TF3 * orientation_i.z;

          b_vector[x_row_index] -= gamma * (TF1 - eigen[force_index_i] * orientation_i.x * QF);
          b_vector[y_row_index] -= gamma * (TF2 - eigen[force_index_i] * orientation_i.y * QF);
          b_vector[z_row_index] -= gamma * (TF3 - eigen[force_index_i] * orientation_i.z * QF);
        }
      }
    }
  }
}

#endif // FIBERS_ASSEMBLE_SYSTEM_1D_KERNEL_
