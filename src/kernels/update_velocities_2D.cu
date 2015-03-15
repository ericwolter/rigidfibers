#ifndef FIBERS_UPDATE_VELOCITIES_2D_KERNEL_
#define FIBERS_UPDATE_VELOCITIES_2D_KERNEL_

#include "constants.cu"

__device__
    void compute_GV_2D(const int j,
                const float4 position_i,
                const float4 orientation_i,
                const float4 position_j,
                const float4 orientation_j,
                const float *coefficients,
                const float4 external_force,
                float *GF
                ) // @TODO better names
{
    for (int quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;

        float4 position_on_fiber_i;
        position_on_fiber_i.x = position_i.x + quadrature_points[quadrature_index_i] * orientation_i.x;
        position_on_fiber_i.y = position_i.y + quadrature_points[quadrature_index_i] * orientation_i.y;
        position_on_fiber_i.z = position_i.z + quadrature_points[quadrature_index_i] * orientation_i.z;

        #pragma unroll
        for (int quadrature_index_j = 0; quadrature_index_j < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_j)
        {
            const float quadrature_point = quadrature_points[quadrature_index_j];
            float4 position_on_fiber_j;
            position_on_fiber_j.x = position_j.x + quadrature_point * orientation_j.x;
            position_on_fiber_j.y = position_j.y + quadrature_point * orientation_j.y;
            position_on_fiber_j.z = position_j.z + quadrature_point * orientation_j.z;

            float4 difference;
            difference.x = position_on_fiber_i.x - position_on_fiber_j.x;
            difference.y = position_on_fiber_i.y - position_on_fiber_j.y;
            difference.z = position_on_fiber_i.z - position_on_fiber_j.z;

            const float d1 = difference.x * difference.x;
            const float d2 = difference.y * difference.y;
            const float d3 = difference.z * difference.z;

            const float invDistance = rsqrtf(d1 + d2 + d3);
            const float invDistance3 = invDistance * invDistance * invDistance;
            const float invDistance5 = invDistance3 * invDistance * invDistance;

            // equation 10
            // Note:    The outer product of a vector with itself is always a symmetric matrix
            //          so to save computation we only compute the upper triangle.
            // TODO calculation can be optimized (i.e. not dividing by distance, simpifing etc.)
            const float K11 = invDistance
                                   + invDistance3 * d1
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * d1);
            const float K22 = invDistance
                                   + invDistance3 * d2
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * d2);
            const float K33 = invDistance
                                   + invDistance3 * d3
                                   + 2.0f * SLENDERNESS * SLENDERNESS * (invDistance3
                                           - 3.0f * invDistance5 * d3);
            const float K12 = invDistance3 * difference.x * difference.y
                                   + 2.0f * SLENDERNESS * SLENDERNESS
                                   * -3.0f * invDistance5 * difference.x * difference.y;
            const float K13 = invDistance3 * difference.x * difference.z
                                   + 2.0f * SLENDERNESS * SLENDERNESS
                                   * -3.0f * invDistance5 * difference.x * difference.z;
            const float K23 = invDistance3 * difference.y * difference.z
                                   + 2.0f * SLENDERNESS * SLENDERNESS
                                   * -3.0f * invDistance5 * difference.y * difference.z;

            const float quadrature_weight = quadrature_weights[quadrature_index_j];

            float4 force_on_fiber_j;
            force_on_fiber_j.x = 0.5f * external_force.x;
            force_on_fiber_j.y = 0.5f * external_force.y;
            force_on_fiber_j.z = 0.5f * external_force.z;

            for (int force_index = 0; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
            {
                const float legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
                int x_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 0;
                int y_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 1;
                int z_row_index = j * NUMBER_OF_TERMS_IN_FORCE_EXPANSION * DIMENSIONS + DIMENSIONS * force_index + 2;

                force_on_fiber_j.x += coefficients[x_row_index] * legendre_polynomial;
                force_on_fiber_j.y += coefficients[y_row_index] * legendre_polynomial;
                force_on_fiber_j.z += coefficients[z_row_index] * legendre_polynomial;
            }

            GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K11 * force_on_fiber_j.x + K12 * force_on_fiber_j.y + K13 * force_on_fiber_j.z);
            GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K12 * force_on_fiber_j.x + K22 * force_on_fiber_j.y + K23 * force_on_fiber_j.z);
            GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K13 * force_on_fiber_j.x + K23 * force_on_fiber_j.y + K33 * force_on_fiber_j.z);
        }
    }
}

__global__ void update_velocities_2D(
    const float4 *positions,
    const float4 *orientations,
    const float *coefficients,
    const float4 *external_forces,
    float4 *translational_velocities,
    float4 *rotational_velocities
)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= NUMBER_OF_FIBERS) return;
    if (j >= NUMBER_OF_FIBERS) return;
    if (i==j) return;

    const float c  = logf(SLENDERNESS * SLENDERNESS * M_E);
    const float d  = -c;

    const float4 position_i = positions[i];
    const float4 orientation_i = orientations[i];

    float4 external_force = external_forces[i];

    const float4 position_j = positions[j];
    const float4 orientation_j = orientations[j];

    float GF[TOTAL_NUMBER_OF_QUADRATURE_POINTS * 3];
    compute_GV_2D(j, position_i, orientation_i, position_j, orientation_j, coefficients, external_force, GF);

    float TF1A0 = 0.0f;
    float TF2A0 = 0.0f;
    float TF3A0 = 0.0f;

    float TF1A1 = 0.0f;
    float TF2A1 = 0.0f;
    float TF3A1 = 0.0f;

    for (int quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
      const float quadrature_weight = quadrature_weights[quadrature_index_i];
      const float legendre_polynomial = legendre_polynomials[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
      const float weighted_polynomial = quadrature_weight * legendre_polynomial;

      TF1A0 += quadrature_weight * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
      TF2A0 += quadrature_weight * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
      TF3A0 += quadrature_weight * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

      TF1A1 += weighted_polynomial * GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
      TF2A1 += weighted_polynomial * GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
      TF3A1 += weighted_polynomial * GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
    }

    const float t1 = orientation_i.x * TF1A1;
    const float t2 = orientation_i.y * TF2A1;
    const float t3 = orientation_i.z * TF3A1;

    atomicAdd(&(translational_velocities[i].x), (0.5f / d) * TF1A0);
    atomicAdd(&(translational_velocities[i].y), (0.5f / d) * TF2A0);
    atomicAdd(&(translational_velocities[i].z), (0.5f / d) * TF3A0);

    atomicAdd(&(rotational_velocities[i].x), (1.5f / d) * (TF1A1 - (orientation_i.x * t1 + orientation_i.x * t2 + orientation_i.x * t3)));
    atomicAdd(&(rotational_velocities[i].y), (1.5f / d) * (TF2A1 - (orientation_i.y * t1 + orientation_i.y * t2 + orientation_i.y * t3)));
    atomicAdd(&(rotational_velocities[i].z), (1.5f / d) * (TF3A1 - (orientation_i.z * t1 + orientation_i.z * t2 + orientation_i.z * t3)));
}

#endif //FIBERS_UPDATE_VELOCITIES_2D_KERNEL_
