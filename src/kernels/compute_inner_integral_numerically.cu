#ifndef FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_
#define FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_

#include "constants.cu"

__device__
void compute_G_numeric(
     const float4 position_i,
     const float4 orientation_i,
     const float4 position_j,
     const float4 orientation_j,
     const int force_index,
     const float4 external_force,
     float *G,
     float *GF,
     const bool debug)
{
    for (int quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        // float CGF0;
        // float CGF1;
        // float CGF2;
        // float YGF0;
        // float YGF1;
        // float YGF2;
        // float TGF0;
        // float TGF1;
        // float TGF2;

        G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;
        G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = 0.0f;

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
            const float legendre_polynomial = legendre_polynomials[quadrature_index_j + force_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS];

//            if(debug && quadrature_index_j==0) {
//                printf("%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t\n",
//                       K11,
//                       K22,
//                       K33,
//                       quadrature_weight,
//                       legendre_polynomial,
//                       K23
////                       2.0f * SLENDERNESS * SLENDERNESS
////                                                          * -3.0f * invDistance5 * difference.x * difference.y,
////                       2.0f * SLENDERNESS * SLENDERNESS
////                                                          * -3.0f * invDistance5 * difference.x * difference.z,
////                       2.0f * SLENDERNESS * SLENDERNESS
////                                                          * -3.0f * invDistance5 * difference.y * difference.z
//                       );
//            }

            // @TEST http://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html#1262
            // Kahan Summation Formula

            // if (quadrature_index_j == 0) {
            //     GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z);
            //     GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z);
            //     GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z);

            //     CGF0 = 0.0f;
            //     CGF1 = 0.0f;
            //     CGF2 = 0.0f;
            // } else {
            //     YGF0 = quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z) - CGF0;
            //     YGF1 = quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z) - CGF1;
            //     YGF2 = quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z) - CGF2;

            //     TGF0 = GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] + YGF0;
            //     TGF1 = GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] + YGF1;
            //     TGF2 = GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] + YGF2;

            //     CGF0 = (TGF0 - GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]) -  YGF0;
            //     CGF1 = (TGF1 - GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]) -  YGF1;
            //     CGF2 = (TGF2 - GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]) -  YGF2;

            //     GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = TGF0;
            //     GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = TGF1;
            //     GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = TGF2;
            // }

            G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K11 * legendre_polynomial;
            G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K22 * legendre_polynomial;
            G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K33 * legendre_polynomial;
            G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K12 * legendre_polynomial;
            G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K13 * legendre_polynomial;
            G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * K23 * legendre_polynomial;

            GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K11 * external_force.x + K12 * external_force.y + K13 * external_force.z);
            GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K12 * external_force.x + K22 * external_force.y + K23 * external_force.z);
            GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] += quadrature_weight * (K13 * external_force.x + K23 * external_force.y + K33 * external_force.z);
        }

//        if(debug) {
//            printf("%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t\n",
//                   G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
//                   G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
//                   G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
//                   G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
//                   G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
//                   G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS]
//                   );
//        }
    }
}

#endif //FIBERS_COMPUTE_INNER_INTEGRAL_NUMERICALLY_KERNEL_
