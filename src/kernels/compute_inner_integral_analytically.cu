#ifndef FIBERS_COMPUTE_INNER_INTEGRAL_ANALYTICALLY_KERNEL_
#define FIBERS_COMPUTE_INNER_INTEGRAL_ANALYTICALLY_KERNEL_

#include "constants.cu"

#define pown2(x) (x*x)
#define pown3(x) (x*x*x)
#define pown4(x) (x*x*x*x)

__device__
void compute_G_analytic(
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
    const int k = force_index + 1;

    float I1[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];
    float I3[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];
    float I5[NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3];

    for (int quadrature_index_i = 0; quadrature_index_i < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++quadrature_index_i)
    {
        float4 position_on_fiber_i;
        position_on_fiber_i.x = position_i.x + quadrature_points[quadrature_index_i] * orientation_i.x;
        position_on_fiber_i.y = position_i.y + quadrature_points[quadrature_index_i] * orientation_i.y;
        position_on_fiber_i.z = position_i.z + quadrature_points[quadrature_index_i] * orientation_i.z;

        // the difference vector between the current point on the fiber and the center
        // of the other fiber
        float4 R0;
        R0.x = position_on_fiber_i.x - position_j.x;
        R0.y = position_on_fiber_i.y - position_j.y;
        R0.z = position_on_fiber_i.z - position_j.z;
        const float b = -2.0f
                             * (R0.x * orientation_j.x
                                + R0.y * orientation_j.y
                                + R0.z * orientation_j.z);
        const float c = R0.x * R0.x
                             + R0.y * R0.y
                             + R0.z * R0.z;
        const float invC = 1.0f / c;
        //const float c_test = ((R0.x * R0.y * R0.z) * (R0.x * R0.y * R0.z)) - 2 * R0.x * R0.y - 2 * R0.x * R0.z - 2 * R0.y * R0.z;

        // if fibers are too far apart we have numerical problems
        // so in order to minimize the effect we inverse the
        // the recursive direction
        // @TODO How/why does this help exactly?
        const float climit = 10.0f;

        const float d = c - 0.25f * b * b;

        const float s_upper = 1.0f;
        const float s_lower = -1.0f;

        const float u_upper = sqrtf(s_upper * s_upper + b * s_upper + c);
        const float u_lower = sqrtf(s_lower * s_lower + b * s_lower + c);

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     float y = 0;
        //     const float x = fabs(2.0f * s_upper + b + 2.0f * u_upper);

        //     for (int i = 0; i < 100; ++i)
        //     {
        //         y = y + 2 * ((x - exp(y))/(x + exp(y)));
        //         printf("%d, %f\n",i,y);
        //     }
        // }

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf(" %.16f\n", 0, position_j.x);
        //     printf(" %.16f\n", 0, c);
        //     printf(" %.16f\n", 0, c_test);
        //     printf(" %.16f\n", 0, u_upper);
        //     printf(" %.16f\n", 0, 2.0f * s_upper + b + 2.0f * u_upper);
        //     printf(" %.16f\n", 0, fabs(2.0f * s_upper + b + 2.0f * u_upper));
        //     printf(" %.16f\n", 0, log(fabs(2.0f * s_upper + b + 2.0f * u_upper)));
        // }

        I1[0] = logf(fabs(2.0f * s_upper + b + 2.0f * u_upper)) - logf(fabs(2.0f * s_lower + b + 2.0f * u_lower));
        I1[1] = u_upper - u_lower + (-b / 2.0f) * I1[0];

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("%d  %.16f\n", 0, I1[0]);
        //     printf("%d  %.16f \n", 1, I1[1]);
        // }

        I3[0] = (d < 1e-7f) ?
                (-2.0f / pown2(2.0f * s_upper + b)) - (-2.0f / pown2(2.0f * s_lower + b)) :
                ((2.0f * s_upper + b) / (2.0f * d * u_upper)) - ((2.0f * s_lower + b) / (2.0f * d * u_lower));
        I3[1] = (-1.0f / u_upper) - (-1.0f / u_lower) - b / 2.0f * I3[0];

        I5[0] = (d < 1e-7f) ?
                (-4.0f / pown4(2.0f * s_upper + b)) - (-4 / pown4(2.0f * s_lower + b)) :
                ((2.0f * s_upper + b) / (6.0f * d * pown3(u_upper))) - ((2.0f * s_lower + b) / (6.0f * d * pown3(u_lower))) + (2.0f / (3.0f * d)) * I3[0];
        I5[1] = (d < 1e-7f) ?
                (-8.0f / (3.0f * pown3(2.0f * s_upper + b))) - (-8.0f / (3.0f * pown3(2.0f * s_lower + b))) - (b / 2.0f) * I5[0] :
                (-(b * s_upper + 2.0f * c) / (6.0f * d * pown3(u_upper))) - (-(b * s_lower + 2.0f * c) / (6.0f * d * pown3(u_lower))) - (b / (3.0f * d)) * I3[0];

        if (c < climit)
        {
            float pow_s_upper = s_upper;
            float pow_s_lower = s_lower;

            for (int n = 2; n < k + 3; ++n)
            {
                I1[n] = (pow_s_upper * u_upper) / n - (pow_s_lower * u_lower) / n
                        + ((1.0f - 2.0f * n) * b) / (2.0f * n) * I1[n - 1] - ((n - 1) * c) / n * I1[n - 2];

                // if (debug && quadrature_index_i == 12 && k == 5)
                // {
                //     printf("%d  %.16f \n", n, I1[n]);
                // }

                I3[n] = I1[n - 2] - b * I3[n - 1] - c * I3[n - 2];

                I5[n] = I3[n - 2] - b * I5[n - 1] - c * I5[n - 2];
                pow_s_upper *= s_upper;
                pow_s_lower *= s_lower;
            }
        }
        else
        {
            float i1n2 = 0.0f;
            float i1n1 = 0.0f;
            float i1n0 = 0.0f;

            float i3n2 = 0.0f;
            float i3n1 = 0.0f;
            float i3n0 = 0.0f;

            float i5n2 = 0.0f;
            float i5n1 = 0.0f;
            float i5n0 = 0.0f;

            // if (debug && quadrature_index_i == 0 && k == 5)
            // {
            //     printf("%d,%f\n", 29, i1n2);
            //     printf("%d,%f\n", 28, i1n1);
            // }
            for (int n = 27; n > 1; --n)
            {
                i1n0 = (n + 2.0f) / ((n + 1.0f) * c) * (
                           ((powf(s_upper, n + 1) * u_upper) / (n + 2.0f)) - ((powf(s_lower, n + 1) * u_lower) / (n + 2.0f))
                           + ((1.0f - 2.0f * (n + 2.0f)) / (2.0f * (n + 2.0f))) * b * i1n1 - i1n2);
                // i1n0 = ((n + 2.0f) / (n + 1.0f)) * invC * (
                //             (__powf(s_upper, n+1) * u_upper - __powf(s_lower, n+1) * u_lower) / (n+2.0f)
                //             + (1.0f/(2.0f*(n+2)) - 1.0f) * b * i1n1 - i1n2
                //         );
                // i1n0 =
                //     (1.0f/(n+1)) * invC * (
                //             (__powf(s_upper, n+1) * u_upper - __powf(s_lower, n+1) * u_lower)
                //             + (-1.5f-n) * b *i1n1 - i1n2 * (n+2)
                //         );

                // if(fabs(test - i1n0) > 1e-7) {
                //     printf("asdfasdf\n");
                // }

                i3n0 = invC * (i1n0 - b * i3n1 - i3n2);
                i5n0 = invC * (i3n0 - b * i5n1 - i5n2);

                // if (debug && quadrature_index_i == 0 && k == 5)
                // {
                //     printf("%d,%f,%f,%f,%f,%f\n", n, i1n0, i1n1, i1n2,
                //            ((powfn(s_upper, n + 1) * u_upper) / (n + 2)) - ((powfn(s_lower, n + 1) * u_lower) / (n + 2)),
                //            ((1.0f - 2.0f * (n + 2.0f)) / (2.0f * (n + 2.0f))) );
                // }

                if (n < NUMBER_OF_TERMS_IN_FORCE_EXPANSION + 3)
                {
                    I1[n] = i1n0;
                    I3[n] = i3n0;
                    I5[n] = i5n0;
                }

                i1n2 = i1n1;
                i1n1 = i1n0;

                i3n2 = i3n1;
                i3n1 = i3n0;

                i5n2 = i5n1;
                i5n1 = i5n0;
            }
        }

        float L01;
        float L03;
        float L05;
        float L13;
        float L15;
        float L23;
        float L25;

        if (k == 0 || k == 1)
        {
            L01 = I1[k];
            L03 = I3[k];
            L05 = I5[k];
            L13 = I3[k + 1];
            L15 = I5[k + 1];
            L23 = I3[k + 2];
            L25 = I5[k + 2];
        }
        else if (k == 2)
        {
            L01 = 0.5f * (3.0f * I1[k] - I1[k - 2]);
            L03 = 0.5f * (3.0f * I3[k] - I3[k - 2]);
            L05 = 0.5f * (3.0f * I5[k] - I5[k - 2]);
            L13 = 0.5f * (3.0f * I3[k + 1] - I3[k - 1]);
            L15 = 0.5f * (3.0f * I5[k + 1] - I5[k - 1]);
            L23 = 0.5f * (3.0f * I3[k + 2] - I3[k]);
            L25 = 0.5f * (3.0f * I5[k + 2] - I5[k]);
        }
        else if (k == 3)
        {
            L01 = 0.5f * (5.0f * I1[k] - 3.0f * I1[k - 2]);
            L03 = 0.5f * (5.0f * I3[k] - 3.0f * I3[k - 2]);
            L05 = 0.5f * (5.0f * I5[k] - 3.0f * I5[k - 2]);
            L13 = 0.5f * (5.0f * I3[k + 1] - 3.0f * I3[k - 1]);
            L15 = 0.5f * (5.0f * I5[k + 1] - 3.0f * I5[k - 1]);
            L23 = 0.5f * (5.0f * I3[k + 2] - 3.0f * I3[k]);
            L25 = 0.5f * (5.0f * I5[k + 2] - 3.0f * I5[k]);
        }
        else if (k == 4)
        {
            L01 = 0.125f * (35.0f * I1[k] - 30.0f * I1[k - 2] + 3.0f * I1[k - 4]);
            L03 = 0.125f * (35.0f * I3[k] - 30.0f * I3[k - 2] + 3.0f * I3[k - 4]);
            L05 = 0.125f * (35.0f * I5[k] - 30.0f * I5[k - 2] + 3.0f * I5[k - 4]);
            L13 = 0.125f * (35.0f * I3[k + 1] - 30.0f * I3[k - 1] + 3.0f * I3[k - 3]);
            L15 = 0.125f * (35.0f * I5[k + 1] - 30.0f * I5[k - 1] + 3.0f * I5[k - 3]);
            L23 = 0.125f * (35.0f * I3[k + 2] - 30.0f * I3[k] + 3.0f * I3[k - 2]);
            L25 = 0.125f * (35.0f * I5[k + 2] - 30.0f * I5[k] + 3.0f * I5[k - 2]);
        }
        else if (k == 5)
        {
            L01 = 0.125f * (63.0f * I1[k] - 70.0f * I1[k - 2] + 15.0f * I1[k - 4]);
            L03 = 0.125f * (63.0f * I3[k] - 70.0f * I3[k - 2] + 15.0f * I3[k - 4]);
            L05 = 0.125f * (63.0f * I5[k] - 70.0f * I5[k - 2] + 15.0f * I5[k - 4]);
            L13 = 0.125f * (63.0f * I3[k + 1] - 70.0f * I3[k - 1] + 15.0f * I3[k - 3]);
            L15 = 0.125f * (63.0f * I5[k + 1] - 70.0f * I5[k - 1] + 15.0f * I5[k - 3]);
            L23 = 0.125f * (63.0f * I3[k + 2] - 70.0f * I3[k] + 15.0f * I3[k - 2]);
            L25 = 0.125f * (63.0f * I5[k + 2] - 70.0f * I5[k] + 15.0f * I5[k - 2]);
        }
        else if (k == 6)
        {
            L01 = 0.0625f * (231.0f * I1[k] - 315.0f * I1[k - 2] + 105.0f * I1[k - 4] - 5.0f * I1[k - 6]);
            L03 = 0.0625f * (231.0f * I3[k] - 315.0f * I3[k - 2] + 105.0f * I3[k - 4] - 5.0f * I3[k - 6]);
            L05 = 0.0625f * (231.0f * I5[k] - 315.0f * I5[k - 2] + 105.0f * I5[k - 4] - 5.0f * I5[k - 6]);
            L13 = 0.0625f * (231.0f * I3[k + 1] - 315.0f * I3[k - 1] + 105.0f * I3[k - 3] - 5.0f * I3[k - 5]);
            L15 = 0.0625f * (231.0f * I5[k + 1] - 315.0f * I5[k - 1] + 105.0f * I5[k - 3] - 5.0f * I5[k - 5]);
            L23 = 0.0625f * (231.0f * I3[k + 2] - 315.0f * I3[k] + 105.0f * I3[k - 2] - 5.0f * I3[k - 4]);
            L25 = 0.0625f * (231.0f * I5[k + 2] - 315.0f * I5[k] + 105.0f * I5[k - 2] - 5.0f * I5[k - 4]);
        }
        else if (k == 7)
        {
            L01 = 0.0625f * (429.0f * I1[k] - 693.0f * I1[k - 2] + 315.0f * I1[k - 4] - 35.0f * I1[k - 6]);
            L03 = 0.0625f * (429.0f * I3[k] - 693.0f * I3[k - 2] + 315.0f * I3[k - 4] - 35.0f * I3[k - 6]);
            L05 = 0.0625f * (429.0f * I5[k] - 693.0f * I5[k - 2] + 315.0f * I5[k - 4] - 35.0f * I5[k - 6]);
            L13 = 0.0625f * (429.0f * I3[k + 1] - 693.0f * I3[k - 1] + 315.0f * I3[k - 3] - 35.0f * I3[k - 5]);
            L15 = 0.0625f * (429.0f * I5[k + 1] - 693.0f * I5[k - 1] + 315.0f * I5[k - 3] - 35.0f * I5[k - 5]);
            L23 = 0.0625f * (429.0f * I3[k + 2] - 693.0f * I3[k] + 315.0f * I3[k - 2] - 35.0f * I3[k - 4]);
            L25 = 0.0625f * (429.0f * I5[k + 2] - 693.0f * I5[k] + 315.0f * I5[k - 2] - 35.0f * I5[k - 4]);
        }
        else if (k == 8)
        {
            L01 = 0.0078125f * (6435.0f * I1[k] - 12012.0f * I1[k - 2] + 6930.0f * I1[k - 4] - 1260.0f * I1[k - 6] + 35.0f * I1[k - 8]);
            L03 = 0.0078125f * (6435.0f * I3[k] - 12012.0f * I3[k - 2] + 6930.0f * I3[k - 4] - 1260.0f * I3[k - 6] + 35.0f * I3[k - 8]);
            L05 = 0.0078125f * (6435.0f * I5[k] - 12012.0f * I5[k - 2] + 6930.0f * I5[k - 4] - 1260.0f * I5[k - 6] + 35.0f * I5[k - 8]);
            L13 = 0.0078125f * (6435.0f * I3[k + 1] - 12012.0f * I3[k - 1] + 6930.0f * I3[k - 3] - 1260.0f * I3[k - 5] + 35.0f * I3[k - 7]);
            L15 = 0.0078125f * (6435.0f * I5[k + 1] - 12012.0f * I5[k - 1] + 6930.0f * I5[k - 3] - 1260.0f * I5[k - 5] + 35.0f * I5[k - 7]);
            L23 = 0.0078125f * (6435.0f * I3[k + 2] - 12012.0f * I3[k] + 6930.0f * I3[k - 2] - 1260.0f * I3[k - 4] + 35.0f * I3[k - 6]);
            L25 = 0.0078125f * (6435.0f * I5[k + 2] - 12012.0f * I5[k] + 6930.0f * I5[k - 2] - 1260.0f * I5[k - 4] + 35.0f * I5[k - 6]);
        }

        G[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G11
            L01
            + L03 * (R0.x * R0.x)
            - L13 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
            + L23 * (orientation_j.x * orientation_j.x)
            + 2.0f * SLENDERNESS * SLENDERNESS * (
                L03
                - 3.0f * L05 * (R0.x * R0.x)
                + 3.0f * L15 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
                - 3.0f * L25 * (orientation_j.x * orientation_j.x)
            );
        G[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G22
            L01
            + L03 * (R0.y * R0.y)
            - L13 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
            + L23 * (orientation_j.y * orientation_j.y)
            + 2.0f * SLENDERNESS * SLENDERNESS * (
                L03
                - 3.0f * L05 * (R0.y * R0.y)
                + 3.0f * L15 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
                - 3.0f * L25 * (orientation_j.y * orientation_j.y)
            );
        G[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G33
            L01
            + L03 * (R0.z * R0.z)
            - L13 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
            + L23 * (orientation_j.z * orientation_j.z)
            + 2.0f * SLENDERNESS * SLENDERNESS * (
                L03
                - 3.0f * L05 * (R0.z * R0.z)
                + 3.0f * L15 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
                - 3.0f * L25 * (orientation_j.z * orientation_j.z)
            );
        G[quadrature_index_i + 3 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G12
            L03 * (R0.x * R0.y)
            - L13 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
            + L23 * (orientation_j.x * orientation_j.y)
            + 2.0f * SLENDERNESS * SLENDERNESS * (
                - 3.0f * L05 * (R0.x * R0.y)
                + 3.0f * L15 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
                - 3.0f * L25 * (orientation_j.x * orientation_j.y)
            );
        G[quadrature_index_i + 4 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G13
            L03 * (R0.x * R0.z)
            - L13 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
            + L23 * (orientation_j.x * orientation_j.z)
            + 2.0f * SLENDERNESS * SLENDERNESS * (
                - 3.0f * L05 * (R0.x * R0.z)
                + 3.0f * L15 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
                - 3.0f * L25 * (orientation_j.x * orientation_j.z)
            );
        G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = //G23
            L03 * (R0.y * R0.z)
            - L13 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
            + L23 * (orientation_j.y * orientation_j.z)
            + 2.0f * SLENDERNESS * SLENDERNESS * (
                - 3.0f * L05 * (R0.y * R0.z)
                + 3.0f * L15 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
                - 3.0f * L25 * (orientation_j.y * orientation_j.z)
            );

        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("a:%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
        //             c,
        //             log(fabs(2.0f * s_upper + b + 2.0f * u_upper)),
        //            I1[0],
        //            I1[1],
        //            I1[2],
        //            I1[3],
        //            I1[4],
        //            I1[5],
        //            I1[6],
        //            I1[7]
        //           );
        // }
        // if (debug && quadrature_index_i == 12 && k == 5)
        // {
        //     printf("b:%f,%f,%f,%f,%f,%f,%f,%f\n",
        //         G[quadrature_index_i + 5 * TOTAL_NUMBER_OF_QUADRATURE_POINTS],
        //            L01,
        //            L03,
        //            L05,
        //            L13,
        //            L15,
        //            L23,
        //            L25
        //           );
        // }

        // if (debug && quadrature_index_i==0 && force_index == 5) {
        //     printf("%f,%f,%f\n",
        //         orientation_j.x * R0.x + R0.x * orientation_j.x,
        //         L13,
        //         (orientation_j.x * R0.x + R0.x * orientation_j.x)*L13);
        // }

        if (force_index == 0)
        {
            L01 = I1[0];
            L03 = I3[0];
            L05 = I5[0];
            L13 = I3[0 + 1];
            L15 = I5[0 + 1];
            L23 = I3[0 + 2];
            L25 = I5[0 + 2];

            const float G11 =
                L01
                + L03 * (R0.x * R0.x)
                - L13 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
                + L23 * (orientation_j.x * orientation_j.x)
                + 2.0f * SLENDERNESS * SLENDERNESS * (
                    L03
                    - 3.0f * L05 * (R0.x * R0.x)
                    + 3.0f * L15 * (orientation_j.x * R0.x + R0.x * orientation_j.x)
                    - 3.0f * L25 * (orientation_j.x * orientation_j.x)
                );
            const float G22 =
                L01
                + L03 * (R0.y * R0.y)
                - L13 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
                + L23 * (orientation_j.y * orientation_j.y)
                + 2.0f * SLENDERNESS * SLENDERNESS * (
                    L03
                    - 3.0f * L05 * (R0.y * R0.y)
                    + 3.0f * L15 * (orientation_j.y * R0.y + R0.y * orientation_j.y)
                    - 3.0f * L25 * (orientation_j.y * orientation_j.y)
                );
            const float G33 =
                L01
                + L03 * (R0.z * R0.z)
                - L13 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
                + L23 * (orientation_j.z * orientation_j.z)
                + 2.0f * SLENDERNESS * SLENDERNESS * (
                    L03
                    - 3.0f * L05 * (R0.z * R0.z)
                    + 3.0f * L15 * (orientation_j.z * R0.z + R0.z * orientation_j.z)
                    - 3.0f * L25 * (orientation_j.z * orientation_j.z)
                );
            const float G12 =
                L03 * (R0.x * R0.y)
                - L13 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
                + L23 * (orientation_j.x * orientation_j.y)
                + 2.0f * SLENDERNESS * SLENDERNESS * (
                    - 3.0f * L05 * (R0.x * R0.y)
                    + 3.0f * L15 * (orientation_j.x * R0.y + R0.x * orientation_j.y)
                    - 3.0f * L25 * (orientation_j.x * orientation_j.y)
                );
            const float G13 =
                L03 * (R0.x * R0.z)
                - L13 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
                + L23 * (orientation_j.x * orientation_j.z)
                + 2.0f * SLENDERNESS * SLENDERNESS * (
                    - 3.0f * L05 * (R0.x * R0.z)
                    + 3.0f * L15 * (orientation_j.x * R0.z + R0.x * orientation_j.z)
                    - 3.0f * L25 * (orientation_j.x * orientation_j.z)
                );
            const float G23 =
                L03 * (R0.y * R0.z)
                - L13 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
                + L23 * (orientation_j.y * orientation_j.z)
                + 2.0f * SLENDERNESS * SLENDERNESS * (
                    - 3.0f * L05 * (R0.y * R0.z)
                    + 3.0f * L15 * (orientation_j.y * R0.z + R0.y * orientation_j.z)
                    - 3.0f * L25 * (orientation_j.y * orientation_j.z)
                );

            GF[quadrature_index_i + 0 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] =
                G11 * external_force.x + G12 * external_force.y + G13 * external_force.z;
            GF[quadrature_index_i + 1 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] =
                G12 * external_force.x + G22 * external_force.y + G23 * external_force.z;
            GF[quadrature_index_i + 2 * TOTAL_NUMBER_OF_QUADRATURE_POINTS] =
                G13 * external_force.x + G23 * external_force.y + G33 * external_force.z;
        }
    }
}

#endif //FIBERS_COMPUTE_INNER_INTEGRAL_ANALYTICALLY_KERNEL_
