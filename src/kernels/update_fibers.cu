#ifndef FIBERS_UPDATE_FIBERS_KERNEL_
#define FIBERS_UPDATE_FIBERS_KERNEL_

#include "constants.cu"

__global__
void update_fibers(
    const float4 *previous_positions,
    const float4 *current_positions,
    float4 *next_positions,
    const float4 *previous_orientations,
    const float4 *current_orientations,
    float4 *next_orientations,
    const float4 *previous_translational_velocities,
    const float4 *current_translational_velocities,
    const float4 *previous_rotational_velocities,
    const float4 *current_rotational_velocities
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    next_positions[i].x = 4.0f / 3.0f * current_positions[i].x
                        - 1.0f / 3.0f * previous_positions[i].x
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_translational_velocities[i].x - previous_translational_velocities[i].x);
    next_positions[i].y = 4.0f / 3.0f * current_positions[i].y
                        - 1.0f / 3.0f * previous_positions[i].y
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_translational_velocities[i].y - previous_translational_velocities[i].y);
    next_positions[i].z = 4.0f / 3.0f * current_positions[i].z
                        - 1.0f / 3.0f * previous_positions[i].z
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_translational_velocities[i].z - previous_translational_velocities[i].z);

    next_orientations[i].x = 4.0f / 3.0f * current_orientations[i].x
                        - 1.0f / 3.0f * previous_orientations[i].x
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_rotational_velocities[i].x - previous_rotational_velocities[i].x);
    next_orientations[i].y = 4.0f / 3.0f * current_orientations[i].y
                        - 1.0f / 3.0f * previous_orientations[i].y
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_rotational_velocities[i].y - previous_rotational_velocities[i].y);
    next_orientations[i].z = 4.0f / 3.0f * current_orientations[i].z
                        - 1.0f / 3.0f * previous_orientations[i].z
                        + 2.0f / 3.0f * TIMESTEP * (2.0f * current_rotational_velocities[i].z - previous_rotational_velocities[i].z);

    float invLen = 1.0f / sqrtf(next_orientations[i].x * next_orientations[i].x
        + next_orientations[i].y * next_orientations[i].y 
        + next_orientations[i].z * next_orientations[i].z);

    next_orientations[i].x *= invLen;
    next_orientations[i].y *= invLen;
    next_orientations[i].z *= invLen;
}

#endif //FIBERS_UPDATE_FIBERS_KERNEL_
