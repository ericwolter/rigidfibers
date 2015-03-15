#ifndef FIBERS_UPDATE_FIBERS_FIRSTSTEP_KERNEL_
#define FIBERS_UPDATE_FIBERS_FIRSTSTEP_KERNEL_

#include "constants.cu"

__global__
void update_fibers_firststep(
    const float4 *current_positions,
    float4 *next_positions,
    const float4 *current_orientations,
    float4 *next_orientations,
    const float4 *translational_velocities,
    const float4 *rotational_velocities
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    next_positions[i].x = current_positions[i].x + TIMESTEP * translational_velocities[i].x;
    next_positions[i].y = current_positions[i].y + TIMESTEP * translational_velocities[i].y;
    next_positions[i].z = current_positions[i].z + TIMESTEP * translational_velocities[i].z;

    next_orientations[i].x = current_orientations[i].x + TIMESTEP * rotational_velocities[i].x;
    next_orientations[i].y = current_orientations[i].y + TIMESTEP * rotational_velocities[i].y;
    next_orientations[i].z = current_orientations[i].z + TIMESTEP * rotational_velocities[i].z;

    float invLen = 1.0f / sqrtf(next_orientations[i].x * next_orientations[i].x
        + next_orientations[i].y * next_orientations[i].y 
        + next_orientations[i].z * next_orientations[i].z);

    next_orientations[i].x *= invLen;
    next_orientations[i].y *= invLen;
    next_orientations[i].z *= invLen;
}

#endif //FIBERS_UPDATE_FIBERS_FIRSTSTEP_KERNEL_
