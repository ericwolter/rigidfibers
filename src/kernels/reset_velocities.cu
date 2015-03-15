#ifndef FIBERS_RESET_VELOCITIES_KERNEL_
#define FIBERS_RESET_VELOCITIES_KERNEL_

#include "constants.cu"

__global__ void reset_velocities(
    const float4 *orientations,
    float4 *translational_velocities,
    float4 *rotational_velocities
)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= NUMBER_OF_FIBERS) return;

    const float c  = logf(SLENDERNESS * SLENDERNESS * M_E);
    const float d  = -c;

    const float4 orientation_i = orientations[i];

    // @TODO Constant external force
    float4 external_force;
    external_force.x = 0.0f;
    external_force.y = 0.0f;
    external_force.z = -1.0f;

    float4 oriented_force;
    oriented_force.x = orientation_i.x * orientation_i.x * external_force.x + orientation_i.x * orientation_i.y * external_force.y + orientation_i.x * orientation_i.z * external_force.z;
    oriented_force.y = orientation_i.x * orientation_i.y * external_force.x + orientation_i.y * orientation_i.y * external_force.y + orientation_i.y * orientation_i.z * external_force.z;
    oriented_force.z = orientation_i.x * orientation_i.z * external_force.x + orientation_i.y * orientation_i.z * external_force.y + orientation_i.z * orientation_i.z * external_force.z;

    translational_velocities[i].x = (0.5f / d) * ((d + 2.0f) * external_force.x + (d - 2.0f) * oriented_force.x);
    translational_velocities[i].y = (0.5f / d) * ((d + 2.0f) * external_force.y + (d - 2.0f) * oriented_force.y);
    translational_velocities[i].z = (0.5f / d) * ((d + 2.0f) * external_force.z + (d - 2.0f) * oriented_force.z);

    rotational_velocities[i].x = 0.0f;
    rotational_velocities[i].y = 0.0f;
    rotational_velocities[i].z = 0.0f;
}
#endif //FIBERS_RESET_VELOCITIES_KERNEL_
