/*
 *  fibers - simulates slender fibers in a fluid.
 *  Copyright (C) 2014  Eric Wolter <eric.wolter@gmx.de>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <iostream>

#include "fiberopt.h"
#include "kernels/constants.cu"
#include "parameters.h"
#include "simulation.h"

int main(int argc, char *argv[])
{
    FiberArgs args = fiberopt(argc, argv,/* help */  1, /* version */ "v0.3.0");

    int nDevices;

    cudaGetDeviceCount(&nDevices);
    std::cout << "**************************************************" << std::endl;
    std::cout << "Devices:" << std::endl;
    for (int i = 0; i < nDevices; ++i)
    {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("  Device Number: %d\n", i);
        printf("    Device name                  : %s\n", prop.name);
        printf("    Memory Clock Rate (KHz)      : %d\n",
               prop.memoryClockRate);
        printf("    Memory Bus Width (bits)      : %d\n",
               prop.memoryBusWidth);
        printf("    Peak Memory Bandwidth (GB/s) : %f\n",
               2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);

        size_t free_byte;
        size_t total_byte;
        cudaMemGetInfo( &free_byte, &total_byte );
        printf("    Total Memory (MB)            : %.0f\n",
               total_byte/1024.0/1024.0);
        printf("    Free Memory (MB)             : %.0f\n",
               free_byte/1024.0/1024.0);
        printf("    Used Memory (MB)             : %.0f\n",
               (total_byte-free_byte)/1024.0/1024.0);
    }
    std::cout << "**************************************************" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  Number of fibers                   : " << NUMBER_OF_FIBERS << std::endl;
    std::cout << "  Number of timesteps                : " << NUMBER_OF_TIMESTEPS << std::endl;
    std::cout << "  Size of timesteps                  : " << TIMESTEP << std::endl;
    std::cout << "  Slenderness                        : " << SLENDERNESS << std::endl;
    std::cout << "  Number of terms in force expansion : " << NUMBER_OF_TERMS_IN_FORCE_EXPANSION << std::endl;
    std::cout << "  Number of quadrature intervals     : " << NUMBER_OF_QUADRATURE_INTERVALS << std::endl;
    std::cout << "  State save interval                : " << STATE_SAVE_INTERVAL << std::endl;
    std::cout << "  Velocity save interval             : " << VELOCITY_SAVE_INTERVAL << std::endl;
#ifdef VALIDATE
    std::cout << "  Validating enabled                 : Yes" << std::endl;
#else
    std::cout << "  Validating enabled                 : No" << std::endl;
#endif //VALIDATE
#ifdef BENCHMARK
    std::cout << "  Benchmarking enabled               : Yes" << std::endl;
#else
    std::cout << "  Benchmarking enabled               : No" << std::endl;
#endif //BENCHMARK
    std::cout << "**************************************************" << std::endl;

    Configuration configuration = Parameters::parseConfigurationFiles(args.layout);

    Simulation simulation(configuration);

    bool running = true;
    unsigned long current_timestep = 0;
    do
    {
        std::cout << "     [CPU]      : Timestep " << current_timestep + 1 << " of " << NUMBER_OF_TIMESTEPS << std::endl;
        simulation.step(current_timestep);

        current_timestep++;

        if (current_timestep >= NUMBER_OF_TIMESTEPS)
        {
            running = false;
        }
    }
    while (running);

    // cleanup
    delete[] configuration.initial_positions;
    delete[] configuration.initial_orientations;
}
