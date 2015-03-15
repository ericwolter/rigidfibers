/*
 *  simulation.cc - contains all logic required for simulating the fibers
 *
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
#include "simulation.h"

#include <stdio.h>

#include <cmath>
#include <ctime>
#include <vector>
#include <sstream>
#include <iomanip>  // used for standard output manipulation (e.g setprecision)

#include "magma.h"

#include "resources.h"

#include "kernels/assemble_system_1D.cu"
#include "kernels/assemble_system_2D.cu"
#include "kernels/assemble_system_3D.cu"
#include "kernels/update_velocities_1D.cu"
#include "kernels/update_velocities_2D.cu"
#include "kernels/reset_velocities.cu"
#include "kernels/update_fibers_firststep.cu"
#include "kernels/update_fibers.cu"
#include "kernels/reset_system.cu"

Simulation::Simulation(Configuration configuration)
{
    configuration_ = configuration;
    performance_ = new Performance();

    global_work_size_ = IntCeil(NUMBER_OF_FIBERS, 256);

    initializeGPUMemory();

    magma_init();

    writeFiberStateToDevice();
    precomputeLegendrePolynomials();
}

Simulation::~Simulation()
{
    magma_finalize();
}

void Simulation::initializeGPUMemory()
{
#ifdef VALIDATE
    checkCuda(cudaMalloc(&gpu_validation_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6 * sizeof(int)));
    checkCuda(cudaMemset(gpu_validation_, -1, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6 * sizeof(int)));
#endif //VALIDATE

    checkCuda(cudaMalloc(&gpu_external_force_, NUMBER_OF_FIBERS * sizeof(float4)));
    float4 *host_external_force = new float4[NUMBER_OF_FIBERS];
    for (size_t m = 0; m < NUMBER_OF_FIBERS; ++m)
    {
      if (m < NUMBER_OF_FIBERS / 2) {
        host_external_force[m].x = 0;
        host_external_force[m].y = 0;
        host_external_force[m].z = -1;
        host_external_force[m].w = 0;
      } else {
        host_external_force[m].x = 0;
        host_external_force[m].y = 0;
        host_external_force[m].z = -0.75;
        host_external_force[m].w = 0;
      }
    }
    std::cout << "[CPU] --> [GPU] : Writing external force..." << std::endl;
    checkCuda(cudaMemcpy(gpu_external_force_, host_external_force, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyHostToDevice));
    delete[] host_external_force;

    checkCuda(cudaMalloc(&gpu_previous_positions_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_positions_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_next_positions_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_previous_orientations_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_orientations_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_next_orientations_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_previous_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_previous_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));
    checkCuda(cudaMalloc(&gpu_current_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4)));

    checkCuda(cudaMalloc(&gpu_a_matrix_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float)));
    checkCuda(cudaMalloc(&gpu_b_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float)));

    checkCuda(cudaMemset(gpu_a_matrix_, 0, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float)));
    checkCuda(cudaMemset(gpu_b_vector_, 0, TOTAL_NUMBER_OF_ROWS * sizeof(float)));
}

void Simulation::writeFiberStateToDevice()
{
    std::cout << "[CPU] --> [GPU] : Writing initial fiber positions..." << std::endl;
    checkCuda(cudaMemcpy(gpu_current_positions_, configuration_.initial_positions, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyHostToDevice));
    std::cout << "[CPU] --> [GPU] : Writing initial fiber orientations..." << std::endl;
    checkCuda(cudaMemcpy(gpu_current_orientations_, configuration_.initial_orientations, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyHostToDevice));
}

void Simulation::readFiberStateFromDevice()
{

}

double Simulation::calculateLegendrePolynomial(double x, unsigned int n)
{
    switch (n)
    {
    case 0:
        return 1;
    case 1:
        return x;
    case 2:
        return (1.0 / 2.0) * (3.0 * pow(x, 2) - 1.0);
    case 3:
        return (1.0 / 2.0) * (5.0 * pow(x, 3) - 3.0 * x);
    case 4:
        return (1.0 / 8.0) * (35.0 * pow(x, 4) - 30.0 * pow(x, 2) + 3.0);
    case 5:
        return (1.0 / 8.0) * (63.0 * pow(x, 5) - 70.0 * pow(x, 3) + 15.0 * x);
    case 6:
        return (1.0 / 16.0) * (231.0 * pow(x, 6) - 315.0 * pow(x, 4) + 105.0 * pow(x, 2) - 5.0);
    case 7:
        return (1.0 / 16.0) * (429.0 * pow(x, 7) - 693.0 * pow(x, 5) + 315.0 * pow(x, 3) - 35.0 * x);
    case 8:
        return (1.0 / 128.0) * (6435.0 * pow(x, 8) - 12012.0 * pow(x, 6) + 6930.0 * pow(x, 4) - 1260.0 * pow(x, 2) + 35.0);
    default:
        std::cerr << "Could not precompute legendre polynomials - n not in range [1..8]: " << n << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Simulation::precomputeLegendrePolynomials()
{
    // These are the precalculated points for a 3rd order gaussian quadrature
    // These can be looked up in the literature
    double p0 = -sqrt(15.0) / 5.0;
    double p1 = 0.0;
    double p2 = sqrt(15.0) / 5.0;

    // These are the correcponding weights also found in the literature
    double w0 = 5.0 / 9.0;
    double w1 = 8.0 / 9.0;
    double w2 = 5.0 / 9.0;

    // Intialize lower bound of the current integral to -1. At the start of the
    // subinterval iteration this is the lowest bound of the overall integral
    double lower_bound = -1.0;

    // Calculate the size of a single subinterval. The overall integral bounds
    // are [-1, 1] so the range is 2, which can simply be divided by the number
    // of subintervals.
    double interval_size = 2.0 / NUMBER_OF_QUADRATURE_INTERVALS;

    float *host_quadrature_points = new float[TOTAL_NUMBER_OF_QUADRATURE_POINTS];
    float *host_quadrature_weights = new float[TOTAL_NUMBER_OF_QUADRATURE_POINTS];
    //  On wikipedia the mapping from [a, b] to [-1, 1] is done with a factor of
    // (b - a) / 2. However in our case b = a + iv, so the factor would simply
    // be iv / 2.
    // Additionally the point as to be shifted by (a + b) / 2, which for us is
    // (a + a + iv) / 2 = (2 * a * iv) / 2.
    // So if we pull out dividing by 2 we arrive at formula below for the point
    // The weight on wikipedia is also scaled by (b - a) / 2, this being iv / 2
    // for us. If we now plug in iv = 2 / NoQI the factor simply becomes
    // 1 / NoQI. So the weights can simply be divided by the number of
    // subintervals as in the formula below
    for (size_t interval_index = 0; interval_index < NUMBER_OF_QUADRATURE_INTERVALS; ++interval_index)
    {
        // @TODO potential micro optimizations as p*, w*, interval_size
        //      number_of_quadrature_intervals are constant they could be
        //      calculated outside the loop, however for clarity we leave and
        //      here right now and precomputing polynomials is not performance
        //      critcal anyway
        // @TODO potential memory savings because weights are the same for each
        //      interval
        size_t interval_start_index = interval_index * NUMBER_OF_QUADRATURE_POINTS_PER_INTERVAL;
        host_quadrature_points[interval_start_index + 0] = (2.0 * lower_bound + interval_size + p0 * interval_size) / 2.0;
        host_quadrature_points[interval_start_index + 1] = (2.0 * lower_bound + interval_size + p1 * interval_size) / 2.0;
        host_quadrature_points[interval_start_index + 2] = (2.0 * lower_bound + interval_size + p2 * interval_size) / 2.0;

        host_quadrature_weights[interval_start_index + 0] = w0 / NUMBER_OF_QUADRATURE_INTERVALS;
        host_quadrature_weights[interval_start_index + 1] = w1 / NUMBER_OF_QUADRATURE_INTERVALS;
        host_quadrature_weights[interval_start_index + 2] = w2 / NUMBER_OF_QUADRATURE_INTERVALS;

        // std::cout << quadrature_points[interval_start_index + 0] << std::endl;
        // std::cout << quadrature_points[interval_start_index + 1] << std::endl;
        // std::cout << quadrature_points[interval_start_index + 2] << std::endl;
        // std::cout << std::endl;
        // std::cout << quadrature_weights[interval_start_index + 0] << std::endl;
        // std::cout << quadrature_weights[interval_start_index + 1] << std::endl;
        // std::cout << quadrature_weights[interval_start_index + 2] << std::endl;

        // Advance to next interval by incrementing the lower bound
        lower_bound += interval_size;
    }

    std::cout << "[CPU] --> [GPU] : Writing precomputed quadrature points..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(quadrature_points, host_quadrature_points, TOTAL_NUMBER_OF_QUADRATURE_POINTS * sizeof(float)));
    std::cout << "[CPU] --> [GPU] : Writing precomputed quadrature weights..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(quadrature_weights, host_quadrature_weights, TOTAL_NUMBER_OF_QUADRATURE_POINTS * sizeof(float)));

    // The output matrix contains the legendre polynomials evaluated at each
    // quadrature point. So for each quadrature point we calculate each
    // legendre polynomial up to the number of terms for the force expansion.
    // The results is a matrix where each row represents a point and each column
    // entry represents a legendre polynomial evaluated at that point.
    // The matrix is in column major order as is the default for GLM and GLSL.
    float *host_legendre_polynomials = new float[NUMBER_OF_TERMS_IN_FORCE_EXPANSION * TOTAL_NUMBER_OF_QUADRATURE_POINTS];
    for (size_t column_index = 0; column_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++column_index)
    {
        for (size_t point_index = 0; point_index < TOTAL_NUMBER_OF_QUADRATURE_POINTS; ++point_index)
        {
            host_legendre_polynomials[point_index + column_index * TOTAL_NUMBER_OF_QUADRATURE_POINTS] = calculateLegendrePolynomial(host_quadrature_points[point_index], column_index + 1);
            // std::cout << legendre_polynomials[point_index + column_index * total_number_of_points] << std::endl;
        }
        // std::cout << std::endl;
    }

    std::cout << "[CPU] --> [GPU] : Writing precomputed legendre polynomials..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(legendre_polynomials, host_legendre_polynomials, NUMBER_OF_TERMS_IN_FORCE_EXPANSION * TOTAL_NUMBER_OF_QUADRATURE_POINTS * sizeof(float)));

    // cleanup
    delete[] host_quadrature_points;
    delete[] host_quadrature_weights;
    delete[] host_legendre_polynomials;

    double *host_double_lambda = new double[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    double *host_double_eigen = new double[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    float *host_lambda = new float[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];
    float *host_eigen = new float[NUMBER_OF_TERMS_IN_FORCE_EXPANSION];

    double c  = log(SLENDERNESS * SLENDERNESS * M_E);
    double d  = -c;
    double e  = 2.0;
    double cc = 1.0;

    host_double_lambda[0] = 2.0;
    host_double_eigen[0] = ((d - e - cc * host_double_lambda[0]) / 2.0) / (d - cc * host_double_lambda[0]);

    host_lambda[0] = host_double_lambda[0];
    host_eigen[0] = host_double_eigen[0];

    for (size_t force_index = 1; force_index < NUMBER_OF_TERMS_IN_FORCE_EXPANSION; ++force_index)
    {
        host_double_lambda[force_index] = host_double_lambda[force_index - 1] + 2.0 / (force_index + 1);
        host_double_eigen[force_index] = ((d - e - cc * host_double_lambda[force_index]) / 2.0) / (d - cc * host_double_lambda[force_index]);

        // do all calulcations in double precision but cast to the correct GPU precision
        host_lambda[force_index] = host_double_lambda[force_index];
        host_eigen[force_index] = host_double_eigen[force_index];
    }

    std::cout << "[CPU] --> [GPU] : Writing precomputed lambda values..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(lambda, host_lambda, NUMBER_OF_TERMS_IN_FORCE_EXPANSION * sizeof(float)));

    std::cout << "[CPU] --> [GPU] : Writing precomputed eigen values..." << std::endl;
    checkCuda(cudaMemcpyToSymbol(eigen, host_eigen, NUMBER_OF_TERMS_IN_FORCE_EXPANSION * sizeof(float)));

    delete[] host_double_lambda;
    delete[] host_double_eigen;

    delete[] host_lambda;
    delete[] host_eigen;
}

void Simulation::step(size_t current_timestep)
{
    (void)current_timestep;

    assembleSystem();
#ifdef VALIDATE
    dumpLinearSystem(current_timestep);
#endif //VALIDATE

    solveSystem();
#ifdef VALIDATE
    dumpSolutionSystem(current_timestep);
#endif //VALIDATE

    updateVelocities();
#ifdef VALIDATE
    dumpVelocities(current_timestep);
#endif //VALIDATE

    updateFibers(current_timestep == 0);

    DoubleSwap(float4*, gpu_previous_translational_velocities_, gpu_current_translational_velocities_);
    DoubleSwap(float4*, gpu_previous_rotational_velocities_, gpu_current_rotational_velocities_);

    TripleSwap(float4*, gpu_previous_positions_, gpu_current_positions_, gpu_next_positions_);
    TripleSwap(float4*, gpu_previous_orientations_, gpu_current_orientations_, gpu_next_orientations_);

#if !defined(BENCHMARK) && !defined(VALIDATE)
    if(STATE_SAVE_INTERVAL > 0 && current_timestep % STATE_SAVE_INTERVAL == 0) {
      saveFibers(current_timestep);
    }
    if(VELOCITY_SAVE_INTERVAL > 0 && current_timestep % VELOCITY_SAVE_INTERVAL == 0) {
      saveVelocities(current_timestep);
    }
#endif

#ifdef VALIDATE
    dumpFibers(current_timestep);
#endif //VALIDATE

#ifdef BENCHMARK
    performance_->exportMeasurements();
#endif //BENCHMARK
}

void Simulation::assembleSystem()
{
    performance_->start("assemble_system");

    checkCuda(cudaMemset(gpu_a_matrix_, 0, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float)));
    checkCuda(cudaMemset(gpu_b_vector_, 0, TOTAL_NUMBER_OF_ROWS * sizeof(float)));

#if defined(FORCE_1D)
    std::cout << "     [GPU]      : Assembling system 1D..." << std::endl;
    assemble_system_1D <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
#ifdef VALIDATE
      gpu_validation_,
#endif //VALIDATE
      gpu_current_positions_,
      gpu_current_orientations_,
      gpu_a_matrix_,
      gpu_b_vector_
    );
#elif defined(FORCE_2D)
    dim3 block_size;
    block_size.x = 8;
    block_size.y = 8;

    dim3 grid_size;
    grid_size.x = (NUMBER_OF_FIBERS + block_size.x-1) / block_size.x;
    grid_size.y = (NUMBER_OF_FIBERS + block_size.y-1) / block_size.y;

    std::cout << "     [GPU]      : Assembling system 2D..." << std::endl;
    assemble_system_2D <<< grid_size, block_size >>> (
#ifdef VALIDATE
      gpu_validation_,
#endif //VALIDATE
      gpu_current_positions_,
      gpu_current_orientations_,
      gpu_external_force_,
      gpu_a_matrix_,
      gpu_b_vector_
    );
#else
    dim3 block_size;
    block_size.x = NUMBER_OF_TERMS_IN_FORCE_EXPANSION;
    block_size.y = 4;
    block_size.z = 4;

    dim3 grid_size;
    grid_size.x = (NUMBER_OF_TERMS_IN_FORCE_EXPANSION + block_size.x-1) / block_size.x;
    grid_size.y = (NUMBER_OF_FIBERS + block_size.y-1) / block_size.y;
    grid_size.z = (NUMBER_OF_FIBERS + block_size.z-1) / block_size.z;

    std::cout << "     [GPU]      : Assembling system 3D..." << std::endl;
    assemble_system_3D <<< grid_size, block_size >>> (
#ifdef VALIDATE
      gpu_validation_,
#endif //VALIDATE
      gpu_current_positions_,
      gpu_current_orientations_,
      gpu_a_matrix_,
      gpu_b_vector_
    );
    {
    cudaError_t cudaerr = cudaDeviceSynchronize();
    if (cudaerr != cudaSuccess)
      printf("kernel launch failed with error \"%s\".\n",
      cudaGetErrorString(cudaerr));
    }
#endif
    performance_->stop("assemble_system");
    performance_->print("assemble_system");
}

void Simulation::solveSystem()
{
#ifdef MAGMA
    std::cout << "     [GPU]      : Solving system (MAGMA + Direct)..." << std::endl;

    magma_int_t *ipiv = NULL;
    magma_int_t info = 0;
    magma_imalloc_cpu( &ipiv, TOTAL_NUMBER_OF_ROWS );

    performance_->start("solve_system");

    magma_sgesv_gpu(TOTAL_NUMBER_OF_ROWS, 1, gpu_a_matrix_, TOTAL_NUMBER_OF_ROWS, ipiv, gpu_b_vector_, TOTAL_NUMBER_OF_ROWS, &info);

    std::cout << "     [GPU]      : Info    : " << info << std::endl;

    performance_->stop("solve_system");
    performance_->print("solve_system");

#else

#ifdef GMRES
    std::cout << "     [GPU]      : Solving system (ViennaCL + GMRES)..." << std::endl;
    viennacl::linalg::gmres_tag custom_solver(1e-5, 1000, 10);
#else
    std::cout << "     [GPU]      : Solving system (ViennaCL + BiCGStab)..." << std::endl;
    viennacl::linalg::bicgstab_tag custom_solver(1e-5, 1000);
#endif // GMRES

    viennacl::matrix_base<float, viennacl::column_major> a_matrix_vienna(gpu_a_matrix_, viennacl::CUDA_MEMORY,
                                    TOTAL_NUMBER_OF_ROWS, 0, 1, TOTAL_NUMBER_OF_ROWS,
                                    TOTAL_NUMBER_OF_ROWS, 0, 1, TOTAL_NUMBER_OF_ROWS);
    viennacl::vector<float> b_vector_vienna(gpu_b_vector_, viennacl::CUDA_MEMORY, TOTAL_NUMBER_OF_ROWS);

    performance_->start("solve_system");

    b_vector_vienna = viennacl::linalg::solve(a_matrix_vienna, b_vector_vienna, custom_solver);

    performance_->stop("solve_system");
    performance_->print("solve_system");

    std::cout << "     [GPU]      : No. of iters : " << custom_solver.iters() << std::endl;
    std::cout << "     [GPU]      : Est. error   : " << custom_solver.error() << std::endl;
#endif //MAGMA

}

void Simulation::updateVelocities()
{
    performance_->start("update_velocities");

#if defined(FORCE_1D)
    std::cout << "     [GPU]      : Updating velocities 1D..." << std::endl;
    update_velocities_1D <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_b_vector_,
        gpu_current_translational_velocities_,
        gpu_current_rotational_velocities_
    );
#else
    std::cout << "     [GPU]      : Resetting velocities 2D..." << std::endl;
    reset_velocities <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
        gpu_current_orientations_,
        gpu_current_translational_velocities_,
        gpu_current_rotational_velocities_
    );

    dim3 block_size;
    block_size.x = 8;
    block_size.y = 8;

    dim3 grid_size;
    grid_size.x = (NUMBER_OF_FIBERS + block_size.x-1) / block_size.x;
    grid_size.y = (NUMBER_OF_FIBERS + block_size.y-1) / block_size.y;

    std::cout << "     [GPU]      : Updating velocities 2D..." << std::endl;
    update_velocities_2D <<< grid_size, block_size >>> (
        gpu_current_positions_,
        gpu_current_orientations_,
        gpu_b_vector_,
        gpu_external_force_,
        gpu_current_translational_velocities_,
        gpu_current_rotational_velocities_
    );
#endif

    performance_->stop("update_velocities");
    performance_->print("update_velocities");
}

void Simulation::updateFibers(bool first_timestep)
{
    std::cout << "     [GPU]      : Updating fibers..." << std::endl;

    // A second order multi-step method
    // @TODO Why? Which one?
    // The first time step is a simple forward euler

    performance_->start("update_fibers");
    if (first_timestep)
    {
        update_fibers_firststep <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
            gpu_current_positions_,
            gpu_next_positions_,
            gpu_current_orientations_,
            gpu_next_orientations_,
            gpu_current_translational_velocities_,
            gpu_current_rotational_velocities_
        );
    }
    else
    {
        update_fibers <<< (NUMBER_OF_FIBERS + 31) / 32, 32 >>> (
            gpu_previous_positions_,
            gpu_current_positions_,
            gpu_next_positions_,
            gpu_previous_orientations_,
            gpu_current_orientations_,
            gpu_next_orientations_,
            gpu_previous_translational_velocities_,
            gpu_current_translational_velocities_,
            gpu_previous_rotational_velocities_,
            gpu_current_rotational_velocities_
        );
    }
    performance_->stop("update_fibers");
    performance_->print("update_fibers");
}

void Simulation::saveFibers(size_t current_timestep)
{
  float4 *p = new float4[NUMBER_OF_FIBERS];
  float4 *o = new float4[NUMBER_OF_FIBERS];

  checkCuda(cudaMemcpy(p, gpu_current_positions_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));
  checkCuda(cudaMemcpy(o, gpu_current_orientations_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));

  std::string executablePath = Resources::getExecutablePath();

  std::stringstream output_path;
  output_path << executablePath << "/XcT_res_" << std::setfill('0') << std::setw(5) << (current_timestep+1) << ".out";

  std::ofstream output_file;

  output_file.open (output_path.str().c_str());

  output_file << NUMBER_OF_FIBERS << std::endl;
  output_file << std::fixed << std::setprecision(8);

  for (size_t row_index = 0; row_index < NUMBER_OF_FIBERS; ++row_index)
  {
    float4 p_value = p[row_index];
    float4 o_value = o[row_index];

    output_file << p_value.x << " ";
    output_file << p_value.y << " ";
    output_file << p_value.z << std::endl;
    output_file << o_value.x << " ";
    output_file << o_value.y << " ";
    output_file << o_value.z << std::endl;
  }
  output_file.close();

  delete[] p;
  delete[] o;
}

void Simulation::saveVelocities(size_t current_timestep)
{
  float4 *t = new float4[NUMBER_OF_FIBERS];
  float4 *r = new float4[NUMBER_OF_FIBERS];

  checkCuda(cudaMemcpy(t, gpu_current_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));
  checkCuda(cudaMemcpy(r, gpu_current_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));

  std::string executablePath = Resources::getExecutablePath();

  std::stringstream output_path;
  output_path << executablePath << "/TcR_res_" << std::setfill('0') << std::setw(5) << (current_timestep+1) << ".out";

  std::ofstream output_file;

  output_file.open (output_path.str().c_str());

  output_file << NUMBER_OF_FIBERS << std::endl;
  output_file << std::fixed << std::setprecision(8);

  for (size_t row_index = 0; row_index < NUMBER_OF_FIBERS; ++row_index)
  {
    float4 t_value = t[row_index];
    float4 r_value = r[row_index];

    output_file << t_value.x << " ";
    output_file << t_value.y << " ";
    output_file << t_value.z << std::endl;
    output_file << r_value.x << " ";
    output_file << r_value.y << " ";
    output_file << r_value.z << std::endl;
  }
  output_file.close();

  delete[] t;
  delete[] r;
}

void Simulation::dumpFibers(size_t current_timestep)
{
    float4 *p = new float4[NUMBER_OF_FIBERS];
    float4 *o = new float4[NUMBER_OF_FIBERS];

    checkCuda(cudaMemcpy(p, gpu_current_positions_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(o, gpu_current_orientations_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::stringstream p_output_path;
    p_output_path << executablePath << "/positions_" << current_timestep << ".out";

    std::stringstream o_output_path;
    o_output_path << executablePath << "/orientations_" << current_timestep << ".out";

    std::ofstream p_output_file;
    std::ofstream o_output_file;

    p_output_file.open (p_output_path.str().c_str());
    o_output_file.open (o_output_path.str().c_str());

    p_output_file << std::fixed << std::setprecision(8);
    o_output_file << std::fixed << std::setprecision(8);

    for (size_t row_index = 0; row_index < NUMBER_OF_FIBERS; ++row_index)
    {
        float4 p_value = p[row_index];
        float4 o_value = o[row_index];

        p_output_file << (p_value.x < 0 ? "     " : "      ") << p_value.x << std::endl;
        p_output_file << (p_value.y < 0 ? "     " : "      ") << p_value.y << std::endl;
        p_output_file << (p_value.z < 0 ? "     " : "      ") << p_value.z << std::endl;

        o_output_file << (o_value.x < 0 ? "     " : "      ") << o_value.x << std::endl;
        o_output_file << (o_value.y < 0 ? "     " : "      ") << o_value.y << std::endl;
        o_output_file << (o_value.z < 0 ? "     " : "      ") << o_value.z << std::endl;
    }
    p_output_file.close();
    o_output_file.close();

    delete[] p;
    delete[] o;
}

void Simulation::dumpLinearSystem(size_t current_timestep)
{
    float *a_matrix = new float[TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS];
    float *b_vector = new float[TOTAL_NUMBER_OF_ROWS];

    checkCuda(cudaMemcpy(a_matrix, gpu_a_matrix_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * sizeof(float), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(b_vector, gpu_b_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::stringstream a_matrix_output_path;
    a_matrix_output_path << executablePath << "/a_matrix_" << current_timestep << ".out";
    std::stringstream b_vector_output_path;
    b_vector_output_path << executablePath << "/b_vector_" << current_timestep << ".out";

    std::ofstream a_matrix_output_file;
    std::ofstream b_vector_output_file;

    a_matrix_output_file.open (a_matrix_output_path.str().c_str());
    b_vector_output_file.open (b_vector_output_path.str().c_str());

    a_matrix_output_file << std::fixed << std::setprecision(8);
    b_vector_output_file << std::fixed << std::setprecision(8);

    for (int row_index = 0; row_index < TOTAL_NUMBER_OF_ROWS; ++row_index)
    {
        for (int column_index = 0; column_index < TOTAL_NUMBER_OF_ROWS; ++column_index)
        {
            float value = a_matrix[row_index + column_index * TOTAL_NUMBER_OF_ROWS];
            if (value < 0)
            {
                a_matrix_output_file << "     " << value;
            }
            else
            {
                a_matrix_output_file << "      " << value;
            }
        }

        float value;
        value = b_vector[row_index];
        if (value < 0)
        {
            b_vector_output_file << "     " << value;
        }
        else
        {
            b_vector_output_file << "      " << value;
        }

        a_matrix_output_file << std::endl;
        b_vector_output_file << std::endl;
    }
    a_matrix_output_file.close();
    b_vector_output_file.close();

    delete[] a_matrix;
    delete[] b_vector;

#ifdef VALIDATE
    // only write the mapping on the first timestep as it will not change
    // throughout the simulation
    if(current_timestep == 0) {
        int *validation = new int[TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6];
        checkCuda(cudaMemcpy(validation, gpu_validation_, TOTAL_NUMBER_OF_ROWS * TOTAL_NUMBER_OF_ROWS * 6 * sizeof(int), cudaMemcpyDeviceToHost));

        std::string mapping_output_path = executablePath + "/current.map";
        std::ofstream mapping_output_file;

        mapping_output_file.open (mapping_output_path.c_str());

        for (int row_index = 0; row_index < TOTAL_NUMBER_OF_ROWS; ++row_index)
        {
            for (int column_index = 0; column_index < TOTAL_NUMBER_OF_ROWS; ++column_index)
            {
                mapping_output_file << "[V]"
                                    << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 0]
                                    << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 1]
                                    << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 2]
                                    << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 3]
                                    << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 4]
                                    << "|" << validation[row_index * 6 + column_index * TOTAL_NUMBER_OF_ROWS * 6 + 5]
                                    << "|" << row_index
                                    << "|" << column_index
                                    << std::endl;
            }
        }
        mapping_output_file.close();
        delete[] validation;
    }
#endif //VALIDATE
}

void Simulation::dumpSolutionSystem(size_t current_timestep)
{
    float *x_vector = new float[TOTAL_NUMBER_OF_ROWS];

    checkCuda(cudaMemcpy(x_vector, gpu_b_vector_, TOTAL_NUMBER_OF_ROWS * sizeof(float), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::stringstream x_vector_output_path;
    x_vector_output_path << executablePath << "/x_vector_" << current_timestep << ".out";

    std::ofstream x_vector_output_file;

    x_vector_output_file.open (x_vector_output_path.str().c_str());

    x_vector_output_file << std::fixed << std::setprecision(8);

    for (int row_index = 0; row_index < TOTAL_NUMBER_OF_ROWS; ++row_index)
    {
        float value = x_vector[row_index];
        if (value < 0)
        {
            x_vector_output_file << "     " << value;
        }
        else
        {
            x_vector_output_file << "      " << value;
        }

        x_vector_output_file << std::endl;
    }
    x_vector_output_file.close();

    delete[] x_vector;
}

void Simulation::dumpVelocities(size_t current_timestep)
{
    float4 *t_vel = new float4[NUMBER_OF_FIBERS];
    float4 *r_vel = new float4[NUMBER_OF_FIBERS];

    checkCuda(cudaMemcpy(t_vel, gpu_current_translational_velocities_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(r_vel, gpu_current_rotational_velocities_, NUMBER_OF_FIBERS * sizeof(float4), cudaMemcpyDeviceToHost));

    std::string executablePath = Resources::getExecutablePath();

    std::stringstream t_vel_output_path;
    t_vel_output_path << executablePath << "/t_vel_" << current_timestep << ".out";
    std::stringstream r_vel_output_path;
    r_vel_output_path << executablePath << "/r_vel_" << current_timestep << ".out";

    std::ofstream t_vel_output_file;
    std::ofstream r_vel_output_file;

    t_vel_output_file.open (t_vel_output_path.str().c_str());
    r_vel_output_file.open (r_vel_output_path.str().c_str());

    t_vel_output_file << std::fixed << std::setprecision(8);
    r_vel_output_file << std::fixed << std::setprecision(8);

    for (size_t row_index = 0; row_index < NUMBER_OF_FIBERS; ++row_index)
    {
        float4 t_value = t_vel[row_index];
        float4 r_value = r_vel[row_index];

        t_vel_output_file << (t_value.x < 0 ? "     " : "      ") << t_value.x << std::endl;
        t_vel_output_file << (t_value.y < 0 ? "     " : "      ") << t_value.y << std::endl;
        t_vel_output_file << (t_value.z < 0 ? "     " : "      ") << t_value.z << std::endl;

        r_vel_output_file << (r_value.x < 0 ? "     " : "      ") << r_value.x << std::endl;
        r_vel_output_file << (r_value.y < 0 ? "     " : "      ") << r_value.y << std::endl;
        r_vel_output_file << (r_value.z < 0 ? "     " : "      ") << r_value.z << std::endl;
    }
    t_vel_output_file.close();
    r_vel_output_file.close();

    delete[] t_vel;
    delete[] r_vel;
}
