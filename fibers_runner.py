#!/usr/bin/python

import runner.fortran_runner as FORTRAN
import runner.cuda_runner as CUDA

import argparse
parser = argparse.ArgumentParser(description='Helper program for running, validating and benchmarking the rigid fibers simulation')
parser.add_argument('--verbose',action='store_true',help='be verbose')

platform_subparsers = parser.add_subparsers(title='supported platforms')

#==================================================
#
# Fortran Parsers
#
#==================================================
parser_fortran = platform_subparsers.add_parser('fortran', help='OpenMP-enabled fortran implementation')
fortran_subparsers = parser_fortran.add_subparsers(title='supported commands')
parser_fortran_run = fortran_subparsers.add_parser('run', help='run the simulation')
parser_fortran_validate = fortran_subparsers.add_parser('validate', help='validate the simulation')
parser_fortran_benchmark = fortran_subparsers.add_parser('benchmark', help='benchmark the simulation')

#--------------------------------------------------
# Fortran Run Command
#--------------------------------------------------
parser_fortran_run.set_defaults(func=FORTRAN.run)
parser_fortran_run.set_defaults(run=True)
parser_fortran_run.set_defaults(validate=False)
parser_fortran_run.set_defaults(benchmark=False)
parser_fortran_run.add_argument('configuration_file',type=str,help='the configuration file')
parser_fortran_run.add_argument('fibers_file',type=str,help='the file containing the initial positions and orientations')
parser_fortran_run.add_argument('--threads',type=int,help='the number of threads to use (default: OpenMP-default)')
fortran_run_solvers = parser_fortran_run.add_argument_group('supported solvers')
fortran_run_solver_exclusive = fortran_run_solvers.add_mutually_exclusive_group(required=True)
fortran_run_solver_exclusive.add_argument('--direct', action='store_true')
fortran_run_solver_exclusive.add_argument('--gmres', action='store_true')
fortran_run_algorithms = parser_fortran_run.add_argument_group('supported algorithms')
fortran_run_algorithm_exclusive = fortran_run_algorithms.add_mutually_exclusive_group(required=True)
fortran_run_algorithm_exclusive.add_argument('--numerical', action='store_true')
fortran_run_algorithm_exclusive.add_argument('--analytical', action='store_true')

#--------------------------------------------------
# Fortran Validate Command
#--------------------------------------------------
parser_fortran_validate.set_defaults(func=FORTRAN.validate)
parser_fortran_validate.set_defaults(run=False)
parser_fortran_validate.set_defaults(validate=True)
parser_fortran_validate.set_defaults(benchmark=False)
parser_fortran_validate.add_argument('configuration_file',type=str,help='the configuration file')
parser_fortran_validate.add_argument('fibers_file',type=str,help='the file containing the initial positions and orientations')
parser_fortran_validate.add_argument('--threads',type=int,help='the number of threads to use (default: OpenMP-default)')
fortran_validate_solvers = parser_fortran_validate.add_argument_group('supported solvers')
fortran_validate_solver_exclusive = fortran_validate_solvers.add_mutually_exclusive_group(required=True)
fortran_validate_solver_exclusive.add_argument('--direct', action='store_true')
fortran_validate_solver_exclusive.add_argument('--gmres', action='store_true')
fortran_validate_algorithms = parser_fortran_validate.add_argument_group('supported_algorithms')
fortran_validate_algorithm_exclusive = fortran_validate_algorithms.add_mutually_exclusive_group(required=True)
fortran_validate_algorithm_exclusive.add_argument('--numerical', action='store_true')
fortran_validate_algorithm_exclusive.add_argument('--analytical', action='store_true')

#--------------------------------------------------
# Fortran Benchmark Command
#--------------------------------------------------
parser_fortran_benchmark.set_defaults(func=FORTRAN.benchmark)
parser_fortran_benchmark.set_defaults(run=False)
parser_fortran_benchmark.set_defaults(validate=False)
parser_fortran_benchmark.set_defaults(benchmark=True)
parser_fortran_benchmark.add_argument('configuration_file',type=str,help='the configuration file')
parser_fortran_benchmark.add_argument('--threads',type=int,help='the number of threads to use (default: OpenMP-default)')
parser_fortran_benchmark.add_argument('--max_rse',type=float,help='the maximum relative standard error (default: 0.2)', default=0.2)
fortran_benchmark_solvers = parser_fortran_benchmark.add_argument_group('supported solvers')
fortran_benchmark_solver_exclusive = fortran_benchmark_solvers.add_mutually_exclusive_group(required=True)
fortran_benchmark_solver_exclusive.add_argument('--direct', action='store_true')
fortran_benchmark_solver_exclusive.add_argument('--gmres', action='store_true')
fortran_benchmark_algorithms = parser_fortran_benchmark.add_argument_group('supported algorithms')
fortran_benchmark_algorithm_exclusive = fortran_benchmark_algorithms.add_mutually_exclusive_group(required=True)
fortran_benchmark_algorithm_exclusive.add_argument('--numerical', action='store_true')
fortran_benchmark_algorithm_exclusive.add_argument('--analytical', action='store_true')

###################################################
#
# CUDA Parsers
#
###################################################
parser_cuda = platform_subparsers.add_parser('cuda', help='CUDA implementation')
cuda_subparsers = parser_cuda.add_subparsers(title='supported commands')
parser_cuda_run = cuda_subparsers.add_parser('run', help='run the simulation')
parser_cuda_validate = cuda_subparsers.add_parser('validate', help='validate the simulation')
parser_cuda_benchmark = cuda_subparsers.add_parser('benchmark', help='benchmark the simulation')

#--------------------------------------------------
# CUDA Run Command
#--------------------------------------------------
parser_cuda_run.set_defaults(func=CUDA.run)
parser_cuda_run.set_defaults(run=True)
parser_cuda_run.set_defaults(validate=False)
parser_cuda_run.set_defaults(benchmark=False)
parser_cuda_run.add_argument('configuration_file',type=str,help='the configuration file')
parser_cuda_run.add_argument('fibers_file',type=str,help='the file containing the initial positions and orientations')
cuda_run_solvers = parser_cuda_run.add_argument_group('supported solvers')
cuda_run_solver_exclusive = cuda_run_solvers.add_mutually_exclusive_group(required=True)
cuda_run_solver_exclusive.add_argument('--magma', action='store_true')
cuda_run_solver_exclusive.add_argument('--gmres', action='store_true')
cuda_run_solver_exclusive.add_argument('--bicgstab', action='store_true')
cuda_run_algorithms = parser_cuda_run.add_argument_group('supported algorithms')
cuda_run_algorithm_exclusive = cuda_run_algorithms.add_mutually_exclusive_group(required=True)
cuda_run_algorithm_exclusive.add_argument('--numerical', action='store_true')
cuda_run_algorithm_exclusive.add_argument('--analytical', action='store_true')
cuda_run_dimensions = parser_cuda_run.add_argument_group('supported parallization dimensions')
cuda_run_dimension_exclusive = cuda_run_dimensions.add_mutually_exclusive_group(required=True)
cuda_run_dimension_exclusive.add_argument('--D1', action='store_true')
cuda_run_dimension_exclusive.add_argument('--D2', action='store_true')
cuda_run_dimension_exclusive.add_argument('--D3', action='store_true')

#--------------------------------------------------
# CUDA Validate Command
#--------------------------------------------------
parser_cuda_validate.set_defaults(func=CUDA.validate)
parser_cuda_validate.set_defaults(run=False)
parser_cuda_validate.set_defaults(validate=True)
parser_cuda_validate.set_defaults(benchmark=False)
parser_cuda_validate.add_argument('configuration_file',type=str,help='the configuration file')
parser_cuda_validate.add_argument('fibers_file',type=str,help='the file containing the initial positions and orientations')
cuda_validate_solvers = parser_cuda_validate.add_argument_group('supported solvers')
cuda_validate_solver_exclusive = parser_cuda_validate.add_mutually_exclusive_group(required=True)
cuda_validate_solver_exclusive.add_argument('--magma', action='store_true')
cuda_validate_solver_exclusive.add_argument('--gmres', action='store_true')
cuda_validate_solver_exclusive.add_argument('--bicgstab', action='store_true')
cuda_validate_algorithms = parser_cuda_validate.add_argument_group('supported algorithms')
cuda_validate_algorithm_exclusive = cuda_validate_algorithms.add_mutually_exclusive_group(required=True)
cuda_validate_algorithm_exclusive.add_argument('--numerical', action='store_true')
cuda_validate_algorithm_exclusive.add_argument('--analytical', action='store_true')
cuda_validate_dimensions = parser_cuda_validate.add_argument_group('supported parallization dimensions')
cuda_validate_dimension_exclusive = cuda_validate_dimensions.add_mutually_exclusive_group(required=True)
cuda_validate_dimension_exclusive.add_argument('--D1', action='store_true')
cuda_validate_dimension_exclusive.add_argument('--D2', action='store_true')
cuda_validate_dimension_exclusive.add_argument('--D3', action='store_true')

#--------------------------------------------------
# CUDA Benchmark Command
#--------------------------------------------------
parser_cuda_benchmark.set_defaults(func=CUDA.benchmark)
parser_cuda_benchmark.set_defaults(run=False)
parser_cuda_benchmark.set_defaults(validate=False)
parser_cuda_benchmark.set_defaults(benchmark=True)
parser_cuda_benchmark.add_argument('configuration_file',type=str,help='the configuration file')
parser_cuda_benchmark.add_argument('--max_rse',type=float,help='the maximum relative standard error (default: 0.2)', default=0.2)
cuda_benchmark_solvers = parser_cuda_benchmark.add_argument_group('supported solvers')
cuda_benchmark_solver_exclusive = cuda_benchmark_solvers.add_mutually_exclusive_group(required=True)
cuda_benchmark_solver_exclusive.add_argument('--magma', action='store_true')
cuda_benchmark_solver_exclusive.add_argument('--gmres', action='store_true')
cuda_benchmark_solver_exclusive.add_argument('--bicgstab', action='store_true')
cuda_benchmark_algorithms = parser_cuda_benchmark.add_argument_group('supported algorithms')
cuda_benchmark_algorithm_exclusive = cuda_benchmark_algorithms.add_mutually_exclusive_group(required=True)
cuda_benchmark_algorithm_exclusive.add_argument('--numerical', action='store_true')
cuda_benchmark_algorithm_exclusive.add_argument('--analytical', action='store_true')
cuda_benchmark_dimensions = parser_cuda_benchmark.add_argument_group('supported parallization dimensions')
cuda_benchmark_dimension_exclusive = cuda_benchmark_dimensions.add_mutually_exclusive_group(required=True)
cuda_benchmark_dimension_exclusive.add_argument('--D1', action='store_true')
cuda_benchmark_dimension_exclusive.add_argument('--D2', action='store_true')
cuda_benchmark_dimension_exclusive.add_argument('--D3', action='store_true')

print "      _      _    _    __ _ _                 "
print "  _ _(_)__ _(_)__| |  / _(_) |__  ___ _ _ ___ "
print " | '_| / _` | / _` | |  _| | '_ \/ -_) '_(_-< "
print " |_| |_\__, |_\__,_| |_| |_|_.__/\___|_| /__/ "
print "       |___/                                  "
print "                                              "

args = parser.parse_args()
args.func(args)
