#ifndef FIBERS_COMMON_H_
#define FIBERS_COMMON_H_

#include <cstdlib>
#include <iostream>

// includes CUDA Runtime
#include <cuda_runtime.h>

#include <assert.h>
#define DEBUG

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << std::endl;
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
// see: https://google-styleguide.googlecode.com/svn/trunk/cppguide.xml?showone=Copy_Constructors#Copy_Constructors
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);             \
    void operator=(const TypeName&)

#define IntCeil(num, divider) ((((num) + (divider) - 1) / (divider)) * (divider))
#define DoubleSwap(t, a, b) { t tmp = a; a = b; b = tmp; }
#define TripleSwap(t, a, b, c) {t tmp = a; a = b; b = c; c = tmp; }
    
#endif // FIBERS_COMMON_H_
