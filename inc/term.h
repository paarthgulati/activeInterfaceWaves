#ifndef TERM_H
#define TERM_H

#include <vector>
#include <fftw3.h>
#include <cuda_runtime.h>
#include <cufft.h>

#include "defines.h"

class field;


class term
{
private:
    int toReal();           // Backwards Fourier transform
    int toComp();           // Forwards Fourier transform
    int normalize();        // Normalize always after toReal()
    int computeProduct();   // Loop through product fields and multiply into real
    int copyComp();         // For when there's only one field
    int applyPrefactors();  // Fourier prefactors
    int applyPres_vector(); // Same as above but multiple
    int precomputePrefactors();
    const int sx, sy;             // system size
    const float dx, dy;           // cell size
    const float stepqx, stepqy;   // smallest wavelength, 

    // fftw stuff
    fftwf_plan plan_forward, plan_backward;

    cufftHandle plan_gpu;
public:
    term(int, int, float, float);
    int prepareDevice();

    bool isCUDA;

    pres prefactors;
    std::vector<pres> prefactors_h;

    pres *prefactors_d;


    float2 *term_real;
    float2 *term_comp;
    float2 *term_real_d;
    float2 *term_comp_d;

    std::vector<field *> product;
    float2 **product_h;
    float2 **product_d;  // Will just contain pointers to fields as an array
                        // Must be of size 'product.size() * sizeof(float2)'
    float *precomp_prefactor_h;
    float *precomp_prefactor_d;

    int multiply_by_i_pre;

    int update();
};

#endif
