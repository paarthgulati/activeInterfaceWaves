/* Short explanation on the way to add fields and terms in this file
 *
 * evolver is a class that merely calls updates on all fields and terms
 * the arguments on its constructor are 
 *
 *      evolver system(x,           sx,             sy,             dx,       dy,       dt);
 *                     Use CUDA | x-system size | y-system size | delta_x | delta_y | delta_t
 *
 * To this evolver we can add fields:
 *
 *      system.createField( name, dynamic );
 *
 * name is a string and dynamic if a boolean that sets whether the field
 * is set in each step through a time derivative or through an equality.
 *
 * To each field we can add terms
 *      
 *      system.createTerm(  field_name, prefactor, {field_1, ..., field_n}  );
 *
 *  This term would be a term of "field_name", with that prefactor, that multiplies
 *  fields field_1 to field_n.
 */ 

#include <cmath>
#include <cstdlib>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <iostream>
#include <ostream>
#include "../../inc/defines.h"
#include "../../inc/evolver.h"
#include "../../inc/field.h"
#include "../../inc/term.h"
#include "parse_input.h"

#ifdef WITHCUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif


void zero_boundaries_y(float2*, int, int);
__global__ void zero_k_y(float2*, int, int);

void phi_boundaries_y(float2*, int, int);
__global__ void phi_k_y(float2*, int, int);

int main (int argc, char* argv[]) 
{
    bool use_GPU = true; 

    if (argc < 2) 
    {
    printf("Usage: model param_file level\n");
    return 1;
    }

    system_parameters sys_par= parser(argv[1]);

    evolver system(use_GPU, sys_par.Nx, sys_par.Ny, sys_par.dx, sys_par.dy, sys_par.dt, sys_par.NSave);

    system.createField("iqxQxx", false);// 0
    system.createField("iqyQxx", false);// 1
    system.createField("iqxQxy", false);// 2
    system.createField("iqyQxy", false);// 3
    system.createField("iqxphi", false);// 4
    system.createField("iqyphi", false);// 5
    system.createField("sigxx", false); // 6
    system.createField("sigxy", false); // 7
    system.createField("vx", false);    // 8
    system.createField("vy", false);    // 9
    system.createField("wxy", false);   // 10
    system.createField("Q2", false);    // 11
    system.createField("Qxx", true);    // 12
    system.createField("Qxy", true);    // 13
    system.createField("phi", true);    // 14

    system.createField("ident", false);    // 15
    system.createField("FgY", false);    // 16
    // system.createField("P_g", false);    // 17
    // system.createField("P_c", false);    // 17

    // CONSTANTS
    // v and Q
    float gamma_fric = sys_par.GammaFric;
    float eta = sys_par.eta;
    float aQ = sys_par.aQ;
    float bQ = sys_par.bQ;
    float kQ = sys_par.KQ;
    float lambda = sys_par.lambda;
    float gamma = sys_par.gammaQ;
    float alpha = sys_par.alpha;
    // phi
    float a = sys_par.aPhi;
    float b = sys_par.bPhi;
    float M = sys_par.MPhi;
    float phi0 = std::sqrt(-a/b);
    float k = sys_par.kappaPhi;
    float ka = sys_par.kappaHatPhi;


    system.fields[8]->hasBC = true;
    system.fields[8]->boundary = zero_boundaries_y;        
    system.fields[9]->hasBC = true;
    system.fields[9]->boundary = zero_boundaries_y;        
  
    system.fields[14]->hasBC = true;
    system.fields[14]->boundary = phi_boundaries_y;        

    // Implicit terms
    system.fields[8]->implicit.push_back({eta, 1, 0, 0, 0});
    system.fields[8]->implicit.push_back({gamma_fric, 0, 0, 0, 0});
    system.fields[9]->implicit.push_back({eta, 1, 0, 0, 0});
    system.fields[9]->implicit.push_back({gamma_fric, 0, 0, 0, 0});
    system.fields[12]->implicit.push_back({-aQ/gamma});
    system.fields[12]->implicit.push_back({-kQ/gamma, 1, 0, 0, 0});
    system.fields[13]->implicit.push_back({-aQ/gamma});
    system.fields[13]->implicit.push_back({-kQ/gamma, 1, 0, 0, 0});
    system.fields[14]->implicit.push_back({-M*a, 1, 0, 0, 0});
    system.fields[14]->implicit.push_back({-M*k, 2, 0, 0, 0});

    //Explicit terms
    system.createTerm("iqxQxx", {{1.0f, 0, 1, 0, 0}}, {"Qxx"});
    system.createTerm("iqyQxx", {{1.0f, 0, 0, 1, 0}}, {"Qxx"});
    system.createTerm("iqxQxy", {{1.0f, 0, 1, 0, 0}}, {"Qxy"});
    system.createTerm("iqyQxy", {{1.0f, 0, 0, 1, 0}}, {"Qxy"});
    system.createTerm("iqxphi", {{1.0f, 0, 1, 0, 0}}, {"phi"});
    system.createTerm("iqyphi", {{1.0f, 0, 0, 1, 0}}, {"phi"});

    system.createTerm("sigxx", {{alpha/2.0f}}, {"Qxx"});
    system.createTerm("sigxy", {{alpha/2.0f}}, {"Qxy"});
    system.createTerm("sigxx", {{alpha/(2.0f*phi0)}}, {"phi", "Qxx"});
    system.createTerm("sigxy", {{alpha/(2.0f*phi0)}}, {"phi", "Qxy"});
    system.createTerm("sigxx", {{-ka/2.0f}}, {"iqxphi", "iqxphi"});
    system.createTerm("sigxx", {{ka/2.0f}}, {"iqyphi", "iqyphi"});
    system.createTerm("sigxy", {{-ka}}, {"iqxphi", "iqyphi"});

    // alc backflow
    system.createTerm("sigxx", {{lambda*aQ}},{"Qxx"});
    system.createTerm("sigxx", {{lambda*bQ}},{"Q2", "Qxx"});
    system.createTerm("sigxx", {{lambda*kQ, 1, 0, 0, 0}},{"Qxx"});
    system.createTerm("sigxy", {{lambda*aQ}},{"Qxy"});
    system.createTerm("sigxy", {{lambda*bQ}},{"Q2", "Qxy"});
    system.createTerm("sigxy", {{lambda*kQ, 1, 0, 0, 0}},{"Qxy"});

    system.createTerm("Qxx", {{lambda, 0, 1, 0, 0}}, {"vx"});
    system.createTerm("Qxx", {{-2.0f}}, {"Qxy", "wxy"});
    system.createTerm("Qxx", {{-bQ/gamma}}, {"Q2", "Qxx"});
    system.createTerm("Qxx", {{-1.0f}}, {"vx", "iqxQxx"});
    system.createTerm("Qxx", {{-1.0f}}, {"vy", "iqyQxx"});

    system.createTerm("Qxy", {{lambda/2, 0, 1, 0, 0}}, {"vy"});
    system.createTerm("Qxy", {{lambda/2, 0, 0, 1, 0}}, {"vx"});
    system.createTerm("Qxy", {{2.0f}}, {"Qxx", "wxy"});
    system.createTerm("Qxy", {{-bQ/gamma}}, {"Q2", "Qxy"});
    system.createTerm("Qxy", {{-1.0f}}, {"vx", "iqxQxy"});
    system.createTerm("Qxy", {{-1.0f}}, {"vy", "iqyQxy"});

    system.createTerm("wxy", {{0.5f, 0, 1, 0, 0}}, {"vy"});
    system.createTerm("wxy", {{-0.5f, 0, 0, 1, 0}}, {"vx"});

    system.createTerm("Q2", {{1.0f}}, {"Qxx", "Qxx"});
    system.createTerm("Q2", {{1.0f}}, {"Qxy", "Qxy"});

    system.createTerm("phi", {{-M*b, 1, 0, 0, 0}}, {"phi", "phi", "phi"});
    system.createTerm("phi", {{-1.0f}}, {"vx", "iqxphi"});
    system.createTerm("phi", {{-1.0f}}, {"vy", "iqyphi"});

    // Terms for vx and vy
    pres iqx = {1.0f, 0, 1, 0, 0};
    pres iqy = {1.0f, 0, 0, 1, 0};
    pres miqy = {-1.0f, 0, 0, 1, 0};
    pres miqy3 = {-1.0f, 0, 0, 3, 2};
    pres iqx3 = {1.0f, 0, 3, 0, 2};
    pres iqx2iqy = {1.0f, 0, 2, 1, 2};
    pres miqxiqy2 = {-1.0f, 0, 1, 2, 2};
    pres iqxiqy2 = {1.0f, 0, 1, 2, 2};
    system.createTerm("vx", {iqx, iqx3, miqxiqy2}, {"sigxx"});
    system.createTerm("vx", {iqy, iqx2iqy, iqx2iqy}, {"sigxy"});
    system.createTerm("vy", {miqy, miqy3, iqx2iqy}, {"sigxx"});
    system.createTerm("vy", {iqx, iqxiqy2, iqxiqy2}, {"sigxy"});

    // Gravitational Terms
    float g = sys_par.gGrav;
    float rho_1 = 2.0f;
    float rho_2 = 1.0f;
   
    system.createTerm("FgY", {{-(g*(rho_1-rho_2)/(2*phi0))}}, {"phi"});
    system.createTerm("FgY", {{-(g*(rho_1-rho_2)/2)}}, {"ident"});
    system.createTerm("vx", {{1.0f, 0, 1, 1, 2}}, {"FgY"});
    system.createTerm("vy", {{-1.0f, 0, 2, 0, 2}}, {"FgY"});

    switch (sys_par.initialConfig)
    {
    case 1:
        // Flat interface initial
        std::srand(sys_par.seed);
        for (int i = 0; i < sys_par.Nx; i++)
        {
            for (int j = 0; j < sys_par.Ny; j++)
            {
                system.fields[14]->real_array[j*sys_par.Nx+i].x = -std::tanh(j-sys_par.Ny/2) + sys_par.phi0_noise * 0.01f * (float)(std::rand() % 200 - 100);;
                system.fields[14]->real_array[j*sys_par.Nx+i].y = 0.0f;

                system.fields[12]->real_array[i*sys_par.Nx+j].x = sys_par.Qxx0_noise * 0.01f * (float)(std::rand() % 200 - 100);
                system.fields[12]->real_array[i*sys_par.Nx+j].y = 0.0f;
                system.fields[13]->real_array[i*sys_par.Nx+j].x = sys_par.Qxy0_noise * 0.01f * (float)(std::rand() % 200 - 100);
                system.fields[13]->real_array[i*sys_par.Nx+j].y = 0.0f;
            }
        }

    break;
    }


    for (int i = 0; i < sys_par.Nx; i++)
    {
        for (int j = 0; j < sys_par.Ny; j++)
        {
            system.fields[15]->real_array[j*sys_par.Nx+i].x =  1.0f;
            system.fields[15]->real_array[j*sys_par.Nx+i].y =  0.0f;
        }
    }


    cudaMemcpy(system.fields[12]->real_array_d, system.fields[12]->real_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[12]->comp_array_d, system.fields[12]->comp_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[13]->real_array_d, system.fields[13]->real_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[13]->comp_array_d, system.fields[13]->comp_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[14]->real_array_d, system.fields[14]->real_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[14]->comp_array_d, system.fields[14]->comp_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[15]->real_array_d, system.fields[15]->real_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    cudaMemcpy(system.fields[15]->comp_array_d, system.fields[15]->comp_array, sys_par.Nx*sys_par.Ny*sizeof(float2), cudaMemcpyHostToDevice);
    system.fields[12]->toComp();
    system.fields[13]->toComp();
    system.fields[14]->toComp();
    system.fields[15]->toComp();


    switch (sys_par.noiseConfig)
    {
    case 1:
    //Conserved Phase Field Noise
        system.fields[14]->isNoisy = true;
        system.fields[14]->noiseType = GaussianWhite;
        system.fields[14]->noise_amplitude = {sys_par.phi_D,1,0,0,0};
        break;
    
    default:
    // No noise in the dynamics
        break;
    }

    
    for (int i = 0; i < system.fields.size(); i++)
    {
        system.fields[i]->prepareDevice();
        system.fields[i]->precalculateImplicit(system.dt);
        system.fields[i]->outputToFile = false;
    }
    // system.fields[11]->outputToFile = true;
    // system.fields[12]->outputToFile = true;
    // system.fields[13]->outputToFile = true;
    system.fields[14]->outputToFile = true;

    int steps = sys_par.NSteps;
    int check = steps/100;
    if (check < 1) check = 1;
    
    if (argv[2] == "v" )
    {
        system.printInformation();
    }


    for (int i = 0; i < steps; i++)
    {
        system.advanceTime();
        if (i % check == 0)
        {
            std::cout << "Progress: " << i/check << "%\r";
            std::cout.flush();
        }
    }

    return 0;
}


void zero_boundaries_y(float2 *real_array, int sx, int sy)
{
    dim3 TPB(32,32);
    dim3 blocks((sx+31)/32, (sy+31)/32);
    zero_k_y<<<blocks, TPB>>>(real_array, sx, sy);
}

__global__ void zero_k_y(float2 *real_array, int sx, int sy)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int j = blockIdx.y * blockDim.y + threadIdx.y;
        int index = j*sx+i;

        if (index < sx*sy )
        {
            if (j < sy/8 || j > sy - sy/8)
                real_array[index].x = 0.0f;
        }
}


void phi_boundaries_y(float2 *real_array, int sx, int sy)
{
    dim3 TPB(32,32);
    dim3 blocks((sx+31)/32, (sy+31)/32);
    phi_k_y<<<blocks, TPB>>>(real_array, sx, sy);
}

__global__ void phi_k_y(float2 *real_array, int sx, int sy)
{
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int j = blockIdx.y * blockDim.y + threadIdx.y;
        int index = j*sx+i;

        if (index < sx*sy )
        {
            if (j < sy/8)
            {
                real_array[index].x = 1.0f;
            }

            if (j > sy - sy/8)
            {
                real_array[index].x = -1.0f;
            }
        }
}

