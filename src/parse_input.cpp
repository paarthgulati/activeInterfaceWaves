
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "../inc/parse_input.h"

system_parameters parser(char* paramsFile){ 
  
  FILE* fparams = fopen(paramsFile, "r");
  if (fparams == NULL) {
    printf("ERROR: cannot open the parameter file!\n");
    exit (EXIT_FAILURE);
  }

  char line [1000];

  system_parameters sys_par;

  int nparams = 0;

  while (fgets(line, sizeof(line), fparams) != NULL) {
    nparams += sscanf(line, "seed = %d", &sys_par.seed);    
    nparams += sscanf(line, "Nx = %d", &sys_par.Nx);    
    nparams += sscanf(line, "Ny = %d", &sys_par.Ny);    
    nparams += sscanf(line, "NSteps = %d", &sys_par.NSteps);    
    nparams += sscanf(line, "NSave = %d", &sys_par.NSave);    
    nparams += sscanf(line, "initialConfig = %d", &sys_par.initialConfig);    
    nparams += sscanf(line, "boundaryConfig = %d", &sys_par.boundaryConfig);    
    nparams += sscanf(line, "noiseConfig = %d", &sys_par.noiseConfig);    
    nparams += sscanf(line, "dx = %f", &sys_par.dx);
    nparams += sscanf(line, "dy = %f", &sys_par.dy);
    nparams += sscanf(line, "dt = %f", &sys_par.dt);
    nparams += sscanf(line, "aQ = %f", &sys_par.aQ);
    nparams += sscanf(line, "bQ = %f", &sys_par.bQ); 
    nparams += sscanf(line, "KQ = %f", &sys_par.KQ);
    nparams += sscanf(line, "gammaQ = %f", &sys_par.gammaQ);
    nparams += sscanf(line, "lambda = %f", &sys_par.lambda);
    nparams += sscanf(line, "alpha = %f", &sys_par.alpha);   
    nparams += sscanf(line, "aPhi = %f", &sys_par.aPhi);   
    nparams += sscanf(line, "bPhi = %f", &sys_par.bPhi);   
    nparams += sscanf(line, "kappaPhi = %f", &sys_par.kappaPhi);   
    nparams += sscanf(line, "kappaHatPhi = %f", &sys_par.kappaHatPhi);   
    nparams += sscanf(line, "MPhi = %f", &sys_par.MPhi);   
    nparams += sscanf(line, "phiAvg = %f", &sys_par.phiAvg);   
    nparams += sscanf(line, "eta = %f", &sys_par.eta);   
    nparams += sscanf(line, "gGrav = %f", &sys_par.gGrav);   
    nparams += sscanf(line, "GammaFric = %f", &sys_par.GammaFric);   
    nparams += sscanf(line, "phi0_noise = %f", &sys_par.phi0_noise);   
    nparams += sscanf(line, "Qxx0_noise = %f", &sys_par.Qxx0_noise);   
    nparams += sscanf(line, "Qxy0_noise = %f", &sys_par.Qxy0_noise);   
    nparams += sscanf(line, "phi_D = %f", &sys_par.phi_D);   
    nparams += sscanf(line, "Qxx_D = %f", &sys_par.Qxx_D);
    nparams += sscanf(line, "Qxy_D = %f", &sys_par.Qxy_D);   
  }

  fclose(fparams);

  if (nparams != 32) {
    printf("ERROR: Incorrect number of parameters! \n");
    exit (EXIT_FAILURE);
  }  
  printf("Read parameters:\n");
  printf("seed = %d; ", sys_par.seed);
  printf("Nx = %d; ", sys_par.Nx);
  printf("Ny = %d; ", sys_par.Ny);
  printf("NSteps = %d; ", sys_par.NSteps);
  printf("NSave = %d; \n", sys_par.NSave);
  printf("intitialConfig = %d; ", sys_par.initialConfig);
  printf("boundaryConfig = %d; ", sys_par.boundaryConfig);
  printf("noiseConfig = %d; \n", sys_par.noiseConfig);
  printf("dy = %.1f; ", sys_par.dy);
  printf("dx = %.1f; ", sys_par.dx);
  printf("dt = %.2f; ", sys_par.dt);
  printf("aQ = %.2f; ", sys_par.aQ);
  printf("bQ = %.2f; ", sys_par.bQ);
  printf("KQ = %.2f; ", sys_par.KQ);
  printf("gammaQ = %.2f; ", sys_par.gammaQ);
  printf("lambda = %.1f; ", sys_par.lambda);
  printf("alpha = %.2f; ", sys_par.alpha);
  printf("aPhi = %.2f; ", sys_par.aPhi);
  printf("bPhi = %.2f; ", sys_par.bPhi);
  printf("kappaPhi = %.2f; ", sys_par.kappaPhi);
  printf("kappaHatPhi = %.2f; ", sys_par.kappaHatPhi);
  printf("MPhi = %.2f; ", sys_par.MPhi);
  printf("phiAvg = %.2f; ", sys_par.phiAvg);
  printf("eta = %.2f; ", sys_par.eta);
  printf("gGrav = %.4f; ", sys_par.gGrav);
  printf("GammaFric = %.4f; \n", sys_par.GammaFric);
  printf("phi0_noise = %.2f;", sys_par.phi0_noise);
  printf("Qxx0_noise = %.2f;", sys_par.Qxx0_noise);
  printf("Qxy0_noise = %.2f;", sys_par.Qxy0_noise);
  printf("phi_D = %.2f;", sys_par.phi_D);
  printf("Qxx_D = %.2f;", sys_par.Qxx_D);
  printf("Qxy_D = %.2f \n", sys_par.Qxy_D);

  return sys_par; 
}