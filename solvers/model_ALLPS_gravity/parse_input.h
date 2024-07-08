struct system_parameters
{
    float dx, dy, dt, aQ, bQ, KQ, gammaQ, lambda, alpha, aPhi, bPhi, kappaPhi, MPhi, kappaHatPhi, eta, GammaFric, gGrav;
    float phi0_noise, Qxx0_noise, Qxy0_noise, phi_D, Qxx_D, Qxy_D;
    int seed, Nx, Ny, NSteps, NSave, initialConfig, noiseConfig;
};

system_parameters parser(char* paramsFile);