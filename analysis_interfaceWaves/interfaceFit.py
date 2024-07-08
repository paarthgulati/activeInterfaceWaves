from scipy.optimize import curve_fit
from multiprocessing import Pool
from scipy.optimize import OptimizeWarning
import warnings
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from tqdm import tqdm 
import pandas as pd 
from scipy import ndimage


def main():
    Nx = 1024
    Ny = 1024
    NSteps = 20_000_000
    NSave = 1000
    dt = 0.001
    recompute_OP = True
    aPhi = -5.0
    gGrav= 0.10
    alp =-1.0

    for KQ in [1.50]:
        inputDir ='../projects/activeInterfaceWaves/Nx_1024_Ny_1024_dx_1_dt_0.001_KQ_'+"{:0.2f}".format(KQ)+'/alpha_'+"{:0.2f}".format(alp)+'_seed_1'

        parameters = np.loadtxt(inputDir+'/sim_parameters.txt', dtype='str')
        if parameters[12][0] == 'KQ' and parameters[15][0] == 'alpha':
            if np.isclose(float(parameters[12][2]), KQ) and np.isclose(float(parameters[15][2]), alp):
                NSteps = int(parameters[3][2])
                NSave = int(parameters[4][2])
                print('inputDir: ' + inputDir)

                xThresh=50
                xRange= np.arange(0,xThresh, 1)
                num_processes = 4  #4 seemed optimal for the system size during benchmakring on the workstation

                nArray=np.array(np.arange(0, int(NSteps), int(NSave)))
                hField = np.zeros([np.size(nArray), int(Nx)])

                outDir = 'heightFieldData'+ inputDir.replace('../projects/activeInterfaceWaves', '', 1)

                if not(os.path.exists(outDir)) or recompute_OP:

                    if os.path.isdir(outDir) == 1:
                        shutil.rmtree(outDir)
                    os.makedirs(outDir)

                    print('outputDir: ' + outDir)

                    for (n, N) in enumerate(tqdm(nArray)):
                        c_n = load_data_file(inputDir, 'phi', N, Nx, Ny)

                        # #Height Field

                        sx = ndimage.sobel(c_n,axis=0,mode='constant')
                        heightField = np.argmin(sx[1:-1,:], axis=0)


                        segment_size = Nx // num_processes

                        with Pool(processes=num_processes) as pool:
                            args_list = [
                                (i, i + segment_size, c_n, heightField, xThresh, xRange) for i in range(0, Nx, segment_size)
                            ]
                            results = pool.map(process_Nx_segment, args_list)

                        hField[n, :] = np.concatenate([segment_result.flatten() for segment_result in results])


                    with open(os.path.join(outDir, 'hField.dat'), 'wb') as f:
                        np.save(f, hField)

                    with open(os.path.join(outDir, 'hField.npy'), 'wb') as f:
                        np.save(f, hField)

                        ############

 


def load_data_file(inputDir, s, N, Nx, Ny):
    filename = inputDir + '/data/'+s+'.csv.' + str(N)
    x = pd.read_csv(filename, delimiter=',').values
    X = x[:, 2].reshape(Ny, Nx)

    return X


def func(x, w, c,d):
    return c*np.tanh((x+d)/w) 

def process_Nx_segment(args):
    seg_start, seg_end, c_n, heightField, xThresh, xRange = args    
    result_for_segment = np.zeros([seg_end-seg_start, 1])
    c_n_array = np.array(c_n)  # Convert once outside the loop

    for j in range(seg_start, seg_end):
        # ... code that processes j ...
        tt = c_n_array[heightField[j]-xThresh//2:heightField[j]+xThresh//2, j]
        
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=OptimizeWarning)
            r = curve_fit(func, xRange, tt)

        result_for_segment[j-seg_start] = heightField[j] -xThresh//2 -r[0][2]
        
    return result_for_segment
    

if __name__ == "__main__":
    main()
