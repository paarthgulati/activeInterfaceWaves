from scipy.optimize import OptimizeWarning
import warnings
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import animation
from matplotlib.animation import PillowWriter
from tqdm import tqdm 
import pandas as pd 


def main():
    Nx = 1024
    Ny = 1024
    NSteps = 20_000_000
    NSave = 1000
    dt = 0.001
    alp =-1.0

    for KQ in [1.50]:
        inputDir ='../projects/activeInterfaceWaves/Nx_1024_Ny_1024_dx_1_dt_0.001_KQ_'+"{:0.2f}".format(KQ)+'/alpha_'+"{:0.2f}".format(alp)+'_seed_1'

        if os.path.exists(inputDir):
            parameters = np.loadtxt(inputDir+'/sim_parameters.txt', dtype='str')
            if parameters[12][0] == 'KQ' and parameters[15][0] == 'alpha':
                if np.isclose(float(parameters[12][2]), KQ) and np.isclose(float(parameters[15][2]), alp):
                    # NSteps = int(parameters[3][2])
                    NSave = int(parameters[4][2])
                    print('inputDir: ' + inputDir)

                    outDir = 'heightFieldData'+ inputDir.replace('../projects/activeInterfaceWaves', '', 1)
                    if os.path.exists(outDir):
                        print('outputDir: ' + outDir)

                        hField = np.load(os.path.join(outDir, 'hField.npy'))
                        
                        fig, ax = plt.subplots()
                        tx = ax.text(0,30,'t={:.1f}'.format(0.0), bbox=dict(boxstyle="round",ec='white',fc='white'))
                        P1=ax.imshow(load_data_file(inputDir, 'phi',0, Nx, Ny), cmap='RdBu_r', origin='lower')
                        cb= fig.colorbar(P1)
                        line, = ax.plot(np.arange(Nx), hField[0,:], 'k', linewidth=1)

                        def animate(t):
                            line.set_ydata(hField[t,:])
                            P1.set_data(load_data_file(inputDir, 'phi',int(t*NSave), Nx, Ny))
                            tx.set_text('t={:.1f}'.format(t*NSave*dt))
                            return fig,

                        ani = animation.FuncAnimation(fig, animate, frames= tqdm(range(NSteps//NSave-1)), interval = 1)
                        writervideo = animation.FFMpegWriter(fps=60)
                        ani.save(os.path.join(outDir, 'hField_phi.avi'), writer=writervideo, dpi=500)


def load_data_file(inputDir, s, N, Nx, Ny):
    filename = inputDir + '/data/'+s+'.csv.' + str(N)
    x = pd.read_csv(filename, delimiter=',').values
    X = x[:, 2].reshape(Ny, Nx)

    return X


if __name__ == "__main__":
    main()
