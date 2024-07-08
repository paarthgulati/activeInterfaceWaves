#!/bin/bash
projectName="activeInterfaceWaves"
solver_name="model_ALLPS_gravity"

rootDir=$(pwd)
log_message="Solver used: ${solver_name}"

################################################################################################
############CHANGE FROM THE DEFAULT VALUES LISTED IN solvers/${solver_name}/default_parameters.txt
seed=1 #Random seed to set up the initial noise perturbation

Nx=1024 #System size
Ny=1024

NSteps=20000000 #Total Number of time steps 
NSave=1000 #Output data files every NSave
dt=0.001 #time step size

kappaPhi=7.50 #Cahn Hilliard Parameters
aPhi=-5.0 #aPhi < 0 for phase separating fluid mixture
bPhi=5.0 #bPhi=|aPhi|
kappaHatPhi=7.50 #Capillary stress coeffcient (explicitly set equal to kappaPhi for modelH)

alpha=-1.00 #Activity of the Active Phase (alpha <0 for extensile activity) 
KQ=1.50 #Elastic Modulus of the Nematic
gGrav=0.10

while getopts ":K:h" option; do
  case $option in
    K)
      KQ=$OPTARG
      ;;
    h)
      echo "Usage: $0 [-K KQ]"
      exit 1
      ;;
  esac
done

DIRTREE=Nx_${Nx}_Ny_${Ny}_dt_${dt}_KQ_${KQ}/alpha_${alpha}_seed_${seed}  #OUTPUT DATA IS SAVED IN A DIRETREE LABELLED AS SUCH. ADAPT BASED ON THE PAREMETERS THAT ARE BEING MODIFIED

################################################################################################
################################################################################################
PROJECTDIR=projects/$projectName
echo "Using solver ${solver_name}."
echo "Working in $PROJECTDIR"

#CHECK FOR PROPER PROJECT INITIALIZATION
if [ ! -d "$PROJECTDIR" ]; then
  echo "$PROJECTDIR does not exist. Initializing the project directory."
  mkdir -p $PROJECTDIR
fi

if [ ! -d "$PROJECTDIR/$DIRTREE" ]; then
  echo "Initializing $PROJECTDIR/$DIRTREE"
  mkdir -p $PROJECTDIR/$DIRTREE
else
  echo "$PROJECTDIR/$DIRTREE already exists! Overwriting."
  rm -rf $PROJECTDIR/$DIRTREE
  mkdir $PROJECTDIR/$DIRTREE
fi

cp solvers/${solver_name}/default_parameters.txt $PROJECTDIR/$DIRTREE/sim_parameters.txt 
mkdir $PROJECTDIR/$DIRTREE/data

cd $PROJECTDIR/$DIRTREE

################################################################################################
################################################################################################
########## $PROJECTDIR/$DIRTREE/sim_parameters.txt, THE DEFAULT PARAMETER VALUES ARE REPLACED BY THE PARMETER VALUES LISTED HERE

sed -i 's/^seed .*$/seed = '${seed}'/' sim_parameters.txt
sed -i 's/^Nx .*$/Nx = '${Nx}'/' sim_parameters.txt
sed -i 's/^Ny .*$/Ny = '${Ny}'/' sim_parameters.txt
sed -i 's/^NSteps .*$/NSteps = '${NSteps}'/' sim_parameters.txt
sed -i 's/^NSave .*$/NSave = '${NSave}'/' sim_parameters.txt
sed -i 's/^KQ .*$/KQ = '${KQ}'/' sim_parameters.txt
sed -i 's/^alpha .*$/alpha = '${alpha}'/' sim_parameters.txt
sed -i 's/^dt .*$/dt = '${dt}'/' sim_parameters.txt
sed -i 's/^kappaHatPhi .*$/kappaHatPhi = '${kappaHatPhi}'/' sim_parameters.txt
sed -i 's/^kappaPhi .*$/kappaPhi = '${kappaPhi}'/' sim_parameters.txt
sed -i 's/^aPhi .*$/aPhi = '${aPhi}'/' sim_parameters.txt
sed -i 's/^bPhi .*$/bPhi = '${bPhi}'/' sim_parameters.txt
sed -i 's/^gGrav .*$/gGrav = '${gGrav}'/' sim_parameters.txt

################################################################################################
################################################################################################

echo "Simulation Initialized! Running."
printf "Solver used: ${solver_name}" >> logf.txt

################################################################################################
################################################################################################
# THE SIMULATION IS NOW RUN USING THE PARAMETERS FILE $PROJECTDIR/$DIRTREE/sim_parameters.txt

${rootDir}/bin/${solver_name} sim_parameters.txt

################################################################################################
################################################################################################
# (OPTIONAL): CREATE A MOVIE FROM THE OUTPUT PHI DATA

# echo "Simulation Finished! Creating avi files."
# if [ ! -d "analysis/movies" ]; then
#   mkdir -p analysis/movies
# fi
# python ${rootDir}/utils/movie_2dfunc.py --Nx ${Nx} --Ny ${Ny} --NSteps ${NSteps} --NSave ${NSave} --dt ${dt} -l "phi"
