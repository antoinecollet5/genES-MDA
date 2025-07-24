import numpy as np
import math, os
import flopy.modflow as mf
import flopy.utils as fu

Obs_file=np.loadtxt('Obs.txt', dtype=float)
True_par_file=np.loadtxt('Par.txt', dtype=float)

os.chdir('Model')

# Define MODFLOW model name and workspace directory
modelname='3Dexample'
workspace='ModelFiles'

# Load the existing MODFLOW model using the name file (.nam)
# 'mf2005dbl.exe' is the MODFLOW 2005 executable (double precision version)
mfMod=mf.Modflow.load(f'{modelname}.nam', model_ws=workspace,exe_name='mf2005dbl.exe')

# Retrieve specific MODFLOW packages from the model object
lpf=mfMod.get_package('lpf')    # LPF: Layer Property Flow package (defines hydraulic properties)
dis=mfMod.get_package('dis')    # DIS: Discretization package (grid structure)
wel=mfMod.get_package('wel')    # WEL: Well package (defines pumping/injection)

# Extract discretization details
delr=dis.delr.array[0]          # Cell width along rows (Dx)
delc=dis.delc.array[0]          # Cell height along columns (Dy)
top=dis.top.array[0,0]          # Top elevation of the model domain
botm=dis.botm.array[-1,-1,-1]   # Bottom elevation of the lowest layer
nlay=dis.nlay                   # Number of vertical layers
nrow=dis.nrow                   # Number of rows
ncol=dis.ncol                   # Number of columns
delv=(top-botm)/nlay            # Cell thickness in vertical direction (Dz)

# Extract coordinates and parameter values from parameter file
X=True_par_file[:,0]
Y=True_par_file[:,1]
Z=True_par_file[:,2]
var_val=True_par_file[:,4]

# Create sorted unique coordinate lists
x_vals=np.sort(np.unique(X))
y_vals=np.sort(np.unique(Y))
z_vals=np.sort(np.unique(Z))

# Find grid indices for each parameter location
z_idx = np.searchsorted(z_vals, Z)[::-1]
x_idx = np.searchsorted(x_vals, X)
y_idx = np.searchsorted(y_vals, Y)[::-1]

os.chdir('..')

# Define function to write a new input field (hydraulic conductivity)
def write_input(par):
    kfield = np.full((nlay, nrow, ncol), np.nan)    # Initialize hydraulic conductivity array
    kfield[z_idx, y_idx, x_idx,] = par              # Assign parameter values at specified indices
    mfMod.remove_package("lpf")                     # Remove existing LPF package
    
    # Create new LPF package with updated conductivity field
    lpf_par=mf.ModflowLpf(mfMod, laytyp=0, chani=-1, layvka=1, hk=kfield,vka=1, hdry=-888.0)        
    lpf_par.write_file()

# Define function to run the MODFLOW model    
def run():
    mfMod.run_model(silent=True)    # Run the model silently

# Define function to extract model output at observation points
def read_output():
    N_obs=Obs_file.shape[0]
    x_obs=Obs_file[:,0]
    y_obs=Obs_file[:,1]
    z_obs=Obs_file[:,2]
    obs_col=np.zeros((N_obs),int)
    obs_raw=np.zeros((N_obs),int)
    obs_lay=np.zeros((N_obs),int)
    for ll in range(0,N_obs):
        obs_col[ll]=math.ceil(x_obs[ll]/delr)
        obs_raw[ll]=nrow-math.floor(y_obs[ll]/delc)
        obs_lay[ll]=nlay-math.floor(z_obs[ll]/delv)
        
    # Read head data from binary output file
    headobj=fu.binaryfile.HeadFile(os.path.join(workspace,modelname+'.hds'))
    head = headobj.get_alldata()
    headobj.close()
    
    # Extract head values at observation locations
    H=np.zeros((N_obs))
    for i in range(0,N_obs):
        H[i]=head[0,obs_lay[i]-1,obs_raw[i]-1,obs_col[i]-1]
    return H