import numpy as np
import math, os
import flopy
import flopy.modflow as mf
import flopy.mt3d as mt
import flopy.utils as fu

Obs_file=np.loadtxt('Obs.txt', dtype=float)
True_par_file=np.loadtxt('Par.txt', dtype=float)

os.chdir('Model') 

# Define MODFLOW and MT3D model names and workspace directory
modelname_flow='ReleaseHistory_flow'
modelname_mt3d='ReleaseHistory_trans'
workspace='ModelFiles'

# Load the existing MODFLOW model
mfMod=mf.Modflow.load(f'{modelname_flow}.nam', model_ws=workspace,exe_name='mf2005dbl.exe')

# Load MT3D-USGS model and link with the MODFLOW model
mtMod=mt.Mt3dms.load(f'{modelname_mt3d}.nam',model_ws=workspace, version='mt3d-usgs',
                exe_name='mt3d-usgs_1.1.0_64.exe', modflowmodel=mfMod)

# Retrieve specific packages from the model object
dis=mfMod.get_package('dis')    # DIS: Discretization package 

# Extract discretization details
delr=dis.delr.array[0]          # Cell width along rows (Dx)
delc=dis.delc.array[0]          # Cell height along columns (Dy)
top=dis.top.array[0,0]          # Top elevation of the model domain
botm=dis.botm.array[-1,-1,-1]   # Bottom elevation of the lowest layer
nlay=dis.nlay                   # Number of vertical layers
nrow=dis.nrow                   # Number of rows
ncol=dis.ncol                   # Number of columns
delv=(top-botm)/nlay            # Cell thickness in vertical direction (Dz)

os.chdir('..')

# Define function to write a new input field (SSM: Source and Sink Mixing Package)
def write_input(par):
    itype = -1       #  Constant concentration
    S_loc=[1, 28, 25]   # Location of source (layer, row, column)
    dtype = np.dtype([('k', '<i8'), ('i', '<i8'), ('j', '<i8'), ('css', '<f4'), ('itype', '<i8')])  # Data type format for stress period data
    ssm_data={}
    nper=len(par)   # Number of stress periods (equal to length of parameter vector)
    
    # Build the stress period dictionary
    for n in range(0,nper,1):
        ssm_data[n]=[S_loc[0]-1]+[S_loc[1]-1]+[S_loc[2]-1]+[round(par[n],5),itype]
        
    # Remove existing SSM package (if any) and create a new one with the updated parameters
    mtMod.remove_package("ssm")
    ssm_par = flopy.mt3d.Mt3dSsm(mtMod, stress_period_data=ssm_data,dtype=dtype)
    ssm_par.write_file()

# Define function to run the MODFLOW model
def run():
    """
    Run the MODFLOW (once) and MT3D models.
    MODFLOW is only run if the Flow Transport Link (.ftl file) is not already present.
    MT3D output is reset before each run to ensure clean results.
    """
    ftl_exists = any(f.lower().endswith('.ftl') for f in os.listdir(workspace))
    
    # Run MODFLOW only if the Flow Transport Link (FTL) file is missing
    if not ftl_exists:
        mfMod.run_model(silent=False)
        
    # Remove previous MT3D output (if any)    
    try:
        os.remove(os.path.join(workspace,'MT3D001.UCN'))
    except:
        pass
    
    # Run MT3D model (silent = True hides console output)
    mtMod.run_model(silent=True) 

# Define function to extract model output at observation points
def read_output():
    N_obs=Obs_file.shape[0]
    x_obs=Obs_file[:,0]
    y_obs=Obs_file[:,1]
    z_obs=Obs_file[:,2]
    time_obs=Obs_file[:,3]
    obs_col=np.zeros((N_obs),int)
    obs_lay=np.zeros((N_obs),int)
    obs_row=np.zeros((N_obs),int)
    for ll in range(0,N_obs):
        obs_col[ll]=math.ceil(x_obs[ll]/delr)
        obs_row[ll]=nrow-math.floor(y_obs[ll]/delc)
        obs_lay[ll]=nlay-math.floor(z_obs[ll]/delv)
        
    # Load concentration data from MT3D output file
    concobj=fu.UcnFile(os.path.join(workspace,'MT3D001.UCN'))
    times_mod = np.array(concobj.get_times())
    conc = concobj.get_alldata()
    concobj.close()
    
    C=np.zeros((N_obs))
    tt=[]
    
    # Extract concentrations at observation points and corresponding times
    for ii in range(0,N_obs):
        tt+=np.where(times_mod==time_obs[ii])
        C[ii]=conc[tt[ii],obs_lay[ii]-1,obs_row[ii]-1,obs_col[ii]-1].item()
    return C