import numpy as np

new_ens='n'            #Do you want to generate a new ensemble? (y/n)
if new_ens=='y':
    ens=500            #Ensemble size
new_err='n'            #Do you want to generate new random errors? (y/n)
N_iter=6               #Number of iterations
alpha_geo=4            #Alpha_geo
w=0                    #Relaxation coefficient
inflation='y'          #Do you want to apply the inflation? (y/n)
if inflation=='y':
    rr=1.01            #Inflation coefficient
localize='y'           #Do you want to apply the localization? (y/n)
if localize=='y':
    loc_space='y'      #Localization in space? (y/n)
    loc_time='n'       #Localization in time? (y/n)
    iter_loc='n'       #Do iterative localization? (y/n)
space_transform='y'    #Do you want to work in transformed space? (y/n)
ref_solution='y'       #Do you have the reference solution? (y/n)
last_iter_pred=0       #Get predictions for the final parameter estimation: 0(no), 1(ensemble mean), 2(all realizations)

# Function to generate the initial ensemble of parameters
def Func_ens(par, ens):
    from Tools import EnsembleGenerator
    N_par=par.shape[0]
    Ensemble=np.zeros((N_par,ens))   
    (Min,Max)=(3,12)
    Ensemble=EnsembleGenerator.ConstantRandom(Min,Max,N_par,ens) 
    return Ensemble

# Function to generate observation error realizations    
def Func_err(Obs, ens):
    from Tools import ErrorGenerator
    N_obs=Obs.shape[0]
    var_y=3e-4       # Assumed observation error variance
    # Generate Gaussian observation errors and their covariance matrix
    eps,R=ErrorGenerator.NormalError(var_y, N_obs, ens)
    return (eps,R)

# Function for forward transformation of parameters    
def forward_transf(xx):
    from Tools import Transformation as T
    xx=T.Log_forward(xx)    # Apply log-transform
    return xx

# Function for backward transformation
def backward_transf(xx):
    from Tools import Transformation as T
    xx=T.Log_backward(xx)   # Apply inverse log-transform
    return xx

# Function to compute localization matrices
def localization(ens_par,par,obs,loc_space,loc_time,iter_loc):
    a_space=50       #correlation length in space
    a_time=3*100     #correlation length in time
    time_par=par[:,3]
    pos_obs=obs[:,0:3]
    time_obs=obs[:,3]
    
    # Determine parameter locations (spatial)
    if iter_loc=='n':
        pos_par=par[:,0:3]
    else:
        ens_m=ens_par.mean(1)
        pos_par=np.tile(ens_m[0:3],(time_par.shape[0],1))
        
    from Tools import Localization as Loc
    # Compute spatial correlation matrices
    [rho_yy_sp,rho_xy_sp,rho_xx_sp]=Loc.SpaceLocal(a_space,pos_par,pos_obs)
    # Compute temporal correlation matrices
    [rho_yy_tm,rho_xy_tm,rho_xx_tm]=Loc.TimeLocal(a_time,time_par,time_obs)
    
    # Combine localization based on flags for space and time
    if loc_space=='y' and loc_time=='n':
        rho_yy=rho_yy_sp
        rho_xy=rho_xy_sp
        rho_xx=rho_xx_sp
    elif loc_space=='n' and loc_time=='y':
        rho_yy=rho_yy_tm
        rho_xy=rho_xy_tm
        rho_xx=rho_xx_tm        
    elif loc_space=='y' and loc_time=='y':
        # Combine by element-wise multiplication
        rho_yy=rho_yy_sp*rho_yy_tm
        rho_xy=rho_xy_sp*rho_xy_tm
        rho_xx=rho_xx_sp*rho_xx_tm 
        
    rho_yy[np.isnan(rho_yy)]=1
    rho_xy[np.isnan(rho_xy)]=1
    rho_xx[np.isnan(rho_xx)]=1
    return (rho_yy,rho_xy,rho_xx)

# Function to compute performance metrics using only observations
def Metrics_obs(Xprev,pred,obs):
    from Tools import Metrics as m
    par=Xprev
    #Root Mean Squared Error between predictions and observations
    RMSE_obs=[m.RMSE(obs.flatten(), pred.mean(1))]  
    #Average Ensemble Spread
    AES=[m.AES(par)] 
    
    metrics_dict ={}
    for variable in ['RMSE_obs','AES','rank_Qxx']:
        metrics_dict[variable]=eval(variable)
    return metrics_dict

# Function to compute performance metrics using both observations and reference parameters
def Metrics_obs_par(Xprev,pred,True_par,obs):
    from Tools import Metrics as m
    #Root Mean Squared Error between predictions and observations
    RMSE_obs=[m.RMSE(obs.flatten(), pred.mean(1))] 
    #Root Mean Squared Error between true and estimated parameter (ensemble mean)
    RMSE_par=[m.RMSE(True_par, Xprev.mean(1))] 
    # Average Ensemble Spread
    AES=[m.AES(Xprev)]  
    #Rank of covariance matrices 
    ens=Xprev.shape[1]
    xm=np.atleast_2d(Xprev.mean(1)).T
    ym=np.atleast_2d(pred.mean(1)).T
    Qx=Xprev-xm*np.ones((1,ens))
    Qy=pred-ym*np.ones((1,ens))
    Qxy=Qx@Qy.T/(ens-1)
    Qyy=Qy@Qy.T/(ens-1)
    Qxx=Qx@Qx.T/(ens-1)
    rank_Qxy=[np.linalg.matrix_rank(Qxy)] 
    rank_Qyy=[np.linalg.matrix_rank(Qyy)] 
    rank_Qxx=[np.linalg.matrix_rank(Qxx)] 
    
    metrics_dict ={}
    for variable in ['RMSE_obs', 'RMSE_par','AES','rank_Qxy','rank_Qyy','rank_Qxx']:
        metrics_dict[variable]=eval(variable)
    return metrics_dict