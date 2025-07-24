import numpy as np
import os

# Function to write the parameter vector to the specific input file
def write_input(par):
    with open('par.txt','wb') as file:
        np.savetxt(file, par, fmt='%.4e', newline='\n') 

# Function to run the forward model
def run():
    par=np.loadtxt('par.txt', dtype=float)
    T=par[0] # Transmissivity in m²/s
    S=par[1] # Storativity (dimensionless)
    
    from Forward_model import well_function
    well_function(T, S)     

# Function to read model output and extract predictions at observation times
def read_output():
    os.chdir('..')
    Obs_file=np.loadtxt('Obs.txt')
    os.chdir('Model')
    t_obs=Obs_file[:,3]
    t,s=np.loadtxt('pred.txt', dtype=float, unpack=True)
    idx=np.where(np.isin(t, t_obs))[0]
    pred=s[idx]
    return pred