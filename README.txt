
						 University of Parma, Italy
							     &
					  Universitat Politènica de València, Spain

						      Valeria Todaro 
						         May 2022               


genES-MDA is an open-source Python software package for the solution of inverse
problems by means of the Ensemble Smoother with Multiple Data Assimilation (ES-MDA).
It can be used with any Python-based platform and does not require an installation. 

This directory contains the genES-MDA software and the Examples folder.
All files in the genES-MDA package need to be downloaded to the same folder preserving the directory hierarchy.


genES-MDA contains the following files:

	ESMDA.py	 : main module, it does not require modifications by the user.
	Mod.py 		 : subordinate module, it requires modifications by the user.
	InputSettings.py : subordinate module, it requires modifications by the user.
	Tools 		 : Python package it does not require modifications by the user.
	Model		 : folder that contains the forward model files, it requires modifications by the user.
	Obs.txt  	 : guidelines to set up the text input file that contains observations.	
	Par.txt 	 : guidelines to set up the text input file that contains parameters.
	Ens.txt		 : guidelines to set up the text input file that contains the initial ensemble of parameters (optional).
	Errors.txt	 : guidelines to set up the text input file that contains the ensemble of measurement errors (optional).
	R.txt 		 : guidelines to set up the text input file that contains the covariance matrix of measurement errors (optional).

The Examples folder includes practical application examples.

For more information, and to cite this tool, please refer to the following paper:

Todaro, V.; D'Oria, M.; Tanda, M.G.; Gómez-Hernández, J.J. genES-MDA: a generic open-source software package to solve inverse problems via the 
Ensemble Smoother with Multiple Data Assimilation. Computers & Geosciences, 2022, 167, 105210. https://doi.org/10.1016/j.cageo.2022.105210.
