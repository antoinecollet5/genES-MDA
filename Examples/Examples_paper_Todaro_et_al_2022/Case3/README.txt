Recovering the release history of two simultaneous point-source groundwater contamination events


Problem Description
-------------------
This test case demonstrates the use of genES-MDA to estimate the time-dependent release history of two point-source groundwater contamination events, based on observed concentration data.
The groundwater flow is assumed known and steady-state, precomputed using MODFLOW 2005. The contaminant transport (forward model) is modeled using MT3DMS.


Dependencies
------------
All required packages and their versions are listed in requirements.txt.

To install the required packages, run:
pip install -r requirements.txt

**Note:** No installation is required for MT3DMS as the executable file (`mt3d-usgs_1.1.0_64.exe`) is already provided within the `Model` subfolder of this example.


Execution
---------
To execute the application example, the Python command ES_MDA should be run using the default settings.


For further details, please refer to:

Todaro, V.; D'Oria, M.; Tanda, M.G.; Gómez-Hernández, J.J. genES-MDA: a generic open-source software package to solve inverse problems via the Ensemble Smoother with Multiple Data Assimilation. Computers & Geosciences, 2022, 167, 105210. https://doi.org/10.1016/j.cageo.2022.105210