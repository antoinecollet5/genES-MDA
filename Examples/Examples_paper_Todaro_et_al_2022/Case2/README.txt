Estimation of a hydraulic conductivity field


Problem Description
-------------------
This example demonstrates the use of genES-MDA to estimate a zoned hydraulic conductivity field from hydraulic head data obtained through a pumping test.
The forward model for groundwater flow is MODFLOW 2005. 


Dependencies
------------
All required packages and their versions are listed in requirements.txt.

To install the required packages, run:
pip install -r requirements.txt

**Note:** No installation is required for MODFLOW 2005 as the executable file (`mf2005dbl.exe`) is already provided within the `Model` subfolder of this example.


Execution
---------
To execute the application example, the Python command ES_MDA should be run using the default settings.


For further details, please refer to:

Todaro, V.; D'Oria, M.; Tanda, M.G.; Gómez-Hernández, J.J. genES-MDA: a generic open-source software package to solve inverse problems via the Ensemble Smoother with Multiple Data Assimilation. Computers & Geosciences, 2022, 167, 105210. https://doi.org/10.1016/j.cageo.2022.105210