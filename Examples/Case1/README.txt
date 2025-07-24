Reverse flow routing


Problem Description
-------------------
This example demonstrates the use of genES-MDA to estimate the inflow hydrograph at an ungauged location in a river reach using the reverse flow routing approach.
The forward model used for simulating the routing process is the Hydrologic Engineering Center River Analysis System (HEC-RAS, version 5.0.7)
This test uses one of the official HEC-RAS example projects, which is freely available for download from the HEC-RAS website:
https://www.hec.usace.army.mil/software/hec-ras/


Dependencies
------------
All required packages and their versions are listed in requirements.txt.

To install the required packages, run:
pip install -r requirements.txt

**Note:** HEC-RAS v.5.0.7 must be installed separately on your system. The software is available free of charge from the USACE website linked above.
	  Ensure that your operating system is set to **English language and regional settings**, as HEC-RAS requires English conventions for correct execution.


Execution
---------
To execute the application example, the Python command ES_MDA should be run using the default settings.


For further details, please refer to:

Todaro, V.; D'Oria, M.; Tanda, M.G.; Gómez-Hernández, J.J. genES-MDA: a generic open-source software package to solve inverse problems via the Ensemble Smoother with Multiple Data Assimilation. Computers & Geosciences, 2022, 167, 105210. https://doi.org/10.1016/j.cageo.2022.105210