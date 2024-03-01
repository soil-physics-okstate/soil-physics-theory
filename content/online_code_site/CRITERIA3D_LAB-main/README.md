# CRITERIA3D Lab
Python version of CRITERIA-3D agro-hydrological model, allows to configure 1D, 2D and 3D hydrological experiments.  
- Numerical solution for three-dimensional water flow in the soil, with management of different soil horizons.  
- Soil water retention curves: Campbell, modified Van Genucthen.  
- Sink/source: evapotranspiration, precipitation, drip irrigation.  
- Boundary conditions: surface runoff, free drainage, prescribed total potential (watertable), no flux.  
- Crop parameters: monthly leaf area index (LAI), maximum crop coefficient (kcmax), readily available water fraction (fRAW), shape factors of the root system.
- Save the output (water content / water potential) at specific points in the domain and assimilate observed water potential values.  


![](https://github.com/ARPA-SIMC/CRITERIA3D_LAB/blob/main/doc/criteria3d.png)

## Requirements
- Cython
- VPython
- numpy  
- scipy  
- pandas

## How to compile solverC.pyx
Requires Cython library and a C compiler, on Windows you can use MinGW or Visual Studio compiler.  
>cd src  
>python cythonSetup.py build_ext --inplace
 
If Cython library is missing, you can execute code set:  
> CYTHON = False

in the **commonConst.py** module.

## Authors
- Fausto Tomei    <ftomei@arpae.it>
- Marco Bittelli  <marco.bittelli@unibo.it>
- Joseph Giovanelli <josephgiovanelli@gmail.com>

## References
Bittelli, M., Campbell, G. S., & Tomei, F. (2015). Soil physics with Python: transport in the soil-plant-atmosphere system. OUP Oxford.
