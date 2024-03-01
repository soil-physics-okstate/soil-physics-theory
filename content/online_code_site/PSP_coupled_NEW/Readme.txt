In this file are listed the updates for the program PSP_coupled.py (Chapter 12.).

1. Updated the closing of the plotting interface.
2. Updated the output printed in hours instead of seconds.
3. Included a verification on program quitting.
4. Included an option to de-activate the output plotting
   to speed up computation. Select the option in the
   file PSP_public.py by choosing isPlotActivated= True
   or isPlotActivated= False
5. Included an option to print output data to file. 
   The data are: water content, water potential and temperature
   as function of time. The user can select the depth for the
   output data.  
   
   The modification included in the file PSP_public.py are listed
   below: 

   isQuitRequired = False                 
   isPlotActivated = True              #An option to avoid plotting the graph                                      
   outputTimeStep = 1                  #[hours] Select the time step for the output
   outputFileName = "output_soil.csv"  # Select the file name for output of water content, potential and temp. 
   outputDepth = 0.1                   #[m] Select the soil depth at which the output data are saved.

6. There are two different weather records that can be used to run the model: weather.dat and SPC_weather.txt. 
   They are weather records from two experimental stations in the Emilia-Romagna region.
   The soil weather.dat is from the station listed in the paper

   2007. Pieri L., M. Bittelli, J. Q. Wu, S. Dun, D. C. Flanagan, P. Rossi Pisa, F. Ventura, and F. Salvatorelli.    
   Using the Water Erosion Prediction Project (WEPP) Model to Simulate Field-Observed Runoff and Erosion in the   
   Apennines Mountain Range, Italy. J. Hydrol., 336, 84-97.  

   while the file SPC_weather.txt is from the station of San Pietro Capofiume (SPC) described in the paper below.

   
   2011. Brocca L., S. Hasenauer, T. Lacava, F. Melone, T. Moramarco, W. Wagner, W. Dorigo, P. Matgen, J.   
   Martínez-Fernández, P. Llorens, J. Latron, C. Martin, M. Bittelli. Soil moisture estimation through ASCAT and 
   AMSR-E sensors: An intercomparison and validation study across Europe, Remote Sens. Environ. 115(12), 
   3390-3408.  

   