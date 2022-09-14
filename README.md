# EBSD-EDS-methodologies
----- GENERAL -----
These scripts were developed to automatically analyse large electron backscatter diffraction (EBSD) 
datasets using the MTEX toolbox for Matlab. 

Explanation of the underlying methodologies are published here:
[1] *****ADD PAPER******

This work is licensed under a Creative Commons Attribution 4.0 International (CC BY 4.0). ******ASK SOPHIE IF CORRECT******
To attribute the material to the creators, please cite the above mentioned paper B. Schulz et al. (2023).

----- INSTALLATION -----
1) Install Matlab 2019b
   (newer versions should work but have not been tested)
2) Download mtex-5.7.0 (https://mtex-toolbox.github.io/download)
   (newer versions should work but have not been tested)
3) start Matlab & change current folder to MTEX folder
4) type 'startup_mtex' into commend window
5) press 'install mtex' on output

----- Important Information -----
Both these files must be generally located in the current Matlab folder:
*) ebsd_analysis.m: Matlab class that includes all methods to analyse the dataset
*) ebsd_execution.m: Executable file that calls the methods

*) The ebsd_execution.m in the current folder includes general explanations to the syntax of some basic methods (e.g. plotting)
*) The tutorials folder include instructions on how to recreate the analyses presented in B. Schulz et al. (2022)

*) Most functions do not directly output into the command window but save the output into analysis.Report

Important: The here developed code was just tested for Ni-base superalloys. For other systems changes 
in the ebsd_analysis.m class may be necessary and bugs may appear. If any issues arise, please contact us!

If there is general questions about MTEX functions, please check out the comprehensive MTEX documentation (https://mtex-toolbox.github.io/Documentation.html)
