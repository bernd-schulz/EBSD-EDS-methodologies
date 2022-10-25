# EBSD-EDS-methodologies
## GENERAL
This data is a companion to the manuscript 'Advancing analytical electron microscopy methodologies to characterise microstructural features in superalloys' submitted to Ultramicroscopy by B. Schulz, N. Haghdadi, T. Leitner, M. Hafok and S. Primig [1].

It includes algorithms that were developed to automatically analyse large electron backscatter diffraction (EBSD) 
datasets using the MTEX toolbox for Matlab. Further, sample datasets and tutorials are provided to illustrate how to apply the developed methodologies. The developed methodologies allow to:
* 3-1) Correct γ/γ' interface artefacts in combined EBSD/Electron dispersive X-ray spectroscopy measurements
* 3-2-1) Classify recrystallised grains based on their location
* 3-2-2) Identify recrystallisation nucleation mechanism at distorted twins
* 3-3-1) Describe the evolution of twin boundaries during hot compression
* 3-3-2) Visualise the coherence and grain boundary planes of Σ3 twins

The numbering of the tutorials corresponds with the section numbers in the paper.

Explanation of the underlying methodologies are published here:
[1] B. Schulz, N. Haghdadi, T. Leitner, M. Hafok, S. Primig, Advancing analytical electron microscopy methodologies to characterise microstructural features in superalloys, Ultramicroscopy. submitted for review (2022).

This work © 2022 by B. Schulz et al. is licensed under a [GNU General Public License v3.0 (GPL-3.0-or-later)](https://www.gnu.org/licenses/gpl-3.0-standalone.html). To attribute the material to the creators, please cite the above mentioned paper B. Schulz et al. [1].

## INSTALLATION
1) Install Matlab 2019b
   (newer versions should work but have not been tested)
2) Download mtex-5.7.0 (https://mtex-toolbox.github.io/download)
   (newer versions should work but have not been tested)
3) Start Matlab & change current folder to MTEX folder
4) Type 'startup_mtex' into commend window
5) Press 'install mtex' on output
6) Change folder to the folder where ebsd_analysis.m & ebsd_execution.m are located
7) Execute ebsd_execution.m

## Important Information
Generally, both these files must be located in the current Matlab folder:
* ebsd_analysis.m: Matlab class that includes all methods to analyse the dataset
* ebsd_execution.m: Executable file that calls the methods
Folder structure:
* The ebsd_execution.m in the current folder includes general explanations to the syntax of some basic methods (e.g. plotting)
* The tutorials folder includes instructions on how to recreate the analyses presented in B. Schulz et al. (2022)
Output:
* Most functions do not directly output into the command window but save the output into analysis.Report

The files include .ang or .ctf EBSD datasets or .m Matlab files.

Important: The here developed code has just been tested on Ni-base superalloys. For other systems changes 
in the ebsd_analysis.m class may be necessary and bugs may appear. If any issues arise, please refer to [1] to contact the corresponding author.

If there is general questions about MTEX functions, please refer to the comprehensive MTEX documentation (https://mtex-toolbox.github.io/Documentation.html)
