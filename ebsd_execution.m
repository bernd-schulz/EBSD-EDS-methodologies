%     Advanced analytical electron microscopy methodologies to characterise microstructural features in superalloys
%     Copyright (C) 2022 B. Schulz et al.
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% IMPORTANT INFORMATION
% Please read README.md before executing this file 
% This file gives an overview of some basic methods, please type following
% to get help with any method:
help ebsd_analysis

% More advanced methodologies, that are also described in the paper [B. Schulz et al. (2023)], can be
% found in the Tutorials folder

%% Reset Matlab
clear;
clc;
close all;

%% Specify Crystal and Specimen Symmetries
% Crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Ni-superalloy', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [4.6 4.6 4.6], 'mineral', 'TiC', 'color', [0.56 0.74 0.56])  
  };

%% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify File Names
% Path to files
pname = [pwd '\Tutorials\3-2-1_RxClassification']; %pwd identifies current folder, but file could be anywhere on PC
% Specify files
fname = [pname '\' 'ebsd_3-2-1_Fig4.ang'];
% Name of File in Report
name = 'RxClassification';

%% Import the Data
% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ang',...
            'convertEuler2SpatialReferenceFrame','setting 2');

%% Script Parameter
minSize = 5; %minGrainSize pixel
minAngle = [1 5]; %[LA-GB HA-GB] degree

%% Create new Object
analysis = ebsd_analysis(ebsd, name)

%% Grain Reconstruction
createGrains(analysis, minAngle, minSize)

% the help command outputs information about the methods, e.g.:
help createGrains

%% Denoise orientation data
cleanData(analysis)

%% Merge Twins
% Merge Rx & deformed grains just within their respective partitions
mergeTwins(analysis, 'Ni-superalloy', 0) %merge all CSL3 grains

%% General plotting functions
help ebsd_analysis/plot

plot(analysis, '', 2,1,0) %IQ map with random GBs & twins

plot(analysis, 'phase', 2,0,0)
plot(analysis, 'phase', 2,0,analysis.Grains) % also shows GBs when grains were merged along CSL3

plot(analysis, 'ipf', 2,0,0) % IPF in direction normal to screen
plot(analysis, 'ipfY', 2,0,0) % IPF in deformation direction (up/down)


%% Boundary Char Analysis between two phases
boundaryChar(analysis, analysis.Grains, 'Ni-superalloy', 'Ni-superalloy')

%% CSL 3, 9 & 27 Output
analyseCSL(analysis, 'Ni-superalloy', analysis.Grains)

%% General grain size analysis 
grainSize(analysis, 'Ni-superalloy')

%% Grain size analysis for Rx/deformed regions
recrystallisedGrains(analysis) % seperates Rx and deformed grains

analyseRx(analysis,'Ni-superalloy')

% saved in analysis.Report
display (analysis.Report)

%% Additionalk plots for deformation microstructures
plot(analysis, '', 3, 1,0) %IQ map with random GBs & twins & subgrains
plot(analysis, 'SubgrainDensity', 2,0,0)

analysis.calculateKAM();
plot(analysis, 'KAM', 2,0,0)
plot(analysis, 'GROD', 2,0,0) 
plot(analysis, 'GOS', 2,0,0)
plot(analysis, 'GAM', 2,0,0)

% analysis.calculateGND() % GND Calculation (high memory usage!)
% plot(analysis, 'GND', 0,0,0)

% Plotting limited to specific regions:
plot(analysis, 'ipfY', 2,0, analysis.RxGrains)
plot(analysis, 'ipfY', 2,0, analysis.DeformedGrains)

%% Documentation
% All steps have been saved in analysis.Report
display(analysis.Report)

%% Save figures, matlab files & report
return %stop script from automatically execuding save

% Will close all open figures!
analysis.save(pname) %pname = save to current path