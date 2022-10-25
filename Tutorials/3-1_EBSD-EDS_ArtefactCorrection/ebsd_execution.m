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

%% Reset Matlab
clear;
clc;
close all;

%% Add folder that includes ebsd_analysis.m to path
% Needed as ebsd_analysis.m is not in same directory as ebsd_execution.m
addpath('../..')

%% Specify Crystal and Specimen Symmetries
% Crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Ni-superalloy', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Ni3Al', 'color', [0.85 0.65 0.13]),...
  crystalSymmetry('m-3m', [4.6 4.6 4.6], 'mineral', 'TiC', 'color', [0.56 0.74 0.56]),...
  };

%% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify File Names
% Path to files
pname = pwd; %pwd identifies current folder, but file could be anywhere on PC
% Specify files
fname = [pname '\' 'ebsd_3-1_Fig1.ctf'];
% Name of File in Report
name = 'EBSD-EDS_ArtefactCorrection';

%% Import the Data
% create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
    'convertEuler2SpatialReferenceFrame');

%% Script Parameter
minSize = 5; %minGrainSize pixel
minAngle = [1 5]; %[LA-GB HA-GB] degree

%% Create new Object
analysis = ebsd_analysis(ebsd, name)

%% Grain Reconstruction
createGrains(analysis, minAngle, minSize)

%% Denoise orientation data
cleanData(analysis)

%% Plot before artefact correction
plot(analysis, '', 0,1,0) %(a) band contrast (BC) map
plot(analysis, 'phase',1,1,0) %(c) Phase with grain boundaries (GBs) & BC in background
plot(analysis, 'ipf', 1,1,0) %(d) IPF with GBs & BC in background

%% Apply y' artefact correction
% Reset denoising to apply correction to as-recorded data
analysis = ebsd_analysis(ebsd, name)
createGrains(analysis, minAngle, minSize)

% Artefact correction properties
datapointLimit = 150; %for 0.2um step size map
pixelThreshold = 1; %in degrees

%First phase: the one that gets removed
cleanEDS(analysis,'Ni-superalloy', datapointLimit, pixelThreshold, 'Ni3Al'); % (e) y' corrected phase map

% Denoise orientation data (after correction)
cleanData(analysis)

plot(analysis, 'ipf', 1,0,0) %(f) IPF with GBs after correction

%% Apply y artefact correction
% Reset denoising to apply correction to as-recorded data
analysis = ebsd_analysis(ebsd, name)
createGrains(analysis, minAngle, minSize)
cleanEDS(analysis,'Ni-superalloy', datapointLimit, pixelThreshold, 'Ni3Al'); % Re-apply the previous y' correction because ebsd_analysis was reseted
close()

% Artefact correction properties
datapointLimit = 15;

%First phase: the one that gets removed
cleanEDS(analysis,'Ni3Al', datapointLimit, pixelThreshold, 'Ni-superalloy'); %(e) y' corrected phase map

%% Remove secondary y'
datapointLimit = 150; %higher limit to remove secondary y'

cleanEDS(analysis,'Ni3Al', datapointLimit, pixelThreshold, 'Ni-superalloy') % (h) phase map after removing secondary y'

% Denoise orientation data (after all corrections)
cleanData(analysis)

%% Documentation
% All steps have been saved in analysis.Report
display(analysis.Report)

%% Save figures, matlab files & report
return %stop script from automatically execuding save

% Will close all open figures!
copyfile ../../ebsd_analysis.m %ebsd_analysis.m must be in same directory as ebsd_execution.m to save it
analysis.save(pname) %pname = save to current path