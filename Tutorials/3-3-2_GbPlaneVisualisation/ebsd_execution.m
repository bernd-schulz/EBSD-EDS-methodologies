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
  crystalSymmetry('m-3m', [4.6 4.6 4.6], 'mineral', 'TiC', 'color', [0.56 0.74 0.56])  
  };

%% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%% Specify File Names
% Path to files
pname = pwd; %pwd identifies current folder, but file could be anywhere on PC
% Specify files
fname = [pname '\' 'ebsd_3-3-2_Fig12-14.ang']; % different dataset than used in paper!
% Name of File in Report
name = 'GBPlaneVisulaisation';

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

%% Denoise orientation data
cleanData(analysis)
 
% %% Merge Twins
% % Merge Rx & deformed grains just within their respective partitions
% mergeTwins(analysis, 'Ni-superalloy', 0) %merge twins within deformed regions
% mergeTwins(analysis, 'Ni-superalloy', analysis.RxGrains) %merge twins within Rx regions

%% Plot phase boundary map with coherent/incoherent CSL3 boundaries (in paper Figure 12)
% method applys extra 95 GB smoothing steps & outputs boundary map

analyseCoherence(analysis, 'Ni-superalloy') % automatically checks entire dataset

%% Plot coherent/incoherent CSL3 of selected grains (in paper Figure 11)
% Select grains of interest
grains = analysis.Grains('Ni-superalloy').smooth(95); % important!!!
g1 = grains('id', 361);
g2 = grains('id', 351);

% Select GB between these grains
gb = g1.boundary(g2);

% Check coherence:
[coherent, incoherent, noCSL3] = analysis.checkCoherence('Ni-superalloy', gb, 10, 0)
% Output
%   coherent: all CSL3 that are coherent
%   incoherent: all CSL3 that are incoherent
%   noCSL3: GBs that are not CSL3

% Plotting the coherent/incoherent CSL3 for the GB between the selected grains
figure
plot(g1)
hold on
plot(g2)
hold on
plot(noCSL3, 'linecolor', 'black', 'lineWidth', 1)
plot(coherent, 'linecolor', 'white', 'lineWidth', 1.5, 'DisplayName', 'coherent CSL3')
plot(incoherent, 'linecolor', 'red', 'lineWidth', 1.5, 'DisplayName', 'incoherent CSL3')
hold off

%% Plot GB planes of selected grains (in paper Figure 13 & 14)

gbdirectionPlanetraceAngle(analysis,'Ni-superalloy', g1, g2, 'incoherent', 1) % to change the planes of interest m, please change them in the method in ebsd_analysis.m
% type 'help gbdirectionPlanetraceAngle' to get information on input
% important: g1/g2 must be smoothed BEFORE input into this method
%
% Output
%   Figure 3 & 7: GB planes with lowest misorientation for the coloured grain
%   Figure 4-6 & 8-10: Angle of the traces of selected GB planes and the GB
%                      segments for {112}, {122} and {110} for the coloured grain


%% Documentation
% All steps have been saved in analysis.Report
display(analysis.Report)

%% Save figures, matlab files & report
return %stop script from automatically execuding save

% Will close all open figures!
copyfile ../../ebsd_analysis.m %ebsd_analysis.m must be in same directory as ebsd_execution.m to save it
analysis.save(pname) %pname = save to current path