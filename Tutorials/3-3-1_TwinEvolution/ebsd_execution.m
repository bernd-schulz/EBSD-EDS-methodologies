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
fname = [pname '\' 'ebsd_3-2-1_Fig4.ang']; % different dataset than used in paper!
% Name of File in Report
name = 'RxAtDistortedTwins';

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

%% Seperate Rx & deformed grains
recrystallisedGrains(analysis); %GOS 2 degrees used

%% Merge Twins
% Merge Rx & deformed grains just within their respective partitions
mergeTwins(analysis, 'Ni-superalloy', analysis.DeformedGrains) %merge twins within deformed regions
mergeTwins(analysis, 'Ni-superalloy', analysis.RxGrains) %merge twins within Rx regions

%% Plot IPF in deformation direction
plot(analysis, 'ipfY', 2,0,0)

%% Analyse distorted twins

twinRotationAxis(analysis, 'Ni-superalloy')
% Output:
%   Figure 2: Phase map with distorted twins coloured in red and CSL3 in
%             white
%   Figure 3: Phase map that shows deviation of distorted twins from 60°
%   Figure 4: boundary misorientation angle distribution for GBs within the respective Rx (black) and deformed regions (blue)
%   Figure 5: (10 (c)) boundary misorientation angle distribution for distorted twins (red) and ?3 boundaries between deformed grains (blue)
%   Figure 6: (10 (d)) fundamental sector projection of the misorientations for distorted twins (red) and ?3 boundaries between deformed grains (blue)
%   Figure 7: boundary misorientation angle distribution for CSL3 boundaries within the respective Rx (black) and deformed regions (blue)

%% Documentation
% All steps have been saved in analysis.Report
display(analysis.Report)

%% Save figures, matlab files & report
return %stop script from automatically execuding save

% Will close all open figures!
copyfile ../../ebsd_analysis.m %ebsd_analysis.m must be in same directory as ebsd_execution.m to save it
analysis.save(pname) %pname = save to current path