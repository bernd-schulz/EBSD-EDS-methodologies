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

%% Analyse DRX grains that nucleated at distorted twins

%First: DRX grains must classified
classifyDRX(analysis, 'TiC', 'Ni-superalloy');
close all; % close output from classifyDRX (not of interest for this tutorial, see 3-2-1)

%Second: Closer analysis of DRX grains that nucleated at distorted twins
[results,misori,gos] = boundaryAnalysisRxOnTwins(analysis, 'Ni-superalloy');
% Output: very small map to show code syntax & thus, not statistical
%         significant results!
%
%   Figure 1: (8 (a)) misorientation angle distribution for high angle
%             (blue), low angle (yellow) GBs and the mean misorienation between both
%              distorted twins at the position of the DRX grain (red)
%   Figure 2: (8 (b)) minimum angle between the closest <111> crystal
%              directions for high angle (blue), low angle (yellow) GBs 
%              and the mean misorienation between both distorted twins 
%              at the position of the DRX grain (red)
%   Figure 3: Boundary planes for CSL3 boundaries between these DRX grains
%             (left) and deformed (right) grains terminate on
%
%   results: 
           %Rows: Twin-Twin | GB1 | GB2
           %Columns: CSL3 fraction | Misori mean | Misori std | Misori std
           %| <111> deviation mean | <111> dev std 
           %| point misori between twin/GB mean | std
           %| percentage of CSL3 which are coherent
%    misorientation: Rx grain id | misorienation def1/def2 | ...
%                    stdv def1/def2 | misori def1/Rx | 
%                    stdv def1/Rx | misori def2/Rx | 
%                    stdv def2/Rx
%    gos: gos def | gos def | gos Rx

%% Documentation
% All steps have been saved in analysis.Report
display(analysis.Report)

%% Save figures, matlab files & report
return %stop script from automatically execuding save

% Will close all open figures!
copyfile ../../ebsd_analysis.m %ebsd_analysis.m must be in same directory as ebsd_execution.m to save it
analysis.save(pname) %pname = save to current path