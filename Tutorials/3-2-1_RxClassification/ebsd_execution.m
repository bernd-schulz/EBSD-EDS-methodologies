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

%% Denoise orientation data
cleanData(analysis)

%% Seperate Rx & deformed grains
recrystallisedGrains(analysis); %GOS 2 degrees used

%% Merge Twins
% Merge Rx & deformed grains just within their respective partitions
mergeTwins(analysis, 'Ni-superalloy', analysis.DeformedGrains) %just show twins within deformed region
%mergeTwins(analysis, 'Ni-superalloy', analysis.RxGrains) %if twins within Rx region are also of interest

%% Plot IPF in deformation direction
plot(analysis, 'ipfY', 2,0,0) % (a) IPF with random GBs & CSL3

%% Classification of DRX grains
%first phase: particle for PSN (if PSN is not of interest or single phase
%             write '' instead)
%second phase: phase that recrystalised
classifyDRX(analysis, 'TiC', 'Ni-superalloy')

% Main Output 
%   Figure 2: (b) classified DRX grains
% Additional Outputs (to show that it is possible to further analyse these
%                     data)
%   Figure 3: Misorientation angle between TiC and recrystallised
%             Ni-superalloy grains(blue) and deformed Ni-superalloy grains (red)
%   Figure 4: Misorientation Axes for TiC and recrystallised
%             Ni-superalloy grains(blue) and deformed Ni-superalloy grains (red)
%             --> Discard as makes little sense to plot
%   Figure 5: Diameter of TiC that lead to PSN (blue) vs. did not lead to
%             PSN (red)
%   Figure 6: Size of DRX grains that are classified as 'DRX at MC' (blue) vs.
%             remaining DRX grains (red)
%   Figure 7: Grain orientation spread (GOS) of DRX grains that are classified as 
%             'DRX at MC' (blue) vs. remaining DRX grains (red)

%% Further Outputs
analysis.boundaryAnalysisRxFractions('Ni-superalloy')

%   Figure 8: Phase map that shows CSL3 boundaries between classified Rx
%             grains & deformed grains
%   Figure 9: Misorientation angle histogram between classified Rx & deformed grains
%             (see Figure 5c/5f in paper)
%   Figure 10: Misorientation axes between deformed grains & DRX at
%              distorted twins
%   Figure 11: Misorientation axes between deformed grains & DRX at
%              random GBs
%   Figure 12: Misorientation axes between deformed grains & DRX at
%              random GBs (CSL3 twin segments removed)
%   Figure 13: Misorientation axes between DRX at random GBs & DRX at
%              random GBs (CSL3 twin segments removed)
%   Figure 14: Misorientation angle histogram between deformed grains & DRX
%              at random GBs (blue); 
%              DRX grains at random GBs & DRX grains at random GBs (red)

%% Documentation
% All steps have been saved in analysis.Report
display(analysis.Report)

%% Save figures, matlab files & report
return %stop script from automatically execuding save

% Will close all open figures!
copyfile ../../ebsd_analysis.m %ebsd_analysis.m must be in same directory as ebsd_execution.m to save it
analysis.save(pname) %pname = save to current path