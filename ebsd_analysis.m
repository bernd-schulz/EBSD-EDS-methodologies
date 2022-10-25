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

classdef ebsd_analysis < handle
    %EBSD Clean-UP & Analyse Class (based on MTEX)
    %MTEX version as mentioned below must be installed!
    %Use in Function order:
    %   1. MTEX function: import ebsd
    %   2. Create object ebsd_analysis(....)
    %   3. createGrains(....)
    %   4. optional: cleanEDS(...)
    %   5. optional: cleanData(....)
    %   6. mergeTwins(....)
    %   7. *OTHER FUNCTIONS*
    %   9. save(...)
    %
    % Refer to B. Schulz et al. (2023) for an overview of the available
    % functions
    
    %% CHANGE LOG
    %v1 25/05/2020: added EBSD/EDS Correction
    %v2 17/11/2020: Update to MTEX 5.5.beta.6; add DRX Analysis
    %v3 05/01/2020: Update to MTEX 5.6.beta.2;
    %v4 25/03/2020: Update to MTEX 5.6.1; change EDAX clean-up (grain
    %               CI-Standardisation + CI>0.1 filter!)
    %v5 03/06/2021: Update to MTEX 5.7.0; clean-up class; new folder saving
    %               structure
    %   20/06/2021: Standardised deformation classes output
    %v6 05/07/2021: Corrected error in grain ID assignment after twinning!
    %               (should have just influenced grain-partition plotting
    %v6.1 12/07/2021: distinguish between grains on prior twin boundaries; do not merge deformed/Rx Twins 
    %     23/07/2021: Added Twin/Rotation Analysis; plot TraceSlipPlane;
    %                 PoleFigure;
    %     29/07/2021: Added Rx-Boundary Char Analysis
    %v6.2 24/09/2021: Added function to analyse Rx-grains on prior-twins
    %     10/03/2022: Fixed issue with IPF plotting in plotGbdirectionPlanetraceAngle
    %v7   06/09/2022: Cleaned up code & added comments
    %v8   13/09/2022: Changed some function arguments
    %                 obj.analysePSN renamed to obj.classifyDRX

    %% Class Properties
    properties
        ScriptVersion = 'v8 (13/09/2022)'; %Version of class & last updated
        Name = 'EBSD Analysis'; %Name of EBSD Scan (placeholder)
        vMTEX = 'MATLAB R2019b, MTEX 5.7.0';
        EbsdOriginal; %Original EBSD Data 
        Ebsd; %Altered EBSD Data (e.g. denoising/EDS Clean-Up)
        EbsdOriginalGrains; %EBSD data before merged grains but after denoising
        Grains; %Contains grains information
        minGrainAngle; %min. angle in degrees to detect grains [LA-GB HA-GB]
        minGrainSize; %min. pixel to detect grains
        Twins; %Contains all twins information
        TwinsDDRX; %Contains twins between Rx & deformed fraction
        MergedGrains; %Contains information of grains where twins are removed
        RxThreshold = 2; %GOS threshold in degrees for Rx/deformed grains
        RxGrains; %Recrystallised (Rx) Grains
        DeformedGrains; %Deformed Grains
        PSNGrains; %DRX at MC carbides
        DRXGrains; %DRX grains (same as RxGrains without PSN)
        DRXGrainsWithinGrain; %DRX inside grains
        DRXGrainsOnGBs; %DRX at random GBs
        DRXGrainsOnPriorTwinGBs; %DRX at distorted twins
        Report; %Stores all processes which have been applied to data
        GrainStatistics; %Stores the diameter of all analysed grains
        minGrainDiameter = 0; %Stores the minGrainDiameter if put into system, which limits grain statistics
        KAM; %Kernal average misorienation (KAM) in rad
        GAM; %Grain average misorientation (GAM) in rad
        GND; %Geometrically necessary dislocations (GNDs) in m^(-2)
        LineProfileCount = 0; %Number of plotted line profiles (for labeling)
        CoherentCSL3Threshold = 10; %max. deviation of GB trace & (111) trace in degrees
        GrainCIStandardisation = 0.1; %Confidence Index (CI) standardisation threshold
    end
    
    methods
        function obj = ebsd_analysis(ebsd, name)
            %% Constructs a class ebsd_analysis with
            % INPUT
            %   ebsd: EBSD data from MTEX with import_wizard(ebsd)!
            %         Check reference convention is correct
            %   name: Name of EBSD scan for Report/Saving Plots with
            %         save-function
            %
            % OUTPUT
            %   sets obj.Ebsd, obj.EbsdOriginal, obj.Report
            
            %% Save ebsd data
            % Saves EBSD scan in property EbsdOriginal
            obj.EbsdOriginal = ebsd;
                        
            % Extracts just indexed points from ebsd
            obj.Ebsd = ebsd('indexed');
            
            %% Add name of scan to report
            %Save used MTEX & Matlab Version to report
            obj.Report = 'EBSD Analaysis Script ('+ string(obj.vMTEX) +'): ' + string(obj.ScriptVersion);
            
            %Save EBSD map name to report
            obj.writeReport('Specimen', name)
          
            %Save EBSD map name to class properties
            name = regexprep(name, '\s+', '');
            obj.Name = name;            
        end
         
        function createGrains(obj, minAngle, minGrainSize)
            %% Create Grains
            % INPUT
            %   minAngle: min. misorienation angle difference to detect grains
            %   minGrainSize: min. pixel to detect grains
            %
            % OUTPUT
            %   sets obj.Grains, obj.Ebsd.grainId
            
            % Save as class properties
            obj.minGrainAngle = minAngle;
            obj.minGrainSize = minGrainSize;
            
            % set boundary angle limit            
            [obj.Grains,obj.Ebsd.grainId,obj.Ebsd.mis2mean] = calcGrains(obj.Ebsd('indexed'),'angle',[minAngle(1)*degree, minAngle(2)*degree]);
            % remove all grains smaller than pixel limit
            obj.Ebsd(obj.Grains(obj.Grains.grainSize <= minGrainSize)) = [];
            % set boundary angle limit after removing small grains
            [obj.Grains,obj.Ebsd.grainId] = calcGrains(obj.Ebsd('indexed'),'angle',[minAngle(1)*degree, minAngle(2)*degree]);
            
            %% Just consider grains with CI>0.1 after Grain Standardisation (17/03/21)
            try
                %% Edax System            	
                %% Grain CI Standardisation & low CI removal
                % Grain CI Standardisation
                obj.Grains.prop.maxCI = grainMean(obj.Ebsd, obj.Ebsd.ci, obj.Grains, @max);
                % Remove Grains with CI < 0.1
                obj.Ebsd(obj.Grains(obj.Grains.prop.maxCI <= obj.GrainCIStandardisation)) = [];
                % Recalculate grains
                [obj.Grains,obj.Ebsd.grainId] = calcGrains(obj.Ebsd('indexed'),'angle', [minAngle(1)*degree, minAngle(2)*degree]);
                
                obj.writeReport('Grain Standardisation', ['Just grains with max(CI)>' num2str(obj.GrainCIStandardisation) ' considered'])        
            catch
                %% Oxford System
                % IS THERE A PARAMTER TO DO THIS? No equivalent known!
            end
                        
            %% smooth boundaries
            obj.Grains = smooth(obj.Grains,5); % 5 iterations
            
            %% Save original EBSD data
            % this property will not be overwritten at any stage (even if
            % grains are merged)
            obj.EbsdOriginalGrains = obj.Ebsd;

            % Save grainIDs of original grains to GBs (for analysis which require un-merged grainIds)
            obj.Grains.prop.originalId = obj.Grains.id;
            
            %% Save Inclusions
            obj.Grains.prop.InclusionBeforeMerge = isInclusion(obj.Grains);

            %% Add to report
            obj.writeReport('Grains Calculated', '')
            obj.writeReport('Min. Grain Size', minGrainSize)
            obj.writeReport('Min. Grain Angle', minAngle(2))
            obj.writeReport('Min. Subgrain Angle', minAngle(1))
            obj.writeReport('-', '-');
                        
            %% Seperate recrystallised & deformed grains
            obj.recrystallisedGrains(); % threshold defined in property obj.RxThreshold
        end
        
        function createGrainsNotIndexed(obj, minAngle, minGrainSize)
            %% Create Grains where notIndexed points are treated as own grains (e.g. for cracks)
            % INPUT
            %   minAngle: min. misorienation angle difference to detect grains
            %   minGrainSize: min. pixel to detect grains
            %
            % OUTPUT
            %   sets obj.Grains, obj.Ebsd.grainId
            
            obj.Ebsd = obj.EbsdOriginal; %
            
            % Save as class properties
            obj.minGrainAngle = minAngle;
            obj.minGrainSize = minGrainSize;
            
            obj.Ebsd(obj.Ebsd.ci<0.1).phase = -1; %Assign all spots with low CI to notIndexed (= Crack)

            % set boundary angle limit            
            [obj.Grains,obj.Ebsd.grainId,obj.Ebsd.mis2mean] = calcGrains(obj.Ebsd,'angle',[minAngle(1)*degree, minAngle(2)*degree]);
            % remove all grains smaller than pixel limit
            obj.Ebsd(obj.Grains(obj.Grains.grainSize <= minGrainSize)) = [];
            % set boundary angle limit after removing small grains
            [obj.Grains,obj.Ebsd.grainId] = calcGrains(obj.Ebsd,'angle',[minAngle(1)*degree, minAngle(2)*degree]);
            
            %% Just consider grains with CI>0.1 after Grain Standardisation (17/03/21)
            try
                %% Edax System            	
                %% Grain CI Standardisation & low CI removal
                % Grain CI Standardisation
                obj.Grains.prop.maxCI = grainMean(obj.Ebsd, obj.Ebsd.ci, obj.Grains, @max);
                % Remove Grains with CI < 0.1
                g = obj.Grains('indexed')
                obj.Ebsd(g(g.prop.maxCI <= 0.1)) = [];
                [obj.Grains,obj.Ebsd.grainId] = calcGrains(obj.Ebsd,'angle', [minAngle(1)*degree, minAngle(2)*degree]);
                
                obj.writeReport('Grain Standardisation', 'Just grains with max(CI)>0.1 considered')        
            catch
                % Oxford System
                % IS THERE A PARAMTER TO DO THIS?
            end
                        
            %% smooth boundaries
            obj.Grains = smooth(obj.Grains,5);
            obj.EbsdOriginalGrains = obj.Ebsd;  % Safe in EBSD Original (in case if obj.Ebsd is merged later)
            
            %% Save grainIDs of original grains to GBs (for analysis which require un-merged grainIds)
            obj.Grains.prop.originalId = obj.Grains.id;
            
            %% Save Inclusions
            obj.Grains.prop.InclusionBeforeMerge = isInclusion(obj.Grains);

            % add to report
            obj.writeReport('Grains Calculated', '')
            obj.writeReport('Min. Grain Size', minGrainSize)
            obj.writeReport('Min. Grain Angle', minAngle(2))
            obj.writeReport('Min. Subgrain Angle', minAngle(1))
            obj.writeReport('-', '-');
            
            obj.recrystallisedGrains();
        end
        
        function resetGrains(obj)
            %% Resets all grain & ebsd information
      
            % Resets EBSD to originally imported data
            obj.Ebsd = obj.EbsdOriginalGrains;
            
            % Delete all grain information
            obj.Grains = obj.Grains;
            obj.Twins = [];
            obj.TwinsDDRX = [];
            obj.MergedGrains = [];
            obj.RxGrains = [];
            obj.DeformedGrains = [];
            obj.PSNGrains = [];
            obj.DRXGrains = [];
            obj.DRXGrainsWithinGrain = [];
            obj.DRXGrainsOnGBs = [];
            obj.DRXGrainsOnPriorTwinGBs = [];
                        
            %% Write report
            obj.writeReport('Grain rest', 'all grain data reseted to as-created')
        end
        
        function cleanData(obj) 
            %% Cleaning functions
            %   Cleans EBSD orientation data using the halfQuadraticFiler
            %   with alpha 0.25
            %   Applied after grains are created, thus does not affect
            %   grain size distribution
            
            % apply filter to indexed EBSD points and fill grains except GB
            % half quadratic filter mantains subgrain structures
            F = halfQuadraticFilter;
            F.alpha = 0.25;
            obj.Ebsd = smooth(obj.Ebsd('indexed'),'fill',F);
            obj.Ebsd = obj.Ebsd('indexed');
            
            obj.EbsdOriginalGrains = obj.Ebsd;
            % Line above can include 'grains' argument to not fill GBs
            % ebsd = smooth(ebsd_original('indexed'),'fill',F,grains);
            
            % add to report
            obj.writeReport('Half Quadratic Denoising with alpha', F.alpha)
            obj.writeReport('-', '-');
        end
        
        function mergeTwins(obj, phase, gPartition)
            %% Detect and merge twins
            % INPUT
            %   phase: name of phase which should be merged
            %   g: grain fraction which should be merged (usually
            %      recrystallised/deformed)
            %      put 0 if it should be applied to all grains
            %
            % OUTPUT
            % sets obj.MergedGrains, obj.TwinsDDRX
            
            brandon = 8.66; %Brandon criterion (deviation from CSL3 that is still considered a twin)
            
            %saves grains properties in local variable in case entire
            if (isempty(obj.MergedGrains) == 1)
                g = obj.Grains;
            else
                g = obj.MergedGrains; 
            end
                
            if strcmp(class(gPartition), 'grain2d')
                %in case just grains within defined g should be merged
                partition = true;
                
                %Reasign name:
                gAll = g;
                %Write property to all grains which should be merged:
                gAll(gPartition.id).prop.MergedFraction = 1;
                
                %Twins within this fraction should be merged:
                g = gAll(gAll.prop.MergedFraction == 1);
                
                %Remaining grains which should not be merged:
                gRest = gAll(gAll.prop.MergedFraction == 0);
            else
                partition = false;
            end
            
            % Find all GB between Phase
            gB = g(phase).boundary(phase,phase);
            
            %% Remove Twins
            if (partition == true)
                %Delete Twins between g-fraction and other fractions
                
                %GBs on outside of partition:
                gBOutline = gB(g(phase).boundary(phase,phase).hasGrain(gRest));
                %GB just within partition:
                gB = gB(~g(phase).boundary(phase,phase).hasGrain(gRest));
                
                %Twin boundaries between partition & outside
                gB3Outline = gBOutline(angle(gBOutline.misorientation,CSL(3,obj.Ebsd(phase).CS)) < brandon*degree);
                
                if (isempty(obj.TwinsDDRX )) %obj.TwinsDDRX includes twin boundaries between two partitions (not merged along these twins)
                    obj.TwinsDDRX = gB3Outline;
                else
                    obj.TwinsDDRX = [obj.TwinsDDRX gB3Outline];
                end
                
                %% saves grains properties in local variable
                if (isempty(obj.MergedGrains) == 1)
                    g = obj.Grains;
                else
                    g = obj.MergedGrains; 
                end                
            end
            
            % select CSL(3) grain boundaries
            if (isempty(gB))
                gB3 = [];
            else
                gB3 = gB(angle(gB.misorientation,CSL(3,obj.Ebsd(phase).CS)) < brandon*degree); %8.66 accoridng to Brendon Criterion (often also 3 degrees used)
            end
            
            %% Save if inclusion prior to merging
            g.prop.InclusionBeforeMerge = double(g.prop.InclusionBeforeMerge);
            
            if (~isempty(gB3))
                % checks if twins of another phase have already been detect
                % before
                if (isempty(obj.Twins) == 1)
                    obj.Twins = gB3;
                else
                    obj.Twins = [obj.Twins gB3];
                end

                %% Merge Twinned Grains             
                %Assign new grains consistent grainIds  
                [merged,parentId] = merge(g,gB3);  
                obj.Ebsd('indexed').grainId = parentId(g.id2ind(obj.Ebsd('indexed').grainId));  

                %% Re-asign property to identify inclusions prior to Merging
                % Could cause problems if inclusion has been merged along twin!
                merged.prop.InclusionBeforeMerge(isnan(merged.prop.InclusionBeforeMerge)) = 0;
                merged.prop.InclusionBeforeMerge = logical( merged.prop.InclusionBeforeMerge);
               
                %% SUM over GOS (Twins area-weighted)
                % extract GOS and area
                grains = g;
                childGOS = grains.GOS;
                childArea = grains.area;

                % compute the weighted averages
                merged.prop.GOS = accumarray(parentId,1:length(grains),size(merged),...
                  @(id) nanmeanWeights(childGOS(id),childArea(id)));

                obj.MergedGrains = merged;
            else
                gB3 = [];
                gB3.segLength = 0;
            end
            
            %% add to report
            if (partition == true)
                obj.writeReport('Twins considered (Partition)', '');
            end
            
            obj.writeReport(['Twins considered (CSL 3 within ' num2str(brandon) ' [Brandon criterion] degrees)'], phase);
            obj.writeReport('Twin-fraction in length.% of total GBs of this phase', sum(gB3.segLength)/sum(g(phase).boundary.segLength)*100);
            obj.writeReport('-', '-');
            
            if (partition == true)
               recrystallisedGrains(obj); 
            end
        end
        
        function analyseCSL(obj, phase, grains)
            %% Looks at 3^n CSL and creates ratio to allow to judge if twins happend before grain growth or after
            %% if before higher ratio of CSL9
            % INPUT
            %   phase: name of phase which should be merged
            %   grains: grain partition that will be investigated
            %
            % OUTPUT
            %   figure
            
            % saves grains properties in local variable
            g = grains;
                       
            %% Detect CSL 3, 9 and 27 within 3 degrees
            % twins in phase
            gB = g(phase).boundary(phase,phase);
            
            % select CSL(3)-CSL(27) grain boundaries
            gB3 = gB(gB.isTwinning(CSL(3,obj.Ebsd(phase).CS), 8.66*degree));
            gB9 = gB(gB.isTwinning(CSL(9,obj.Ebsd(phase).CS), 5.00*degree));
            gB27 = gB(gB.isTwinning(CSL(27,obj.Ebsd(phase).CS), 2.89*degree)); %2.89 (15/sqrt(27))
            
            % Alternative: gB3 = gB(angle(gB.misorientation,CSL(3,obj.Ebsd(phase).CS)) < 8.66*degree);
                       
            %% Plot grains with special GBs
            f = plot(obj, 'phase', 1, 0, 0);
            f.Name = ['Phase_with_CSL_in_ ' phase];
            hold on
            plot(gB27,'lineColor','yellow','linewidth',1,'DisplayName','CSL 27')
            hold on
            plot(gB9,'lineColor','red','linewidth',1,'DisplayName','CSL 9')
            hold on
            plot(gB3,'lineColor','white','linewidth',1,'DisplayName','CSL 3')
            hold on
            %plot(gb3_coherent,'lineColor','green','linewidth',1,'DisplayName','CSL 3 (Coherent)')
            hold off
            
            %% Overall Proportion
            % Lenght Proportion
            lProp3 = sum(gB3.segLength)/sum(g(phase).boundary.segLength)*100;
            lProp9 = sum(gB9.segLength)/sum(g(phase).boundary.segLength)*100;
            lProp27 = sum(gB27.segLength)/sum(g(phase).boundary.segLength)*100;

            % Number Fraction (NOT CORRECT! Different implementation
            % needed)
            nProp3 = length(gB3)/length(g(phase).boundary)*100;
            nProp9 = length(gB9)/length(g(phase).boundary)*100;
            nProp27 = length(gB27)/length(g(phase).boundary)*100;
            
            %% Plot number & length fraction of different CSL boundaries
            figure('Name', ['CSL_NumberFraction_ ' phase]);
            scatter([1,2,3], [nProp3, nProp9, nProp27], 'DisplayName', 'Number Fraction')
            hold on
            scatter([1,2,3], [lProp3, lProp9, lProp27], 'DisplayName', 'Length Fraction')
            legend show
            xlabel(['\Sigman type boundary'])
            ylabel('Fraction [%]')
            set(gca,'xtick',[1,2,3])
            set(gca,'xticklabel',[3,9,27])
            
            %% Write Report
            obj.writeReport('-','-')
            obj.writeReport('CSL ANALYSIS', phase)
            obj.writeReport('Length Proportion CSL3/total [%]', lProp3)
            obj.writeReport('Length Proportion CSL9/total [%]', lProp9)
            obj.writeReport('Length Proportion CSL27/total [%]', lProp27)
            obj.writeReport('Number Proportion CSL3/total [%]', nProp3)
            obj.writeReport('Number Proportion CSL9/total [%]', nProp9)
            obj.writeReport('Number Proportion CSL27/total [%]', nProp27)
        end

        function limitGrains(obj, phase, minD)
            %% Exclude grains below minDiamter
            % INPUT
            %   phase: name of phase of interest
            %   minD: smallest diameter that should be considered
            %
            % OUTPUT
            %   update obj.Ebsd to exclude grains
            
            % Check if Twins have been detected
            if (isempty(obj.MergedGrains) == 1)
                grains = obj.Grains;
            else
                grains = obj.MergedGrains; 
            end
             
            obj.Ebsd(grains(grains.equivalentRadius < minD/2)) = [];
            
            obj.minGrainDiameter = minD;
        end
        
        function recrystallisedGrains(obj)
            %% Seperates DRX (RxGrains) from not-DRX grains (deformedGrains)
            
            threshold = obj.RxThreshold; % in degrees
            
            if (isempty(obj.MergedGrains) == 1)
                grains = obj.Grains;
            else
                grains = obj.MergedGrains;
            end
            
            %% Classification based on GOS
            obj.RxGrains = grains(grains.prop.GOS./degree < threshold);
            obj.DeformedGrains = grains(grains.prop.GOS./degree >= threshold);

            obj.writeReport('Recrystallisation GOS Threshold [°]', threshold)                       
        end
        
        function analyseRx(obj, matrix)
            %% Analyse Rx fraction of matrix phase
            % INPUT
            %   matrix: name of phase of interest
            %
            % OUTPUT
            %   histogram
            
            %% If RX Function has not been applied yet, apply it
            % Seperate Recrystallised/Deformed Grains
            if (isempty(obj.RxGrains) == true && isempty(obj.DeformedGrains) == true)
                obj.recrystallisedGrains();
            end
            
            %% Calculate properties of Rx/deformed partition
            rxAreaFraction = sum(obj.RxGrains(matrix).area)/sum(obj.MergedGrains(matrix).area);
            rxD = mean(abs(obj.RxGrains(matrix).equivalentRadius))*2;
            rxDstd = std(abs(obj.RxGrains(matrix).equivalentRadius))*2;
            rxMedian = median(abs(obj.RxGrains(matrix).equivalentRadius))*2;
            
            deformedD = mean(abs(obj.DeformedGrains(matrix).equivalentRadius))*2;
            deformedDstd = std(abs(obj.DeformedGrains(matrix).equivalentRadius))*2;
            deformedMedian = median(abs(obj.DeformedGrains(matrix).equivalentRadius))*2;
        
            %% Histogram of matrix/matrix GB (before merging)
            rxGrainsWithTwins = obj.Grains(obj.Grains.prop.GOS./degree < obj.RxThreshold);
            deformedGrainsWithTwins = obj.Grains(obj.Grains.prop.GOS./degree >= obj.RxThreshold);

            AngleRx = rxGrainsWithTwins.boundary(matrix,matrix).misorientation.angle./degree;
            try
                AngleDeformed = deformedGrainsWithTwins.boundary(matrix,matrix).misorientation.angle./degree;
            catch
                AngleDeformed = [];
            end
           
            %% Plot in Histogram
            step = 1;
            fBoundaryChar = figure('Name', ['BoundaryCharacter_' matrix '_RxVsDeformed']);
            fBoundaryChar = histogram(AngleRx, 0:step:65, 'Normalization','probability');
            hold on
            fBoundaryChar = histogram(AngleDeformed, 0:step:65, 'Normalization','probability'); 
            legend('Rx Grains', 'Deformed Grains')
            xlabel(['Misorientation [°]'])
            ylabel('Relative Frequency')
            hold off
            
            %% Write Report
            obj.writeReport('-','-')
            obj.writeReport(['Rx Area Fraction ' matrix], rxAreaFraction)
            obj.writeReport('Rx Grains Diameter [um]',rxD)
            obj.writeReport('Rx Grains Diameter Std',rxDstd)
            obj.writeReport('Rx Grains Diameter Median [um]',rxMedian)
            obj.writeReport('GOS Mean [°]', mean(abs(rxGrainsWithTwins.prop.GOS./degree)) )
            obj.writeReport('GOS Std [°]', std(abs(rxGrainsWithTwins.prop.GOS./degree)) )
            
            if (~isnan(deformedD))
                obj.writeReport('Deformed Grains Diameter [um]',deformedD)
                obj.writeReport('Deformed Grains Diameter Std',deformedDstd)
                obj.writeReport('Deformed Grains Diameter Median [um]',deformedMedian)
                obj.writeReport('GOS Mean [°]', mean(abs(deformedGrainsWithTwins.prop.GOS./degree)) )
                obj.writeReport('GOS Std [°]', std(abs(deformedGrainsWithTwins.prop.GOS./degree)) )
            end
        end
        
        function twinRotationAxis(obj, phase)
            %% Stereographic projection of rotation axis of twins
            % INPUT
            %   phase: name of phase of interest
            %
            % OUTPUT
            %   figures
            
            %% Define deformed(Def)/recrystallised(Rex) grains (g) and grain boundaries (gB)
            gDef = obj.Grains(obj.Grains.prop.GOS./degree >= obj.RxThreshold);
            gRex = obj.Grains(obj.Grains.prop.GOS./degree < obj.RxThreshold);
            gBRex = gRex(phase).boundary(phase, phase); 
            gBDef = gDef(phase).boundary(phase,phase);       
            
            %% CSL 3GBs 
            %Rx twins
            sigma = 3;
            % Just consideres Twins between Rx grains, not between Rx &
            % Deformed grains (takes long to calculate!)
            gB3 = gBRex(gBRex.isTwinning(CSL(sigma,obj.Ebsd(phase).CS), (15/sqrt(sigma))*degree));
            gB3Rex = gB3(~gB3.hasGrain(gDef)); %Twins just within rx fraction
            
            % GB Deformed
            % select CSL(3) grain boundaries
            gB3 = gBDef(gBDef.isTwinning(CSL(sigma,obj.Ebsd(phase).CS), (15/sqrt(sigma))*degree));
           
            grainIds = unique(gB3.grainId,'rows');
            twinsDef = [];
            
            wb = waitbar(0, 'Number of CSL3 GBs analysed:');
            
            % The following loop to extract distorted twins could be replaced 
            % with matrix calculation, which would be faster but requires 
            % a lot of memory!
            
            % Extract all distorted twins based on CSL3 between deformed
            % grains
            for i=1:length(grainIds) 
                id = grainIds(i, :);
                
                waitbar(i/length(grainIds), wb)
                
                if isempty(twinsDef)
                    twinsDef = gDef(gDef.id==id(1)).boundary(gDef (gDef.id==id(2)));
                else         
                    twinsDef = [twinsDef gDef(gDef.id==id(1)).boundary(gDef (gDef.id==id(2)))];
                end
            end
            
            close(wb)
            gB3Def = gB3(~gB3.hasGrain(gRex)); %information on twins within deformed grains which still follow twin boundary conditions within Brandon
            
           %% CSL 9GBs 
           %Rx twins
           sigma = 9;
           gB9 = gBRex(gBRex.isTwinning(CSL(sigma,obj.Ebsd(phase).CS), (15/sqrt(sigma))*degree));
           
           if (~isempty(gB9))
                gB9Rex = gB9(~gB9.hasGrain(gDef)); %Twins just within rx fraction
           else
               gB9Rex = gB9;
           end
           
           % GB Deformed
           % select CSL(9) grain boundaries
           gB9 = gBDef(gBDef.isTwinning(CSL(sigma,obj.Ebsd(phase).CS), (15/sqrt(sigma))*degree));

            if(~isempty(gB9))
                gB9Def = gB9(~gB9.hasGrain(gRex)); %information on twins within deformed grains which still follow twin boundary conditions within Brandon
            else
                gB9Def = gB9;
            end
            
          %% CSL 27GBs 
          %Rx grains
          gB27 = gBRex(gBRex.isTwinning(CSL(sigma,obj.Ebsd(phase).CS), (15/sqrt(sigma))*degree));
          
          if(~isempty(gB27))
            gB27Rex = gB27(~gB27.hasGrain(gDef)); %Twins just within rx fraction
          else
              gB27Rex = gB27;
          end
          
          %Deformed grains
          gB27 = gBDef(gBDef.isTwinning(CSL(sigma,obj.Ebsd(phase).CS), (15/sqrt(sigma))*degree));

          if(~isempty(gB27))
              gB27Def = gB27(~gB27.hasGrain(gRex)); %information on twins within deformed grains which still follow twin boundary conditions within Brandon
          else
              gB27Def = gB27;
          end           
            
            %% Plots: Comparison Rx & Deformed grains
            %GB map with distorted twins marked
            f1 = plot(obj, 'phase', 2,0,0)
            f1.Name = 'TwinRot_PriorTwinsVsCSL3_Boundaries';
            hold on
            plot(twinsDef,'LineColor', 'red','linewidth', 2, 'DisplayName','Distorted twins') %= prior-twins
            plot(gB3Def,'LineColor', 'lineColor','white','linewidth',1.5)
            
            %GB map with deviation from CSL3 of distorted twins marked
            f2 = plot(obj, 'phase', 2,0,0)
            f2.Name = 'TwinRot_PriorTwinsVsCSL3_Boundaries_MisoriDev';
            hold on
            if(~isempty(twinsDef))
                plot(twinsDef,((60-abs(twinsDef.misorientation.angle./degree))),'linewidth', 2, 'DisplayName','Distorted twins')
            end
            mtexColorbar('title', 'Misorientation deviation from 60° [°]')
            setColorRange([0, 20])
            mtexColorMap LaboTeX
            
            colorList = obj.getColorList();
            
            %Histogram that shows misorienation angle
            figure('Name', ['TwinRot_DefVsRx_BoundaryCharacter_' phase]);
            edges = 0:2.5:65;
            h1 = histcounts(gBRex(~gBRex.hasGrain(gDef)).misorientation.angle./degree, edges, 'Normalization', 'probability');
            h2 = histcounts(gBDef(~gBDef.hasGrain(gRex)).misorientation.angle./degree, edges, 'Normalization', 'probability');
            b1 = bar(edges(1:end-1),[h1; h2]',1);
            hold on
            [x1,p1] = obj.calcRandomMisorientation(edges(1:end-1), h1); %sum(h1)=sum(h2)=sum(h3)
            hold on
            plot(x1,p1, '--', 'linewidth',1, 'color', 'black')
            legend([{'GBs between recrystallised grains'}; {'GBs between deformed grains'}; {'Random distribution'}])
            xlabel(['Misorientation [°]'])
            ylabel('Relative Frequency [-]')
            legend show
            xlim([-1 edges(end)])
            b1(1).FaceColor = colorList(1,:);
            b1(2).FaceColor = colorList(3,:);

            if (~isempty(gB3Def))
                %Histogram that shows misorienation angle
                figure('Name', ['TwinRot_PriorTwinsVsCSL3_BoundaryCharacter_' phase]);
                edges = 0:2.5:65;
                h1 = histcounts(twinsDef.misorientation.angle/degree, edges);
                h2 = histcounts(gB3Def.misorientation.angle/degree, edges); %, 'Normalization', 'probability' does not work in this case
                b2 = bar(edges(1:end-1),[h1; h2]', 1);
                hold on
                [x1,p1] = obj.calcRandomMisorientation(edges(1:end-1), h1); %sum(h1)=sum(h2)=sum(h3)
                hold on
                plot(x1,p1, '--', 'linewidth',1, 'color', 'black')
                legend([{'Distorted twins between deformed grains'};{'CSL 3 between deformed grains'}; {'Random distribution'}])
                xlabel(['Misorientation [°]'])
                ylabel('Absolute Frequency [-]')
                legend show
                xlim([-1 edges(end)])
                b2(1).FaceColor = colorList(2,:);
                b2(2).FaceColor = colorList(3,:);

                %Plot that shows rotation axes
                f2 = figure('Name', ['TwinRot_PriorTwinsVsCSL3_MisorientationAxis_' phase]);
                hold on
                plot(twinsDef.misorientation.axis, 'fundamentalRegion', 'MarkerColor', colorList(2,:), 'DisplayName', 'Distorted twins between deformed grains')
                hold on
                plot(gB3Def.misorientation.axis, 'fundamentalRegion', 'MarkerColor', colorList(3,:), 'DisplayName', 'CSL 3 between deformed grains')
                legend show
                f2.Name = ['TwinRot_PriorTwinsVsCSL3_MisorientationAxis_' phase];
            end
                
                %% Comprison CSL3 Boundaries
                figure('Name', ['TwinRot_CSL3_RxVsDeformedGrains_BoundaryCharacter_' phase]);
                edges = 50:1:61; %in paper 61 as upper limit used --> maybe saver to use 62 to not miss any!
                h1 = histcounts(gB3Rex.misorientation.angle/degree, edges, 'Normalization', 'probability');
                if(~isempty(gB3Def))
                    h2 = histcounts(gB3Def.misorientation.angle/degree, edges, 'Normalization', 'probability');
                else
                    h2 = histcounts([], edges, 'Normalization', 'probability');
                end
                b3 = bar(edges(1:end-1),[h1; h2]',1);
                hold on
                legend([{'CSL3 between recrystallised grains'}; {'CSL3 between deformed grains'};])
                xlabel(['Misorientation [°]'])
                ylabel('Relative Frequency [-]')
                legend show
                xlim([edges(1) edges(end)])
                b3(1).FaceColor = colorList(1,:);
                b3(2).FaceColor = colorList(3,:);
                            
            %% Write Report: 
            obj.writeReport('-','-')
            obj.writeReport('Twin Rotation Axis Analysis (CSL within fractions)',phase)
            %Twin fraction: Deformed vs Rx grains
            %CSL9 & 27: Within Deformed vs Rx grains
                obj.writeReport('Sigma3 Num. % Recrystallised',gB3Rex.length/gBRex.length*100)
                obj.writeReport('Sigma3 Length % Recrystallised',sum(gB3Rex.segLength)/sum(gBRex.segLength)*100)
                obj.writeReport('Sigma9 Num. % Recrystallised',gB9Rex.length/gBRex.length*100)
                obj.writeReport('Sigma9 Length % Recrystallised',sum(gB9Rex.segLength)/sum(gBRex.segLength)*100)
                obj.writeReport('Sigma27 Num. % Recrystallised',gB27Rex.length/gBRex.length*100)
                obj.writeReport('Sigma27 Length % Recrystallised',sum(gB27Rex.segLength)/sum(gBRex.segLength)*100)

                obj.writeReport('Sigma3 Num. % Deformed',gB3Def.length/gBDef.length*100)
                obj.writeReport('Sigma3 Length % Deformed',sum(gB3Def.segLength)/sum(gBDef.segLength)*100)
                obj.writeReport('Sigma9 Num. % Deformed',gB9Def.length/gBDef.length*100)
                obj.writeReport('Sigma9 Length % Deformed',sum(gB9Def.segLength)/sum(gBDef.segLength)*100)
                obj.writeReport('Sigma27 Num. % Deformed',gB27Def.length/gBDef.length*100)
                obj.writeReport('Sigma27 Length % Deformed',sum(gB27Def.segLength)/sum(gBDef.segLength)*100)

                obj.writeReport('Deformed Sigma3/priorTwin Num. %',gB3Def.length/twinsDef.length*100)
                obj.writeReport('Deformed Sigma3/priorTwin Length %',sum(gB3Def.segLength)/sum(twinsDef.segLength)*100)%Twin fraction: Deformed vs Rx grains
                
            %Twin Density
                obj.writeReport('Twin Density: Deformed Sigma3 [um^-1]',sum(gB3Def.segLength)/sum(gDef(phase).area))
                obj.writeReport('Twin Density: Deformed Prior Twins [um^-1]',sum(twinsDef.segLength)/sum(gDef(phase).area))
                obj.writeReport('Twin Density: Recrystallised Sigma3 [um^-1]',sum(gB3Rex.segLength)/sum(gRex(phase).area))
        end
        
        function c = getColorList(obj)
            %% Color List based on origin color4line
            % OUTPUT
            %   c: colour codes for plots
            
            c = [[81 81 81]; %black
                [241 64 64]; %red
                [26 111 223]; % blue
                [55 173 107] %green: do not use with red
                [177 119 222] %purple
                [204 153 0] %orca
                [0 203 204] %turqouise
                [125 78 78] %dark brown
                [142 142 0] %after this hard to see for colourblind people
                [251 101 1]
                [102 153 204]
                [111 184 2]]; 
            
            c = 1/255.* c;

        end
        
        function boundaryAnalysisRxFractions(obj, phase)
            %% Analysis of boundaries of different classified Rx grains with deformed grains
            %% (applied after function classifyDRX has been used to classify grains)
            % INPUT
            %   phase: name of phase of interest
            %
            % OUTPUT
            %   figures
            %
            % LEGEND (Variable names)
            %   Rx = RxGrains
            %   Def = DeformedGrains
            %   Psn = PSNGrains; %PSN Grains
            %   Drx = DRXGrains; %DRX grains (Rx without PSN)
            %   DrxIg = DRXGrainsWithinGrain;
            %   DrxGb = DRXGrainsOnGBs;
            %   DrxTb = DRXGrainsOnPriorTwinGBs;        

            %% Rx GBs with Deformed Grains(Nucleation Mechanism: Twin vs SIBM)
            %DDRX boundaries with deformed Grains:
            gbDef_DrxGb = obj.DRXGrainsOnGBs.boundary(obj.DeformedGrains);
            gbDef_DrxGb = gbDef_DrxGb(phase, phase);
            %Of these boundaries, these are twinned:
            gbDef_DrxGbCSL3 = gbDef_DrxGb(gbDef_DrxGb.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
            gbDef_DrxGbCSL9 = gbDef_DrxGb(gbDef_DrxGb.isTwinning(CSL(9,obj.Ebsd(phase).CS), (15/sqrt(9))*degree));
            gbDef_DrxGbCSL27 = gbDef_DrxGb(gbDef_DrxGb.isTwinning(CSL(27,obj.Ebsd(phase).CS), (15/sqrt(27))*degree));
            %these are not twinned:
            gbDef_DrxGbNT = gbDef_DrxGb(~gbDef_DrxGb.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
                       
            %%PSN boundaries with deformed Grains: 
            %gbDef_Psn = obj.PSNGrains.boundary(obj.DeformedGrains);
            
            %DRX Grains within grains & deformed grains:
            gbDef_DrxIg = obj.DRXGrainsWithinGrain.boundary(obj.DeformedGrains);
            if (~isempty(gbDef_DrxIg))
                gbDef_DrxIgbCSL3 = gbDef_DrxIg(gbDef_DrxIg.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
                gbDef_DrxIgbCSL9 = gbDef_DrxIg(gbDef_DrxIg.isTwinning(CSL(9,obj.Ebsd(phase).CS), (15/sqrt(9))*degree));
                gbDef_DrxIgbCSL27 = gbDef_DrxIg(gbDef_DrxIg.isTwinning(CSL(27,obj.Ebsd(phase).CS), (15/sqrt(27))*degree));

                gbDef_DrxIgNT = gbDef_DrxIg(~gbDef_DrxIg.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
            else
                gbDef_DrxIgbCSL3 = [];
                gbDef_DrxIgbCSL3.length = 0;
                gbDef_DrxIgbCSL3.segLength = 0;
                
                gbDef_DrxIgbCSL9 = gbDef_DrxIgbCSL3;
                gbDef_DrxIgbCSL27 = gbDef_DrxIgbCSL3;

                gbDef_DrxIgNT = [];
                gbDef_DrxIgNT.length = 0;
                gbDef_DrxIgNT.segLength = 0;
            end

            %DRX Grains on prior Twin boundaries & deformed grains:
            gbDef_DrxTb = obj.DRXGrainsOnPriorTwinGBs.boundary(obj.DeformedGrains);
            if (~isempty(gbDef_DrxTb))
                gbDef_DrxTbCSL3 = gbDef_DrxTb(gbDef_DrxTb.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
                gbDef_DrxTbCSL9 = gbDef_DrxTb(gbDef_DrxTb.isTwinning(CSL(9,obj.Ebsd(phase).CS), (15/sqrt(9))*degree));
                gbDef_DrxTbCSL27 = gbDef_DrxTb(gbDef_DrxTb.isTwinning(CSL(27,obj.Ebsd(phase).CS), (15/sqrt(27))*degree));

                gbDef_DrxTbNT = gbDef_DrxTb(~gbDef_DrxTb.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
            else
                gbDef_DrxTbCSL3 = [];
                gbDef_DrxTbCSL3.length = 0;
                gbDef_DrxTbCSL3.segLength = 0;
                
                gbDef_DrxTbCSL9 = gbDef_DrxTbCSL3;
                gbDef_DrxTbCSL27 = gbDef_DrxTbCSL3;
                
                gbDef_DrxTbNT = [];
                gbDef_DrxTbNT.length = 0;
                gbDef_DrxTbNT.segLength = 0;
            end
            
            %PSN Grains
            if(~isempty(obj.PSNGrains))
                gbDef_DrxPsn = obj.PSNGrains.boundary(obj.DeformedGrains);
                gbDef_DrxPsnCSL3 = gbDef_DrxPsn(gbDef_DrxPsn.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
                gbDef_DrxPsnCSL9 = gbDef_DrxPsn(gbDef_DrxPsn.isTwinning(CSL(9,obj.Ebsd(phase).CS), (15/sqrt(9))*degree));
                gbDef_DrxPsnCSL27 = gbDef_DrxPsn(gbDef_DrxPsn.isTwinning(CSL(27,obj.Ebsd(phase).CS), (15/sqrt(27))*degree));
            else
                gbDef_DrxPsn.length = 0;
                gbDef_DrxPsnCSL3 = [];
                gbDef_DrxPsnCSL3.length = 0;
                gbDef_DrxPsnCSL3.segLength = 0;
                
                gbDef_DrxPsnCSL9 = gbDef_DrxPsnCSL3;
                gbDef_DrxPsnCSL27 = gbDef_DrxPsnCSL3;

            end
            
            %% Boundaries within Drx (GB) Fraction - no twins as excluded 
            %GBs just between recrystallised DDRX grains
            gbDrxGb_DrxGb = obj.DRXGrainsOnGBs.boundary(phase,phase);
            gbDrxGb_DrxGb = gbDrxGb_DrxGb(~gbDrxGb_DrxGb.hasGrain(obj.DeformedGrains));
            
            %% Plotting Phase Map           
            colorList = 1/255.* [[0 158 115]; [240 228 66]; [230 159 0]; [213 94 0]];

            gbDef_Rx = obj.RxGrains.boundary(obj.DeformedGrains);
            gB3 = gbDef_Rx(gbDef_Rx.isTwinning(CSL(3,obj.Ebsd(phase).CS),(15/sqrt(3))*degree));
            gB9 = gbDef_Rx(gbDef_Rx.isTwinning(CSL(9,obj.Ebsd(phase).CS), (15/sqrt(9))*degree));
            gB27 = gbDef_Rx(gbDef_Rx.isTwinning(CSL(27,obj.Ebsd(phase).CS), (15/sqrt(27))*degree));
            
            % Plot CSL3 boundaries between different Rx fractions and
            % deformed grains
            f = plot(obj, 'phase', 2,0,0)
            hold on
            try
                plot(gbDef_DrxGbCSL3, 'linecolor', 1/255*[0 158 115], 'lineWidth', 3, 'DisplayName', 'CSL3 btw. DRX (on GBs)/Deformed')
            end
            
            try
                plot(gbDef_DrxTbCSL3, 'linecolor', 1/255*[230 159 0], 'lineWidth', 3, 'DisplayName', 'CSL3 btw. DRX (on Twins)/Deformed')
            end
            
            try
                plot(gbDef_DrxIgbCSL3, 'linecolor', 1/255*[240 228 66], 'lineWidth', 3, 'DisplayName', 'CSL3 btw. DRX (Within Grains)/Deformed')
            end
            
            try
                plot(gbDef_DrxPsnCSL3, 'linecolor', 1/255*[213 94 0], 'lineWidth', 3, 'DisplayName', 'CSL3 btw. DRX (PSN)/Deformed')
            end
                hold off
            f.Name = ['RxBoundaryChar_' phase '_Twins-RXvsDef'];
            
            %% Plot Histogram with Misorientation of RX grains to deformed
            
            figure('Name', ['RxBoundaryChar_' phase]);
            edges = 0:5:65;
            h1 = histcounts(gbDef_DrxGb(phase,phase).misorientation.angle/degree, edges, 'Normalization', 'probability');
            if (gbDef_DrxIg.length > 0)
                h2 = histcounts(gbDef_DrxIg(phase,phase).misorientation.angle/degree, edges, 'Normalization', 'probability');
            else
                h2 = histcounts([],edges);
            end
            
            if (gbDef_DrxTb.length > 0)
                h3 = histcounts(gbDef_DrxTb(phase,phase).misorientation.angle/degree, edges, 'Normalization', 'probability');
            else
                h3 = histcounts([],edges);
            end
            
            if (gbDef_DrxPsn.length > 0)
                h4 = histcounts(gbDef_DrxPsn(phase,phase).misorientation.angle/degree, edges, 'Normalization', 'probability');
            else
                h4 = histcounts([],edges);
            end
            b = bar(edges(1:end-1),[h1; h2; h3; h4]',1);
            
            hold on
            [x1,p1] = obj.calcRandomMisorientation(edges(1:end-1), h1); %both histograms give the same output
            plot(x1,p1, '--', 'linewidth',1, 'color', 'black', 'DisplayName', 'Random distribution')
            legend([{'Deformed/DRX Grains (on GB)'}; {'Deformed/DRX Grains (not on GBs)'};{'Deformed/DRX Grains (on prior Twins)'};{'Deformed/DRX Grains (PSN)'};])
            xlabel(['Misorientation [°]'])
            ylabel('Relative Frequency [-]')
            legend show
            xlim([edges(1) edges(end)])
            b(1).FaceColor = colorList(1,:);
            b(2).FaceColor = colorList(2,:);
            b(3).FaceColor = colorList(3,:);
            b(4).FaceColor = colorList(4,:);
            
            %% IPF Maps of GBs
            %DDRX vs DDRX(without Twins): shift to [011] as in twins in
            %deformed sample -> indiciation that there was more twinning as
            %nucleation mechanism but grains rotated away from it
            
            if (gbDef_DrxTb.length > 0)
                f = figure;
                plot (calcDensity( gbDef_DrxTb(phase,phase).misorientation.axis,10*degree), 'contourf')
                %mtexColorMap LaboTeX
                mtexColorbar
                %setColorRange([0,1.8])
                f.Name = ['RxBoundaryIPF_' phase '_OnTwins-DRXvsDef'];
            end
            
            f = figure;
            plot (calcDensity( gbDef_DrxGb(phase,phase).misorientation.axis,10*degree), 'contourf')
            %mtexColorMap LaboTeX
            mtexColorbar
            %setColorRange([0,1.8])
            f.Name = ['RxBoundaryIPF_' phase '_Def-GBDRX'];
            
            f = figure;
            plot (calcDensity( gbDef_DrxGbNT(phase,phase).misorientation.axis,10*degree), 'contourf')
            %mtexColorMap LaboTeX
            mtexColorbar
            %setColorRange([0,1.8])
            f.Name = ['RxBoundaryIPF_' phase '_Def-GBDRXwithoutTwins'];
            
            f = figure
            plot (calcDensity(gbDrxGb_DrxGb.misorientation.axis,'halfwidth',10*degree), 'contourf')
            %mtexColorMap LaboTeX
            mtexColorbar
            %setColorRange([0,1.8])
            f.Name = ['RxBoundaryIPF_' phase '_GBDRX-GBDRXwithoutTwins'];

            %Get just DRX grains which border other DRX grains
            % Misorientation DRX/DRX and DRX/Deformed Grains (without twins):
            figure('Name', ['RxBoundaryChar_' phase '_DRXonGBs_withoutTwins']);
            h1 = histogram(gbDef_DrxGbNT(phase,phase).misorientation.angle/degree,0:1:65,'Normalization', 'probability', 'DisplayName', 'Deformed/DRX Grains (on GB)')
            hold on;
            histogram(gbDrxGb_DrxGb(phase,phase).misorientation.angle/degree,0:1:65,'Normalization', 'probability', 'DisplayName', 'DRXGrains/DRX Grains (on GB)')
            xlabel(['Misorientation [°]'])
            ylabel('Relative Frequency')
            legend show
            [x1,p1] = obj.calcRandomMisorientation(h1.BinEdges(1:end-1), h1.Values); %both histograms give the same output
             hold on
             plot(x1,p1, '--', 'linewidth',1, 'color', 'black', 'DisplayName', 'Random distribution')
            
            %% Write Report
            %Fraction of Twinning
            obj.writeReport('-', '-');
            obj.writeReport('Analysis Rx Boundary', phase);
            obj.writeReport('Deformed vs. DRX-GB GB Fraction [num.]', gbDef_DrxGb.length/obj.DRXGrainsOnGBs.boundary.length)
            obj.writeReport('Deformed vs. DRX-GB GB Fraction [length.]', sum(gbDef_DrxGb.segLength)/sum(obj.DRXGrainsOnGBs.boundary.segLength))
            obj.writeReport('Deformed vs. DRX-GB CSL3 [num.%]', gbDef_DrxGbCSL3.length/gbDef_DrxGb.length*100); %num
            obj.writeReport('Deformed vs. DRX-GB CSL3 [length.%]', sum(gbDef_DrxGbCSL3.segLength)/sum(gbDef_DrxGb.segLength)*100); %length
            obj.writeReport('Deformed vs. DRX-GB CSL9 [num.%]', gbDef_DrxGbCSL9.length/gbDef_DrxGb.length*100); %num
            obj.writeReport('Deformed vs. DRX-GB CSL9 [length.%]', sum(gbDef_DrxGbCSL9.segLength)/sum(gbDef_DrxGb.segLength)*100); %length
            obj.writeReport('Deformed vs. DRX-GB CSL27 [num.%]', gbDef_DrxGbCSL27.length/gbDef_DrxGb.length*100); %num
            obj.writeReport('Deformed vs. DRX-GB CSL27 [length.%]', sum(gbDef_DrxGbCSL27.segLength)/sum(gbDef_DrxGb.segLength)*100); %length
            
            %Within Grains
            obj.writeReport('Deformed vs. DRX-withinGrains GB Fraction [num.]', gbDef_DrxIg.length/obj.DRXGrainsWithinGrain.boundary.length)
            obj.writeReport('Deformed vs. DRX-withinGrains GB Fraction [length.]', sum(gbDef_DrxIg.segLength)/sum(obj.DRXGrainsWithinGrain.boundary.segLength))
            obj.writeReport('Deformed vs. DRX-withinGrains CSL3 [num.%]', gbDef_DrxIgbCSL3.length/gbDef_DrxIg.length*100); %some of TB still included
            obj.writeReport('Deformed vs. DRX-withinGrains CSL3 [length.%]', sum(gbDef_DrxIgbCSL3.segLength)/sum(gbDef_DrxIg.segLength)*100); %length  
            obj.writeReport('Deformed vs. DRX-withinGrains CSL9 [num.%]', gbDef_DrxIgbCSL9.length/gbDef_DrxIg.length*100); %some of TB still included
            obj.writeReport('Deformed vs. DRX-withinGrains CSL9 [length.%]', sum(gbDef_DrxIgbCSL9.segLength)/sum(gbDef_DrxIg.segLength)*100); %length  
            obj.writeReport('Deformed vs. DRX-withinGrains CSL27 [num.%]', gbDef_DrxIgbCSL27.length/gbDef_DrxIg.length*100); %some of TB still included
            obj.writeReport('Deformed vs. DRX-withinGrains CSL27 [length.%]', sum(gbDef_DrxIgbCSL27.segLength)/sum(gbDef_DrxIg.segLength)*100); %length  
                            
             %Fraction of twins 
            obj.writeReport('Deformed vs. DRX-priorTwins GB Fraction [num.]', gbDef_DrxTb.length/obj.DRXGrainsOnPriorTwinGBs.boundary.length)
            obj.writeReport('Deformed vs. DRX-priorTwins GB Fraction [length.]', sum(gbDef_DrxTb.segLength)/sum(obj.DRXGrainsOnPriorTwinGBs.boundary.segLength))  
            obj.writeReport('Deformed vs. DRX-priorTwins CSL3 [num.%]', gbDef_DrxTbCSL3.length/gbDef_DrxTb.length*100); %looks like easier to twin for them!
            obj.writeReport('Deformed vs. DRX-priorTwins CSL3 [length.%]', sum(gbDef_DrxTbCSL3.segLength)/sum(gbDef_DrxTb.segLength)*100); %length
            obj.writeReport('Deformed vs. DRX-priorTwins CSL9 [num.%]', gbDef_DrxTbCSL9.length/gbDef_DrxTb.length*100); %looks like easier to twin for them!
            obj.writeReport('Deformed vs. DRX-priorTwins CSL9 [length.%]', sum(gbDef_DrxTbCSL9.segLength)/sum(gbDef_DrxTb.segLength)*100); %length
            obj.writeReport('Deformed vs. DRX-priorTwins CSL27 [num.%]', gbDef_DrxTbCSL27.length/gbDef_DrxTb.length*100); %looks like easier to twin for them!
            obj.writeReport('Deformed vs. DRX-priorTwins CSL27 [length.%]', sum(gbDef_DrxTbCSL27.segLength)/sum(gbDef_DrxTb.segLength)*100); %length

            %Fraction of twins PSN
            if(~isempty(obj.PSNGrains))
                obj.writeReport('Deformed vs. DRX-PSN GB Fraction [num.]', gbDef_DrxPsn.length/obj.PSNGrains.boundary.length)
                obj.writeReport('Deformed vs. DRX-PSN GB Fraction [length.]', sum(gbDef_DrxPsn.segLength)/sum(obj.PSNGrains.boundary.segLength))  
                obj.writeReport('Deformed vs. DRX-PSN CSL3 [num.%]', gbDef_DrxPsnCSL3.length/gbDef_DrxPsn.length*100); %looks like easier to twin for them!
                obj.writeReport('Deformed vs. DRX-PSN CSL3 [length.%]', sum(gbDef_DrxPsnCSL3.segLength)/sum(gbDef_DrxPsn.segLength)*100); %length
                obj.writeReport('Deformed vs. DRX-PSN CSL9 [num.%]', gbDef_DrxPsnCSL9.length/gbDef_DrxPsn.length*100); %looks like easier to twin for them!
                obj.writeReport('Deformed vs. DRX-PSN CSL9 [length.%]', sum(gbDef_DrxPsnCSL9.segLength)/sum(gbDef_DrxPsn.segLength)*100); %length
                obj.writeReport('Deformed vs. DRX-PSN CSL27 [num.%]', gbDef_DrxPsnCSL27.length/gbDef_DrxPsn.length*100); %looks like easier to twin for them!
                obj.writeReport('Deformed vs. DRX-PSN CSL27 [length.%]', sum(gbDef_DrxPsnCSL27.segLength)/sum(gbDef_DrxPsn.segLength)*100); %length
            else
                obj.writeReport('Deformed vs. DRX-PSN GB Fraction', 'no PSN grains')
            end
        end
        
        function [results, misorientation, gos] = boundaryAnalysisRxOnTwins(obj, phase)
           %% Information on the misorientation of prior-twins, which nucleated Rx grains
           % just information on grains which lay between exactly 2 deforme grains
           %
           % OUTPUT
           %    results: 
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
           
           gDef = obj.Grains(obj.Grains.prop.GOS./degree >= obj.RxThreshold); 
           gDef = gDef(phase);
           
           originalIds = obj.DRXGrainsOnPriorTwinGBs.prop.originalId;
           originalIds = originalIds(~isnan(originalIds));
           grains = obj.Grains('id', originalIds);  %Problem: surrounding deformed grains already merged
           ebsd = obj.Ebsd;
                
%%         Selecting Rx grains on prior-twins which are wedged in between two grains
           misorientation = [];
           gos = [];
           n = []; %get n to be class grain2d
           gbDrxDef1 = [];
           gbDrxDef2 = [];
           moriDef1Def2 = [];
           moriDef1Def2_mean = [];
           grainsDef1 = [];
           grainsDef2 = [];
           grainsDRX = [];
           
           % could be solved via matrix calculation instead of loops but
           % not enough memory!
           for i=1:grains.length
               g = grains(i);
               gB = g.boundary(gDef);  
               
               gId = unique(gB.grainId);
               gId(gId == g.id) = [];
               
               %% Check that just 2 deformed grains, which are touching each other
               if(length(gId) == 2 && ~isempty(gDef('id', gId(1)).boundary(gDef('id', gId(2)))))

                   %For coherency analysis
                   if (isempty(grainsDRX))
                       grainsDRX = grains(i);
                   else
                       grainsDRX = [grainsDRX grains(i)];
                   end  
                   
                   ori = ebsd('id', gB.ebsdId).orientations;

                   aGb = gB(gB.grainId==gId(1));
                   bGb = gB(gB.grainId==gId(2));
                   
                   gB.prop.ori1 = ori(:,1);
                   gB.prop.ori2 = ori(:,2);

                   %Misorientation
                   a1 = mean(gB(gB.grainId==gId(1)).prop.ori1); %Deformed grain 1
                   a2 = mean(gB(gB.grainId==gId(1)).prop.ori2); %Rx grain 1
                   b1 = mean(gB(gB.grainId==gId(2)).prop.ori1); %Deformed grain 2
                   b2 = mean(gB(gB.grainId==gId(2)).prop.ori2); %Rx grain 1 (other side)

                   %Misorientation
                   a2b2 = angle(a2,b2)/degree; % angle between both deformed grains
                   a1a2 = angle(a1,a2)/degree;% angle between Rx & deformed 1
                   b1b2 = angle(b1,b2)/degree;% angle between Rx & deformed 2

                   %<111> Axis deviation
                   mori_a2b2_d = inv(a2) .* b2;
                   a2b2_d = min(angle(mori_a2b2_d.symmetrise * Miller(1,1,1, ebsd(phase).CS), Miller(1,1,1, ebsd(phase).CS)))./degree;       
                   mori1 = inv(a1) .* a2;
                   a1a2_d = min(angle(mori1.symmetrise * Miller(1,1,1, ebsd(phase).CS), Miller(1,1,1, ebsd(phase).CS)))./degree;
                   mori2 = inv(b1) .* b2;
                   b1b2_d = min(angle(mori2.symmetrise * Miller(1,1,1, ebsd(phase).CS), Miller(1,1,1, ebsd(phase).CS)))./degree;

                   gos1 = gDef(gDef.id == gId(1)).GOS/degree;
                   gos2 = gDef(gDef.id == gId(2)).GOS/degree;
                   gos3 = g.GOS/degree;
                   
                   if (a1a2 > b1b2) %program better -> no double code!
                       misorientation = [misorientation; g.id a2b2 a2b2_d a1a2 a1a2_d b1b2 b1b2_d];
                       gos =  [gos; gos1 gos2 gos3];
                       
                       if (isempty(gbDrxDef1))
                           gbDrxDef1 = aGb;
                           gbDrxDef2 = bGb;
                       else
                           gbDrxDef1 = [gbDrxDef1 aGb];
                           gbDrxDef2 = [gbDrxDef2 bGb];
                       end
                       
                       % For coherency check!
                       if (isempty(grainsDef1))
                            grainsDef1 =  gDef('id', gId(1));
                            grainsDef2 =  gDef('id', gId(2));
                       else 
                           grainsDef1 = [grainsDef1 gDef('id', gId(1))];
                           grainsDef2 = [grainsDef2 gDef('id', gId(2))];
                       end
                   else
                       misorientation = [misorientation; g.id a2b2 a2b2_d b1b2 b1b2_d a1a2 a1a2_d];
                       gos = [gos; gos2 gos1 gos3];
                       
                       if (isempty(gbDrxDef1))
                           gbDrxDef1 = bGb;
                           gbDrxDef2 = aGb;
                       else
                           gbDrxDef1 = [gbDrxDef1 bGb];
                           gbDrxDef2 = [gbDrxDef2 aGb];
                       end  
                       
                       % For coherence check!
                       if (isempty(grainsDef1))
                            grainsDef1 =  gDef('id', gId(2));
                            grainsDef2 =  gDef('id', gId(1));
                       else 
                           grainsDef1 = [grainsDef1 gDef('id', gId(2))];
                           grainsDef2 = [grainsDef2 gDef('id', gId(1))];
                       end
                   end
                   
                   %% Mean Orientation of grains
                   meanDef1 = gDef('id', gId(1)).meanOrientation;
                   meanDef2 = gDef('id', gId(2)).meanOrientation;
                   
                   mori_mean = inv(meanDef1) .* meanDef2;
                   mean_a = mori_mean.angle./degree;
                   mean_d = min(angle(mori_mean.symmetrise * Miller(1,1,1, ebsd(phase).CS), Miller(1,1,1, ebsd(phase).CS)))./degree;       
   
                   if (isempty(moriDef1Def2))
                       moriDef1Def2 = mori_a2b2_d;
                       moriDef1Def2_mean = [mean_a mean_d];
                   else
                       moriDef1Def2 = [moriDef1Def2 mori_a2b2_d];
                       moriDef1Def2_mean = [moriDef1Def2_mean; mean_a mean_d];
                   end
               else
                   %misorientation = [misorientation; g.id 0 0 0 0 0 0]
                   if (isempty(n))
                       n = g;
                   else
                       n = [n g];
                   end
               end
           end
           
           %% Treatment as segments:
           %Misorientation above takes averages of each GB (helps to find relation between LAGB/HAGB/deformed -> none found in 1:1 comparison)
           %to treat every GB segment seperatelty (just possible for DRX/deformed):
           
           % High angle GBs
           if (isempty(gbDrxDef1))
               results = [];
               obj.writeReport('-','-');
               obj.writeReport('DRX on distorted twins',phase);
               obj.writeReport('N/A','No DRX grains on distorted twins between two deformed grains');
               return
           end
           gbDrxDef1_misorientation(:,1) = gbDrxDef1.misorientation.angle./degree;
           
           ori = ebsd('id', gbDrxDef1.ebsdId).orientations;
           gbDrxDef1.prop.ori1 = ori(:,1);
           gbDrxDef1.prop.ori2 = ori(:,2);
           mori = inv(gbDrxDef1.prop.ori1) .* gbDrxDef1.prop.ori2;
           gbDrxDef1_misorientation(:,2) = min(angle(mori.symmetrise .* Miller(1,1,1, ebsd(phase).CS), Miller(1,1,1, ebsd(phase).CS)))./degree;       
           
           % Low angle GBs
           gbDrxDef2_misorientation(:,1) = gbDrxDef2.misorientation.angle./degree;
           
           ori = ebsd('id', gbDrxDef2.ebsdId).orientations;
           gbDrxDef2.prop.ori1 = ori(:,1);
           gbDrxDef2.prop.ori2 = ori(:,2);
           mori = inv(gbDrxDef2.prop.ori1) .* gbDrxDef2.prop.ori2;
           gbDrxDef2_misorientation(:,2) = min(angle(mori.symmetrise .* Miller(1,1,1, ebsd(phase).CS), Miller(1,1,1, ebsd(phase).CS)))./degree; 
           
           %% Output        
           colorList = obj.getColorList;
           %Histogram: Misorientation of different GBs          
           edges = 0:5:65;
           h1 = histcounts(misorientation(:,2), edges, 'Normalization', 'probability');
           h2 = histcounts(gbDrxDef1_misorientation(:,1), edges, 'Normalization', 'probability');
           h3 = histcounts(gbDrxDef2_misorientation(:,1), edges, 'Normalization', 'probability');
           figure('Name', ['DrxOnTwins__' phase '_Misorientation']);
           b1 = bar(edges(1:end-1),[h1; h2; h3]', 1);
           [x1,p1] = obj.calcRandomMisorientation(edges(1:end-1), h1); %sum(h1)=sum(h2)=sum(h3)
           hold on
           plot(x1,p1, '--', 'linewidth',1, 'color', 'black')
           legend([{'Distorted twin 1/distorted twin 2'}; {'DRX grain/distorted twin 1'}; {'DRX grain /distorted twin 2'}; {'Random distribution'}])
           xlabel(['Misorientation [°]'])
           ylabel('Relative Frequency [-]')
           legend show
           xlim([-1 edges(end)])
           b1(1).FaceColor = colorList(2,:);
           b1(2).FaceColor = colorList(3,:);
           b1(3).FaceColor = colorList(6,:);
           
           %Histogram: Deviation from <111>
           edges = 0:2:30;
           h1 = histcounts(misorientation(:,3), edges, 'Normalization', 'probability');
           %h2 = histcounts(misorientation(:,5), edges, 'Normalization', 'probability');
           h2 = histcounts(gbDrxDef1_misorientation(:,2), edges, 'Normalization', 'probability');
           %h3 = histcounts(misorientation(:,7), edges, 'Normalization', 'probability');
           h3 = histcounts(gbDrxDef2_misorientation(:,2), edges, 'Normalization', 'probability');
           figure('Name', ['DrxOnTwins__' phase '_DeviationFrom111']);
           b2 = bar(edges(1:end-1),[h1; h2; h3]', 1);
           [x2,p2] = obj.calcRandomAxisDeviation(111, edges(1:end-1), h1); %sum(h1)=sum(h2)=sum(h3)
           hold on
           plot(x2,p2, '--', 'linewidth',1, 'color', 'black')
           legend([{'Distorted twin 1/distorted twin 2'}; {'DRX grain/distorted twin 1'}; {'DRX grain /distorted twin 2'}; {'Random distribution'}])
           xlabel(['Deviation from <111> rotation axis [°]'])
           ylabel('Relative Frequency [-]')
           legend show
           xlim([-1 edges(end)])
           b2(1).FaceColor = colorList(2,:);
           b2(2).FaceColor = colorList(3,:);
           b2(3).FaceColor = colorList(6,:);

%            %Histogram: Deviation from <111>
%            edges = 0:2:30;
%            h1 = histcounts(misorientation(:,3), edges, 'Normalization', 'probability');
%            %h2 = histcounts(misorientation(:,5), edges, 'Normalization', 'probability');
%            h2 = histcounts(moriDef1Def2_mean(:,2), edges, 'Normalization', 'probability');
%            figure('Name', ['DrxOnTwins__' phase '_MeanPriorTwin111']);
%            b3 = bar(edges(1:end-1),[h1; h2]',1);
%            [x2,p2] = obj.calcRandomAxisDeviation(111, edges(1:end-1), h1); %sum(h1)=sum(h2)=sum(h3)
%            hold on
%            plot(x2,p2, '--', 'linewidth',1, 'color', 'black')
%            legend([{'Distorted twin (GB average at DRX position)'}; {'Distorted twin (mean grain average)'}; {'Random distribution'}])
%            xlabel(['Deviation from <111> rotation axis [°]'])
%            ylabel('Relative Frequency [-]')
%            legend show
%            xlim([-1 edges(end)])
%            b3(1).FaceColor = colorList(2,:);
%            b3(2).FaceColor = colorList(1,:);

                
           %% Percentage of CSL3 GBs:
           results = [];
           
           %Results structure:
           %Rows: Twin-Twin | GB1 | GB2
           %Columns: CSL3 fraction | Misori mean | Misori std | Misori std
           %| <111> deviation mean | <111> dev std 
           %| point misori between twin/GB mean | std
           %| percentage of CSL3 which are coherent
           
           limit = 8.66; %brandon for CSL3
           results(1,1) = sum((abs(60-misorientation(:,2)) < limit) .* (misorientation(:,3) < limit))/length(misorientation)*100;
           results(2,1) = sum(gbDrxDef1(angle(gbDrxDef1.misorientation,CSL(3,obj.Ebsd(phase).CS)) < 8.66*degree).segLength)/sum(gbDrxDef1.segLength)*100;
           results(3,1) = sum(gbDrxDef2(angle(gbDrxDef2.misorientation,CSL(3,obj.Ebsd(phase).CS)) < 8.66*degree).segLength)/sum(gbDrxDef2.segLength)*100;

           obj.writeReport('-','-');
           obj.writeReport('DRX on distorted twins',phase);
           obj.writeReport('CSL3 Fraction [average GB %]: distorted twin 1/distorted twin 2',results(1,1));
           obj.writeReport('CSL3 Fraction [length. %}: DRX grain/distorted twin 1',results(2,1));
           obj.writeReport('CSL3 Fraction [length. %}: DRX grain/distorted twin 2',results(3,1));
           
           %% Mean values
           % Betweewen prior Twins
           results(1,2) = mean(misorientation(:,2)); %Misorientation mean
           results(1,3) = std(misorientation(:,2)); %Misori Std
           results(1,4) = mean(misorientation(:,3)); %111 deviation mean
           results(1,5) = std(misorientation(:,3)); %111 deviation std
           
           %Between prior-twin & Rx-GB 1
           results(2,2) = mean(misorientation(:,4));
           results(2,3) = std(misorientation(:,4));
           results(2,4) = mean(misorientation(:,5));
           results(2,5) = std(misorientation(:,5));
           
           %Between prior-twin & Rx-GB 2
           results(3,2) = mean(misorientation(:,6));
           results(3,3) = std(misorientation(:,6));
           results(3,4) = mean(misorientation(:,7));
           results(3,5) = std(misorientation(:,7));
           
           %% Deviation of each point
           %Between Twin & Rx-GB 1
           results(1,6) = 0;
           results(1,7) = 0;
           
           gb1 = abs(misorientation(:,2) - misorientation(:,4));
           results(2,6) = mean(gb1);
           results(2,7) = std(gb1);
           
           %Between Twin & Rx-GB 2
           gb2 = abs(misorientation(:,2) - misorientation(:,6));
           results(3,6) = mean(gb2);
           results(3,7) = std(gb2); 
           
%            %% Plot misorientation axis
%            f1 = figure('Name', ['DrxOnTwins__' phase '_MisorientationAxis']);
%            hold on
%            plot(gbDrxDef1.misorientation.axis, 'fundamentalRegion','MarkerColor', colorList(3,:),'MarkerSize',2, 'DisplayName', 'DRX grain/distorted twin 1')
%            hold on
%            plot(gbDrxDef2.misorientation.axis, 'fundamentalRegion','MarkerColor', colorList(6,:),'MarkerSize',2, 'DisplayName', 'DRX grain/distorted twin 2')
%            legend show
%            f1.Name = ['TwinRot_PriorTwinsVsCSL3_MisorientationAxis_' phase];
           
           %% Check angle to lattice planes
           %Smooth GBs
           smoothedGrains = obj.Grains.smooth(20);
           grainsDRXSm = smoothedGrains('id', grainsDRX.id);
           grainsDef1Sm = smoothedGrains('id', grainsDef1.id);
           grainsDef2Sm = smoothedGrains('id', grainsDef2.id);
           
           %clear def1_gb def1_minAngle1 def1_minAngle2 def2_gb def2_minAngle1 def2_minAngle2
           %identity = 'allMinusCoherent';
           identity = 'incoherent'; %just GB planes of incoherent CSL3 are analysed
           
           clear def1_gb def2_gb def1_minAngle1 def2_minAngle1 def1_minAngle2 def2_minAngle2 def1_gb_AllGBs def2_gb_AllGBs
           
           w0 = waitbar(0, 'Grains Analysed');
           for i=1:length(grainsDRXSm)
               w0 = waitbar(i/length(grainsDRXSm));
               [def1_gb{i}, def1_minAngle1{i}, def1_minAngle2{i}, m] = obj.gbdirectionPlanetraceAngle(phase, grainsDRXSm(i), grainsDef1Sm(i), identity, false);
               [def2_gb{i}, def2_minAngle1{i}, def2_minAngle2{i}, m] = obj.gbdirectionPlanetraceAngle(phase, grainsDRXSm(i), grainsDef2Sm(i), identity, false);
                def1_gb_AllGBs{i} = grainsDRXSm(i).boundary(grainsDef1Sm(i));
                def2_gb_AllGBs{i} = grainsDRXSm(i).boundary(grainsDef2Sm(i));
           end
           close(w0);
           %Combine cells into long file:
           def1_gb = cat(1, def1_gb{:});
           def1_minAngle1 = cat(1, def1_minAngle1{:});
           def1_minAngle2 = cat(1, def1_minAngle2{:});
           def1_gb_AllGBs = cat(1, def1_gb_AllGBs{:});
           def2_gb = cat(1, def2_gb{:});
           def2_minAngle1 = cat(1, def2_minAngle1{:});
           def2_minAngle2 = cat(1, def2_minAngle2{:});
           def2_gb_AllGBs = cat(1, def2_gb_AllGBs{:});
           
%            obj.plotGbdirectionPlanetraceAngle(phase, grainsDRXSm, grainsDef1Sm, def1_gb, def1_minAngle1, def1_minAngle2, m, identity)
%            obj.plotGbdirectionPlanetraceAngle(phase, grainsDRXSm, grainsDef2Sm, def2_gb, def2_minAngle1, def2_minAngle2, m, identity)
           
           %% Deformed vs deformed
           %identity = 'allMinusCoherent';
           
           idA = [grainsDef1Sm.id grainsDef2Sm.id];
           [idB, iA, iB] = unique(idA, 'rows');
           
           idB_1 = [];
           idB_2 = [];
           
           for i=1:length(idB(:,1))
               
               %check if all Twin-DRX grains have high misorientation to
               %the same deformed grain
               k = idB(i,1) == idB(:,2);
               if (any(idB(k,1) == idB(i,2)))
                   %If different Rx-grains have misorientation to different
                   %grains
                   idB_1 = [idB_1; idB(i,1) idB(i,2)];
               else
                   idB_2 = [idB_2; idB(i,1) idB(i,2)];
               end
           end

           idB_1 = sort(idB_1.').';
           idB_1  = unique(idB_1, 'rows');
           
           obj.writeReport('Fraction of where DRX grains have high/low angle GB aligned',length(idB_2)/(length(idB_1)+length(idB_2))) %percentage where all misoriented DRX have high misorientation GB to same prior-twin

           clear defDef_gb_1 defDef_minAngle1_1 defDef_minAngle2_1 defDef_gb_1_AllGBs defDef_gb_2_AllGBs
           
           if (~isempty(idB_1))
               w1 = waitbar(0, 'Grains Analysed');
               for i=1:length(idB_1(:,1))
                   w1 = waitbar(i/length(idB_1(:,1)), w1);
                   [defDef_gb_1{i}, defDef_minAngle1_1{i}, defDef_minAngle2_1{i}, m] = obj.gbdirectionPlanetraceAngle(phase, smoothedGrains('id', idB_1(i,1)), smoothedGrains('id', idB_1(i,2)), identity, false);
                    defDef_gb_1_AllGBs{i} = smoothedGrains('id', idB_1(i,1)).boundary(smoothedGrains('id', idB_1(i,2)));
               end
               close(w1);
           
               defDef_gb_1 = cat(1, defDef_gb_1{:});
               defDef_minAngle1_1 = cat(1, defDef_minAngle1_1{:});
               defDef_minAngle2_1 = cat(1, defDef_minAngle2_1{:});
               defDef_gb_1_AllGBs = cat(1, defDef_gb_1_AllGBs{:});

%                obj.plotGbdirectionPlanetraceAngle(phase, smoothedGrains('id', idB_1(:,1)), smoothedGrains('id', idB_1(:,2)), defDef_gb_1, defDef_minAngle1_1, defDef_minAngle2_1, m, identity)
           else
               defDef_gb_1 = [];
               defDef_minAngle1_1 = [];
               defDef_minAngle2_1 = [];
               defDef_gb_1_AllGBs = [];
           end
           
           
           clear defDef_gb_2 defDef_minAngle1_2 defDef_minAngle2_2
           w2 = waitbar(0, 'Grains Analysed');
           for i=1:length(idB_2(:,1))
               w2 = waitbar(i/length(idB_2), w2)
               [defDef_gb_2{i}, defDef_minAngle1_2{i}, defDef_minAngle2_2{i}, m] = obj.gbdirectionPlanetraceAngle(phase, smoothedGrains('id', idB_2(i,1)), smoothedGrains('id', idB_2(i,2)), identity, false);
                defDef_gb_2_AllGBs{i} = smoothedGrains('id', idB_2(i,1)).boundary(smoothedGrains('id', idB_2(i,2)));

           end
           close(w2);
           
           defDef_gb_2 = cat(1, defDef_gb_2{:});
           defDef_minAngle1_2 = cat(1, defDef_minAngle1_2{:});
           defDef_minAngle2_2 = cat(1, defDef_minAngle2_2{:});
           defDef_gb_2_AllGBs = cat(1, defDef_gb_2_AllGBs{:});
           
%            obj.plotGbdirectionPlanetraceAngle(phase, smoothedGrains('id', idB_2(:,1)), smoothedGrains('id', idB_2(:,2)), defDef_gb_2, defDef_minAngle1_2, defDef_minAngle2_2, m, identity)

           %% Plot fraction
           colorList = getColorList(obj);
           %plotColors = {[1 1 1]; [1 0 0]; 	[1 1 0]; [0.4660 0.6740 0.1880]; [0.4940 0.1840 0.5560]};
           plotColors = {[1 1 1] ; colorList(5,:); 	colorList(7,:); colorList(8,:); colorList(10,:)};
           
           [~, defDef_minAngle1_1_id] =  min(defDef_minAngle1_1, [], 2); %Grains 1
           [~, defDef_minAngle2_1_id] =  min(defDef_minAngle2_1, [], 2); %Grains 1
           [~, defDef_minAngle1_2_id] =  min(defDef_minAngle1_2, [], 2); %Grains 1
           [~, defDef_minAngle2_2_id] =  min(defDef_minAngle2_2, [], 2); %Grains 1
           
           if strcmp(identity, 'incoherent')
               %based on CSL3
               defDef_gb_1_basis = [defDef_gb_1 obj.checkCoherence(phase, defDef_gb_1_AllGBs, obj.CoherentCSL3Threshold, 0)];
               defDef_gb_2_basis = [defDef_gb_2 obj.checkCoherence(phase, defDef_gb_2_AllGBs, obj.CoherentCSL3Threshold, 0)];
               def1_gb_basis = [def1_gb obj.checkCoherence(phase, def1_gb_AllGBs, obj.CoherentCSL3Threshold, 0)];
               def2_gb_basis = [def2_gb obj.checkCoherence(phase, def2_gb_AllGBs, obj.CoherentCSL3Threshold, 0)];
           else
               %based on all
               defDef_gb_1_basis = defDef_gb_1_AllGBs;
               defDef_gb_2_basis = defDef_gb_2_AllGBs;
               def1_gb_basis = def1_gb_AllGBs;
               def2_gb_basis = def2_gb_AllGBs;
           end
           
           gbchar_defDef_2x(1) = sum(obj.checkCoherence(phase, defDef_gb_1_AllGBs, obj.CoherentCSL3Threshold, 0).segLength)/sum(defDef_gb_1_basis.segLength)*100;
           gbchar_defDef_2_sym(1) = gbchar_defDef_2x(1);
           
           gbchar_defDef_11(1) = sum(obj.checkCoherence(phase, defDef_gb_2_AllGBs, obj.CoherentCSL3Threshold, 0).segLength)/sum(defDef_gb_2_basis.segLength)*100;
           gbchar_defDef_12(1) = gbchar_defDef_11(1);
           gbchar_defDef_1_sym(1) = gbchar_defDef_11(1);
           
           for im = 1:(length(m))
               if(~isempty(defDef_gb_1))
                   gbchar_defDef_21(im+1) = sum(defDef_gb_1(defDef_minAngle1_1_id == im).segLength)/sum(defDef_gb_1_basis.segLength)*100;
                   gbchar_defDef_22(im+1) = sum(defDef_gb_1(defDef_minAngle2_1_id == im).segLength)/sum(defDef_gb_1_basis.segLength)*100;
                   gbchar_defDef_2x(im+1) = (sum(defDef_gb_1(defDef_minAngle1_1_id == im).segLength)+sum(defDef_gb_1(defDef_minAngle2_1_id == im).segLength))/(2*sum(defDef_gb_1_basis.segLength))*100;
                   
                   sym = (defDef_minAngle1_1_id == defDef_minAngle2_1_id); %is symmetric?
                   symmetric = sym .* defDef_minAngle1_1_id; %just ids of symmetric GBs
                   gbchar_defDef_2_sym(im+1) = (sum(defDef_gb_1(symmetric == im).segLength)+sum(defDef_gb_1(symmetric == im).segLength))/(2*sum(defDef_gb_1_basis.segLength))*100;;      
                   clear sym symmetric;     
               else
                   gbchar_defDef_21(im+1)=0;
                   gbchar_defDef_22(im+1)=0;
                   gbchar_defDef_2x(im+1)=0;
                   gbchar_defDef_2_sym(im+1)=0;
               end
                     
               if(~isempty(defDef_gb_2))
                   gbchar_defDef_11(im+1) = sum(defDef_gb_2(defDef_minAngle1_2_id == im).segLength)/sum(defDef_gb_2_basis.segLength)*100;
                   gbchar_defDef_12(im+1) = sum(defDef_gb_2(defDef_minAngle2_2_id == im).segLength)/sum(defDef_gb_2_basis.segLength)*100;

                   sym = (defDef_minAngle1_2_id == defDef_minAngle2_2_id); %is symmetric?
                   symmetric = sym .* defDef_minAngle1_2_id; %just ids of symmetric GBs
                   gbchar_defDef_1_sym(im+1) = sum(defDef_gb_2(symmetric == im).segLength)/sum(defDef_gb_2_basis.segLength)*100;      
                   clear sym symmetric;  
               else
                   gbchar_defDef_11(im+1)=0;
                   gbchar_defDef_12(im+1)=0;
                   gbchar_defDef_1_sym(im+1) = 0;
               end
           end
%            figure
%            b1 = bar([gbchar_defDef_12; gbchar_defDef_11; gbchar_defDef_2x], 1);
%            
%            set(gcf, 'Name', ['DrxOnTwins_' phase '_PriorTwinGBPlanes'])
%            str = {'HA-GB to DRX'; 'LA-GB to DRX'; 'Mixed DRX'};
%            set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
%            ylabel('Fraction of GB segments [%]')
%
%            for i=1:length(b1)
%                b1(i).FaceColor = plotColors{i};
%                
%                xtips = b1(i).XEndPoints;
%                ytips = b1(i).YEndPoints;
%                
% %                %Add symmetry info
% %                hold on
% %                s = scatter(xtips, [gbchar_defDef_1_sym(i) gbchar_defDef_1_sym(i) gbchar_defDef_2_sym(i)] ,'.','w')
% %                s.DisplayName = 'Symmetrical';
%                
%                if (i == 1)
%                    b1(i).DisplayName = 'Coherent CSL3';
%                    labels =  'coh.';
%                else
%                    b1(i).DisplayName = char(m(i-1));
%                    labels =  char(m(i-1));
%                    
%                end
%              
%                text(xtips,ytips,labels,'HorizontalAlignment','center',...
%                        'VerticalAlignment','bottom', 'FontSize', 7)
%           end
%           legend show
           
           
           [~, def1_minAngle1_id] =  min(def1_minAngle1, [], 2); %DRX HAGB
           [~, def1_minAngle2_id] =  min(def1_minAngle2, [], 2); %Deformed HAGB
           [~, def2_minAngle1_id] =  min(def2_minAngle1, [], 2); %DRX LAGB
           [~, def2_minAngle2_id] =  min(def2_minAngle2, [], 2); %Deformed LAGB
           
           %HA-GBs
           gbchar_def1_1(1) = sum(obj.checkCoherence(phase, def1_gb_AllGBs, obj.CoherentCSL3Threshold, 0).segLength)/sum(def1_gb_basis.segLength)*100;
           gbchar_def1_2(1) = gbchar_def1_1(1);
           gbchar_def1_sym(1) = gbchar_def1_1(1);
           %LA-GBs
           gbchar_def2_1(1)  = sum(obj.checkCoherence(phase, def2_gb_AllGBs, obj.CoherentCSL3Threshold, 0).segLength)/sum(def2_gb_basis.segLength)*100;
           gbchar_def2_1(isnan(gbchar_def2_1))= 0;
           gbchar_def2_2(1) = gbchar_def2_1(1);
           gbchar_def2_sym(1) = gbchar_def1_1(1);

           for im = 1:length(m)
               if (length(def1_gb) > 0)
                    %For DRX Grain
                    gbchar_def1_1(im+1) = sum(def1_gb(def1_minAngle1_id == im).segLength)/sum(def1_gb_basis.segLength)*100;
                    %For Deformed Grain
                    gbchar_def1_2(im+1) = sum(def1_gb(def1_minAngle2_id == im).segLength)/sum(def1_gb_basis.segLength)*100;
                    
                    sym = (def1_minAngle1_id == def1_minAngle2_id); %is symmetric?
                    symmetric = sym .* def1_minAngle2_id; %just ids of symmetric GBs
                    gbchar_def1_sym(im+1) = sum(def1_gb(symmetric == im).segLength)/sum(def1_gb_basis.segLength)*100;    
                    clear sym symmetric;     
               else
                    gbchar_def1_1(im+1) = 0;
                    gbchar_def1_2(im+1) = 0;
                    gbchar_def1_sym(im+1) = 0;
               end
               
               if (length(def2_gb) > 0)
                   gbchar_def2_1(im+1) = sum(def2_gb(def2_minAngle1_id == im).segLength)/sum(def2_gb_basis.segLength)*100;
                   gbchar_def2_2(im+1) = sum(def2_gb(def2_minAngle2_id == im).segLength)/sum(def2_gb_basis.segLength)*100;
               
                    sym = (def2_minAngle1_id == def2_minAngle2_id); %is symmetric?
                    symmetric = sym .* def2_minAngle2_id; %just ids of symmetric GBs
                    gbchar_def2_sym(im+1) = sum(def2_gb(symmetric == im).segLength)/sum(def2_gb_basis.segLength)*100;                   
                    clear sym symmetric;     
               else
                   gbchar_def2_1(im+1) = 0;
                   gbchar_def2_2(im+1) = 0;
                   gbchar_def1_sym(im+1) = 0;
               end
           end
           f2 = figure;
           %b2 = bar([gbchar_def1_1; gbchar_def1_2; gbchar_def2_1; gbchar_def2_2]);
           b2 = bar([gbchar_def1_1; gbchar_def1_2], 1);
           
           set(gcf, 'Name', ['DrxOnTwins_' phase '_DrxDeformedGBPlanes'])
           str = {'HA-GB DRX'; 'HA-GB Def.'};
           %str = {'HA-GB DRX'; 'HA-GB Def.'; 'LA-GB DRX'; 'LA-GB Def.'};
           set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
           ylabel('Fraction of GB segments [%]')
           
           yl = get(gca, 'Ylim'); %Maybe remove -> mostly for paper!
           if (yl(2) < 55)
               set(gca, 'ylim', [0 55])
           end

           for i=1:length(b2)
               b2(i).FaceColor = plotColors{i};
               
               xtips = b2(i).XEndPoints;
               ytips = b2(i).YEndPoints;
               
%                %Add symmetry info
%                hold on
%                s = scatter(xtips, [gbchar_def1_sym(i) gbchar_def1_sym(i)] ,'.','w')
%                s.DisplayName = 'Symmetrical';
               
               if (i == 1)
                   b2(i).DisplayName = 'Coherent CSL3';
                   labels =  'coh.';
               else
                   b2(i).DisplayName = char(m(i-1));
                   labels =  char(m(i-1));
               end
               
%                text(xtips,ytips,labels,'HorizontalAlignment','center',...
%                         'VerticalAlignment','bottom', 'FontSize', 11)

           end
           legend show

        end
        
        function calcRandomAxisDistribution(obj, phase)
            cs = obj.Ebsd(phase).CS;
            
            %% Random IPF
            figure
            calcAxisDistributionFunction(cs,cs,'antipodal')
            mtexColorbar
            
%             % compute pairs of orientations to be used to compute axis
%             % distributions in specimen coordinates
%             [ori1,ori2] = calcMisorientation(obj.Ebsd(phase))
%             plot(axis(ori1,ori2),'contourf')
            
            
        end
        
        function [x,p] = calcRandomAxisDeviation(obj, axis, normX, normY)
            %% Distribution of closest rotation axis for random orientations
            % Script included to calculate random distribution, but for
            % <111> already calculated & just values used
            %
            % INPUT
            %   axis: rotation axis of interest (e.g. 111)
            %   normX, normY: values to noramlise distribution
            %
            % OUTPUT
            %   x: min. axis deviation angle
            %   p: probability/frequency
            
            
            %Normalization:
            total = trapz(normX, normY);
            
            if (axis == 111)
                 %% Script used to generate random Data:
%                 %Calculated before as it takes a long time!
%                 l = 50000;
%                 odf = uniformODF(cs)
%                 ori1 =  odf.discreteSample(l)
%                 ori2 =  odf.discreteSample(l)
%                 mori = inv(ori1) .* ori2;
% 
%                 minAxis = [];
%                 for i=1:l
%                     minAxis = [minAxis min(angle(mori(i).symmetrise * Miller(1,1,1, cs), Miller(1,1,1, cs))./degree)]; 
%                 end
% 
%                 %Faster but memory intensive:
%                 %minAxis = min(angle(mori.symmetrise .* Miller(1,1,1, cs), Miller(1,1,1, cs)), [], 1)./degree;
% 
%                 histogram(minAxis)
%                 hold on
%                 y = histfit(minAxis, 200, 'kernel')
%                 x = p(2).XData
%                 y = p(2).YData

                  %% Random Data for min<111> Axis deviation
                  x = [0.201821794725010,0.538856016617115,0.875890238509220,1.21292446040133,1.54995868229343,1.88699290418554,2.22402712607764,2.56106134796975,2.89809556986185,3.23512979175396,3.57216401364606,3.90919823553817,4.24623245743027,4.58326667932238,4.92030090121448,5.25733512310659,5.59436934499869,5.93140356689080,6.26843778878290,6.60547201067501,6.94250623256712,7.27954045445922,7.61657467635133,7.95360889824343,8.29064312013554,8.62767734202764,8.96471156391975,9.30174578581185,9.63878000770396,9.97581422959606,10.3128484514882,10.6498826733803,10.9869168952724,11.3239511171645,11.6609853390566,11.9980195609487,12.3350537828408,12.6720880047329,13.0091222266250,13.3461564485171,13.6831906704092,14.0202248923013,14.3572591141934,14.6942933360855,15.0313275579776,15.3683617798697,15.7053960017619,16.0424302236540,16.3794644455461,16.7164986674382,17.0535328893303,17.3905671112224,17.7276013331145,18.0646355550066,18.4016697768987,18.7387039987908,19.0757382206829,19.4127724425750,19.7498066644671,20.0868408863592,20.4238751082513,20.7609093301434,21.0979435520355,21.4349777739276,21.7720119958197,22.1090462177119,22.4460804396040,22.7831146614961,23.1201488833882,23.4571831052803,23.7942173271724,24.1312515490645,24.4682857709566,24.8053199928487,25.1423542147408,25.4793884366329,25.8164226585250,26.1534568804171,26.4904911023092,26.8275253242013,27.1645595460934,27.5015937679855,27.8386279898776,28.1756622117697,28.5126964336618,28.8497306555540,29.1867648774461,29.5237990993382,29.8608333212303,30.1978675431224,30.5349017650145,30.8719359869066,31.2089702087987,31.5460044306908,31.8830386525829,32.2200728744750,32.5571070963671,32.8941413182592,33.2311755401513,33.5682097620434];
                  p = [19.9989983778619,29.5438437749637,40.7980644235201,53.2337306192514,66.2921850330632,79.5017176310369,92.5533270023048,105.339980530299,117.939996127775,130.540104914990,143.290801099921,156.171605777173,168.975628414106,181.462300253301,193.566750803709,205.498113721878,217.600325924686,230.044988573204,242.595082304092,254.660601811778,265.610788550006,275.157183011493,283.520228666289,291.288621187365,299.085766837371,307.269287207189,315.826345314664,324.508760587003,333.072465688041,341.437877196448,349.689500300560,357.952567766839,366.280997129273,374.656861696822,383.033665026392,391.370571323792,399.626627102191,407.726321709210,415.580434334522,423.140649059232,430.444731005858,437.577539780925,444.566619760692,451.310470683175,457.604944808931,463.226241780510,468.024974752324,472.001888811426,475.360190784558,478.455174372820,481.621293738237,484.967426944004,488.292447581874,491.140710508484,492.940874596892,493.114245091593,491.198624022791,487.004078223534,480.696009649962,472.679428984565,463.328615414911,452.723928934308,440.609650269665,426.625566572570,410.611254204565,392.765598567767,373.594726499904,353.728214457827,333.765283295961,314.143699991855,295.058930366064,276.411608770199,257.906353656713,239.296152149786,220.607838847878,202.199249681479,184.609463201469,168.316091644929,153.535776845645,140.178339006782,127.927497307991,116.454102230383,105.591104544306,95.3779234591615,85.9591626149511,77.4401979410320,69.8170353919479,62.9980196557996,56.8399638510479,51.1675169293893,45.7959537959163,40.5864957299092,35.5144781519853,30.6916664733566,26.2983592384533,22.4715252187245,19.2146902756588,16.4000226545657,13.8498742518876,11.4375827817777];
                  p = p/trapz(x,p); %Noramlise to 1!
            end
            
            p = p*total;
        end
        
        function [x,p] = calcRandomMisorientation(obj, normX, normY)
            %% Calculated grain orientation distribution for random cubic grains
            % Mackenzie distribution
            %
            % INPUT
            %   axis: rotation axis of interest (e.g. 111)
            %   normX, normY: values to noramlise distribution
            %
            % OUTPUT
            %   x: min. axis deviation angle
            %   p: probability/frequency
            
            
            %% if other than cube:
            %cs = obj.Ebsd(phase).CS;
            %[p,x] =calcAngleDistribution(cs,cs)
            
            %Normalization:
            total = trapz(normX, normY);
            
            %% if cube:          
            %%Angle in degrees
            x = [0; 
                5
                10
                15
                20
                25
                30
                35
                40
                45
                50
                55
                60
                60.2
                60.4
                60.6
                60.72
                61.0
                61.4
                61.8
                62.2
                62.6
                62.799];
            
            p = [0.00000
                0.00051
                0.00203
                0.00454
                0.00804
                0.01249
                0.01786
                0.02411
                0.03119
                0.03905
                0.03167
                0.02201
                0.01015
                0.00856
                0.00695
                0.00533
                0.00434
                0.00283
                0.00151
                0.00070
                0.00024
                0.00003
                0.00000];
            
            p = p*total;
            
        end
        
        function [coherent, incoherent, noCSL3] = checkCoherence(obj, phase, gb, limCohInco, segMean)
            %% Check <111> Twin Coherency
            % GB plane trace analysis to check if CSL3 segments are coherent
            % twins
            %
            % INPUT
            %   gb: GB segments of interest
            %   limCohInco: max. allowed deviation of GB plane trace to GB
            %       segment direction for it to still be coherent (usually 10
            %       degrees)
            %   segMean: number of segments that should be averaged to get
            %       a direction value (usually 0 if enough GB smoothing
            %       used)
            %
            % OUTPUT
            %   coherent: CSL3 GB segments that are coherent
            %   incoherent: CSL3 GB segments that are incoherent
            %   noCSL3: GB segments that are not CSL3
             
             %% Parameters              
             display('Reminder: apply smoothing to gb (grains.smooth(20)) before applying!');
             
             if (isempty(gb))
                 coherent = [];
                 coherent.segLength = 0;
                 incoherent = [];
                 incoherent.segLength = 0;
                 noCSL3 = [];
                 noCSL3.segLength = 0;
                 return;
             end
             
             noCSL3 = gb(~gb.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree));
             gb = gb(gb.isTwinning(CSL(3,obj.Ebsd(phase).CS), (15/sqrt(3))*degree)); %Just CSL3 (brandon 
             gb = gb('indexed'); %includes all CSL3 segments
             
             if (isempty(gb))
                coherent = gb;
                incoherent = gb;
                return;
             end
             
             %% Trace of <111>
             gbOri = obj.Ebsd('id',gb.ebsdId).orientations;
             
             %if just 1 gb segment -> wrong way around
             if (gbOri.size == [2 1])
                gbOri = transpose(gbOri);
             end
             
             gbAxes = axis(gbOri(:,1),gbOri(:,2), 'antipodal'); %rotation axis in specimen coordinates
             
             %% Trace/Direction of GB 
             %(not GB normal!)
             
             if (segMean == 0)
                 gbDir = gb.direction;
             else
                 gbDir = gb.calcMeanDirection(segMean);
             end
             
             %% Coherent vs. incoherent
             %https://mtex-toolbox.github.io/TiltAndTwistBoundaries.html
             %in mtex wiki: twist/tilt explanation wrong way around?
             
             ang = angle(gbDir, gbAxes)./degree;
             coherent = gb(abs(90-ang) < limCohInco); %90° between GB trace & rotation axis
             incoherent = gb(abs(90-ang) >= limCohInco);       
             
             % Pretty much same results:
%              ang = dot(gbDir, gbAxes)./degree;
%              coherent = gb(ang < limCohInco); %90° between GB trace & rotation axis
%              incoherent = gb(ang >= limCohInco);  

%             % Alternative (same results as first alternative)
%              gbnormal = cross(zvector, gbDir);
%              axp = cross(gbAxes, gbDir);
%              %ang = angle(gbnormal, axp)./degree;
%              ang = dot(gbDir, gbAxes)./degree;
%              
%              coherent = gb(abs(ang) < limCohInco); %90° between GB trace & rotation axis
%              incoherent = gb(abs(ang) >= limCohInco);       

        end
        
        function analyseCoherence(obj, phase)
            %% Analysis of coherent CSL3 segments
            %   Auttomatically smoothes GBs & execudes obj.checkCoherence
            %   and does additional analysis
            %
            % OUTPUT
            %   figure
                         
             %Additional 95 smoothing (= 100 iterations in total)
             gb = obj.Grains(phase).smooth(95).boundary(phase,phase);
             
             [coherent, incoherent, noCSL3] = obj.checkCoherence(phase, gb, obj.CoherentCSL3Threshold, 0);
             
             cohFraction = sum(coherent.segLength)/sum(gb.segLength);
             
             %% Plots    
             f1 = plot(obj, 'phase', 1,0,0)
             hold on
             plot(noCSL3, 'linecolor', 'black', 'lineWidth', 1)
             plot(coherent, 'linecolor', 'white', 'lineWidth', 1.5, 'DisplayName', 'coherent CSL3')
             plot(incoherent, 'linecolor', 'red', 'lineWidth', 1.5, 'DisplayName', 'incoherent CSL3')
             f1.Name = [phase '_TwinCoherency'];
        end
        
        function [gb, minAngle1, minAngle2, m] = gbdirectionPlanetraceAngle(obj, phase, grains1, grains2, identity, opt)
            %% GB plane trace analysis
            % checks for closest plane trace that is parallel to GB
            % direction
            % important: grains1/grains2 should be sufficiently smoothed
            %   before applying function!
            %
            % INPUT
            %     grains1: grain or grains of interest (partition 1)
            %     grains2: grain or grains of interes
            %     identity : 'all', if for all GB between 1/2, 
            %                'incoherent' for just incoherent ones
            %                'allPlusCoherent', same as all but mark coherent CSL3
            %                'allMinusCoherent', mark coh. CSL3 but don't calcualte planes
            %     opt : '1', if plots & documentation shoud be written
            %
            % OUTPUT
            %   gb: GB segments
            %   minAngle1: angle between GB direction & GB plane of grain1 (columns
            %       are m, rows are GB segments)
            %   minAngle2: angle between GB direction & GB plane of grain2 (columns
            %       are m, rows are GB segments)
            %   m: checked Miller indices (defined in function)
                                   
            %% Smooth GBs
            %actually +smooth(5) as grains have already been smoothed
            grains1 = grains1;
            grains2 = grains2;
            
            %% GB Segments!
            gb = grains1.boundary(grains2);          
            gb = gb(phase, phase);
            
            %gb = obj.Twins;
            %apply grains.smooth before as segMean smoothing has issues with assigning GB segments to grains after
            %do not avearage segment length (0) as problems with assigning
            %to GBs after
            
            %defines on which types of GBs (identity) the function should
            %be applied on 
            if (strcmp(identity, 'allPlusCoherent') || (strcmp(identity, 'incoherent')) || (strcmp(identity, 'allMinusCoherent')))
                [coh, inco, noCSL3] = obj.checkCoherence(phase, gb,obj.CoherentCSL3Threshold,0); %threshold between coh/incoherent is 10 --> check if it makes sense!
                
                if (strcmp(identity, 'incoherent'))
                    gb = inco;
                end
                
                if (strcmp(identity, 'allMinusCoherent'))
                    gb = [inco noCSL3];
                end
            end

            %% EBSD Data
            ebsd = obj.Ebsd('id', gb.ebsdId);
            cs = ebsd(phase).CS;
            ori = ebsd.orientations;
            
            %in case there are orientations with just 1 GB segment (puts
            %them wrong way around)
            if (ori.size == [2 1])
                ori = transpose(ori);
            end
           
            %Miller Indecis to check!
            m = Miller({1,1,2},{1,2,2},{1,1,0},cs);
            
            % cycle through all major lattice planes
            minAngle1 = []; %columns are m, rows are GB segments
            minAngle2 = [];
            
            w = waitbar(0,'GB segments analysed'); 
            
            % could be implemented as matrix instead (but uses a lot of
            % memory)
            for igb = 1:length(gb)
                 w = waitbar(igb/length(gb));
                 
                for im = 1:length(m)
                    plane = ori(igb,1) * m(im).symmetrise;
                    %plot(plane, 'upper', 'projection', 'stereo', 'grid', 'grid_res', 5*degree,'DisplayName',char(m(im)))
                    %hold on
                    %legend({},'location','SouthEast','FontSize',13);
                    minAngle1(igb,im) = min(90-angle(gb(igb).direction, plane)./degree); %gb.prop.minAngle 
                    
                    plane = ori(igb,2) * m(im).symmetrise;
                    minAngle2(igb,im) = min(90-angle(gb(igb).direction, plane)./degree); %gb.prop.minAngle 
                                        
                   %minAngle2(igb,im) =min(angle(gb(igb).direction, plane)./degree) %gb.prop.minAngle 
                end
                plane2 = ori(igb,1) * m(2).symmetrise;
                plane3 = ori(igb,1) * m(3).symmetrise;
            end
            
            close(w)
            
            %% Plotting (seperate class!)
            if (opt == true)
               obj.plotGbdirectionPlanetraceAngle(phase, grains1, grains2, gb, minAngle1, minAngle2, m, identity) 
            end
        end
        
        function plotGbdirectionPlanetraceAngle(obj, phase, grains1, grains2, gb, minAngle1, minAngle2, m, identity)
           %% GB plane trace visualisation 
            % visualises results of obj.gbdirectionPlanetraceAngle
            %
            % INPUT
            %     grains1: grain or grains of interest (partition 1)
            %     grains2: grain or grains of interes
            %     identity : 'all', if for all GB between 1/2, 
            %                'incoherent' for just incoherent ones
            %                'allPlusCoherent', same as all but mark coherent CSL3
            %                'allMinusCoherent', mark coh. CSL3 but don't calcualte planes
            %   gb: GB segments
            %   minAngle1: angle between GB direction & GB plane of grain1 (columns
            %       are m, rows are GB segments)
            %   minAngle2: angle between GB direction & GB plane of grain2 (columns
            %       are m, rows are GB segments)
            %   m: checked Miller indices (defined in function)
            
            %% Parameter
            ipf = false; %put IPF behind fraction of grains of which GB plane traces are plotted!
            
            %% Coherence
            if (strcmp(identity, 'allPlusCoherent') || (strcmp(identity, 'incoherent')) || strcmp(identity, 'allMinusCoherent'))
                [coh, inco] = obj.checkCoherence(phase, grains1.boundary(grains2),obj.CoherentCSL3Threshold,0); %threshold between coh/incoherent is 10 --> check if it makes sense!
            end
            
            if (strcmp(identity, 'allPlusCoherent') || (strcmp(identity, 'incoherent')) || strcmp(identity, 'allMinusCoherent'))
                obj.writeReport('Fraction Analysed [%]',sum(inco.segLength)/(sum(coh.segLength)+sum(inco.segLength))*100);
            elseif(strcmp(identity, 'all'))
                obj.writeReport('Fraction Analysed [%]','100');
            end
            
            %% Plot for both grain partitions extra
            for i=1:2
                if(i == 1)
                    [val, id] =  min(minAngle1, [], 2); %Grains 1
                    g = grains1;
                elseif(i == 2)
                    [val, id] =  min(minAngle2, [], 2); %Grains 2
                    g = grains2;
                end
                
                %% Plane with smallest deviation
                f1(i) = figure;
                hold on
                 
                total = [grains1 grains2];
                   plot(total.boundary, 'lineWidth', 1)
                   hold on
                   
                if (ipf == true)
                    total = [grains1 grains2];
                    plot(total.boundary, 'lineWidth', 1)
                    hold on
                    
                    ipfKey = ipfHSVKey(obj.Ebsd(phase).CS);
                    ipfKey.inversePoleFigureDirection = vector3d.Y; %OR Z!
                    
                    ebsd = obj.EbsdOriginalGrains(g);
                    
                    colors = ipfKey.orientation2color(ebsd.orientations);

                    plot(ebsd, colors, 'FaceAlpha', 0.5)
                else
                    plot(g, 'lineWidth', 1)
                end
                hold on
                
                if(isempty(gb))
                   close
                   continue; 
                end
                
                colorList = getColorList(obj);
               %plotColors = {[1 1 1]; [1 0 0]; 	[1 1 0]; [0.4660 0.6740 0.1880]; [0.4940 0.1840 0.5560]};
               plotColors = {colorList(5,:); 	colorList(7,:); colorList(8,:); colorList(10,:)};
           
                if (~isempty(gb(id == 1)))
                    plot(gb(id == 1), 'linecolor', plotColors{1},'linewidth', 3.5,'DisplayName',char(m(1))); %111
                    hold on
                end
                if (~isempty(gb(id == 2)))
                    plot(gb(id == 2), 'linecolor', plotColors{2},'linewidth', 3.5,'DisplayName',char(m(2))); %112
                    hold on
                end
                if (~isempty(gb(id == 3)))
                    plot(gb(id == 3), 'linecolor', plotColors{3},'linewidth', 3.5,'DisplayName',char(m(3)));%122
                    hold on
                end
                if (~isempty(gb(id == 4)))
                    plot(gb(id == 4), 'linecolor', plotColors{4},'linewidth', 3.5,'DisplayName',char(m(4))); %110
                    hold on
                end
                if (strcmp(identity, 'allPlusCoherent') || (strcmp(identity, 'incoherent')) || strcmp(identity, 'allMinusCoherent'))
                   hold on
                   plot(coh, 'linecolor', 'black','linewidth', 3); %111
                   plot(coh, 'linecolor', 'white','linewidth', 2,'DisplayName','Coherent (111) CSL3'); %111
                end
                f1(i).Name = ['GBPlaneTraceAngle_Grain' int2str(i)];                %f1(i).Name = ['GBPlaneTraceAngle_Grain' int2str(i)];
                legend show
                
                %% Plot deviation of each plane along GBs
                for im = 1:length(m)
                    f2(i,im) = figure;
                    plot(g)
                    hold on;
                    plot(gb, minAngle1(:,im), 'linewidth', 3.5,'DisplayName',char(m(im))) %111
                    mtexColorbar('title',['Angle between GB trace and ' char(m(im)) ' trace [°]'])
                    mtexColorMap('LaboTeXColorMap')
                    setColorRange([0 10])
                    f2(i,im).Name = ['GBPlaneTraceAngle_Grain' int2str(i) '_' char(m(im))];
                end
                %% Add to Report
                %Fractions of different m (Automated)
                %Average deviation for differnet m
                                
                obj.writeReport('GB Trace vs. Plane Trace Analysis',['Grain' int2str(i)]);

                 for im = 1:length(m)
                     obj.writeReport(['% GB segments with smallest angle to ' char(m(im))], sum(gb(id == im).segLength)/sum(gb.segLength)*100);
                 end
                 
                 for im = 1:length(m) %plot afterwards!
                     obj.writeReport(['Mean Deviation ' char(m(im))], mean(minAngle1(:,im)));
                     obj.writeReport(['Std Deviation ' char(m(im))], std(minAngle1(:,im)));
                     obj.writeReport(['Min Deviation ' char(m(im))], min(minAngle1(:,im)));
                     obj.writeReport(['Max Deviation ' char(m(im))], max(minAngle1(:,im)));  
                 end
            end
        end
                    
        function correspondingLatticePlane(obj, maxHKL)
            %% NOT YET FINISHED (test)
            
            f = gcf;
            %COLORS = {'darkblue' 'fuchsia' 'purple' 'coral'}; %more?
            COLORS = [];
            EBSD = [];
            ORI = [];
            CS = [];
            n = 2;
            
            for i=1:n
                hold on;
               [x,y] = ginput(1); 
                ebsd = obj.Ebsd(x,y);

                cS = crystalShape.cube(ebsd.CS);
                scaling = 5; % scale the crystal shape to have a nice size

                ipfKey = ipfHSVKey(ebsd);
                COLORS(i,:) = ipfKey.orientation2color(ebsd.orientations);

                plot(x, y, -50, ebsd.orientations * cS * scaling, 'FaceColor', COLORS(i,:), 'FaceAlpha', 0.6); %change

                ori = ebsd.orientations;
                cs = ebsd.CS;
                
                EBSD = [EBSD ebsd];
                CS = [CS cs];
                ORI = [ORI ori];
            end

             ori1 = EBSD(2).orientations;
             ori2 = EBSD(1).orientations;
             mori = inv(ori1) * ori2;
             m = Miller({1,1,1},{1,0,0}, {1,1,2},{1,2,2},cs);
             % cycle through all major lattice planes
             g = figure;
             for im = 1:length(m)
                 % plot the lattice planes of grains 85 with respect to the
                 % reference frame of grain 74
                 plot(mori * m(im).symmetrise,'MarkerSize',10,...
                     'DisplayName',char(m(im)),'figSize','large','noLabel','upper', 'projection','stereo', 'grid', 'grid_res', 5*degree)
                 hold all
             end
             hold off
             % mark the corresponding lattice planes in the twin
             mm = round(unique(mori*m.symmetrise,'noSymmetry'),'maxHKL',maxHKL);
             annotate(mm,'labeled','MarkerSize',5,'figSize','large','textAboveMarker') %TRY: ,'upper', 'projection','stereo', 'grid', 'grid_res', 5*degree
             % show legend
             legend({},'location','SouthEast','FontSize',13);
             hold on
        end
        
        function f = plot(obj, type, gb, iq, partition)
            %% Plot function to create consistant looking plots
            % see MTEX documentation for details on plotting
            %
            % INPUT
            %   type: type of plot created ('phase' or 'ipf' or ''
            %         or 'SubgrainDensity' or 'KAM' or 'GAM' or ...)
            %   gb: '0' no boundaries, '1' just GBs, '2' GBs and Twins,
            %       '3' GBs & Twins % Subgrains
            %   iq: '0' no IQ map in background, '1' IQ map in background
            %   partition: 'char' for phase name
            %              grain2d for grain partition
            %              ebsd for ebsd partition
            % OUTPUT
            %   figure f
            
            % Creates new figure
            alpha = 1;   
            f = figure('Name', 'figure');
            
            %% Define EBSD partition
            ebsdOriginal = obj.EbsdOriginal; %as-recorded EBSD i.e. unindexed points included
            
            originalGrains = obj.Grains; %grains before merging
            ebsdOriginalGrains = obj.EbsdOriginalGrains; %denoised EBSD with grainId from obj.Grains
            
            if strcmp(class(partition), 'char')
                
                ebsd = obj.Ebsd(partition);
                ebsdOriginal = obj.EbsdOriginal(partition); 

                if (isempty(obj.MergedGrains) == 1)
                     grains = obj.Grains(partition);
                else
                     grains = obj.MergedGrains(partition);
                end
                originalGrains = obj.Grains(partition);
                ebsdOriginalGrains = obj.EbsdOriginalGrains(partition);
                
            elseif strcmp(class(partition), 'grain2d')
                grains = partition;
                originalGrains = grains;
                
                ebsd = [];
                
                for i=1:real(length(grains))     
                    ebsd = [ebsd obj.Ebsd(obj.Ebsd.grainId == grains(i).id)];                   
                    hold on; 
                end
                
                if strcmp(type, 'GND')
                   disp('GND cannot be automatically calculated with grain partition');
                   close;
                   return;
                end
            elseif strcmp(class(partition), 'EBSD')
                ebsd = partition;
                
                if (isempty(obj.MergedGrains) == 1)
                     grains = obj.Grains;
                else
                     grains = obj.MergedGrains;
                end
                
                grains = grains(ebsd.grainId)
                
            else
                ebsd = obj.Ebsd;
                if (isempty(obj.MergedGrains) == 1)
                     grains = obj.Grains;
                else
                     grains = obj.MergedGrains;
                end

            end
            
            %% Plotting of IQ Map (in Background)
            if (iq == true)
                try
                    % Oxford System
                    iq = ebsdOriginal.prop.bc;
                catch
                    % Edax System
                    iq  = ebsdOriginal.prop.iq;
                end
                
                %Plot IQ Map
                plot(ebsdOriginal, iq)
                mtexColorMap black2white;
                alpha = 0.6;
                hold on;
                f.Name = 'ImageQuality';
            end
            
            %% Plotting of Phase/IPF Map
            if (strcmp(type, 'phase'))
                % Plotting of Phase Map
                plot(grains, ebsd('indexed'),'noboundary','FaceAlpha',alpha); 
                hold on;
                f.Name = 'Phase';
                
            elseif (strcmp(type, 'ipf'))
                % Plotting of IPF
                
                % gets alls phases saved in EBSD
                a = size(ebsd.mineralList);
                a = a(2);
                
                % runs through all phases, and plot their IPF seperatly
                % the legend can be plot via plot(ipfKey)
                for x = 2:(a+1)
                    % Skips first one because it is 'notIndexed'
                    try
                        hold on;
                        f.Name = 'IPF';
                        phase = ebsd.mineralList{x};
                        
                        ipfKey = ipfHSVKey(ebsd(phase));
                        colors = ipfKey.orientation2color(ebsd(phase).orientations);
                        
                        plot (ebsd(phase), colors,'FaceAlpha',alpha)  
                    catch
                        % in case there is a phase in data where no indexed
                        % points are found
                        continue
                    end
                    
                end
                
                 elseif (strcmp(type, 'ipfY'))
                % Plotting of IPF
                
                % gets alls phases saved in EBSD
                a = size(ebsd.mineralList);
                a = a(2);
                
                % runs through all phases, and plot their IPF seperatly
                % the legend can be plot via plot(ipfKey)
                for x = 2:(a+1)
                    % Skips first one because it is 'notIndexed'
                    try
                        hold on;
                        f.Name = 'IPF_Y';
                        phase = ebsd.mineralList{x};
                        
                        ipfKey = ipfHSVKey(ebsd(phase));
                        ipfKey.inversePoleFigureDirection = vector3d.Y;
                        colors = ipfKey.orientation2color(ebsd(phase).orientations);
                        
                        plot (ebsd(phase), colors,'FaceAlpha',alpha)  
                    catch
                        % in case there is a phase in data where no indexed
                        % points are found
                        continue
                    end
                    
                end
%                 figure;
%                 plot(ipfKey);

            elseif (strcmp(type, 'ipdfY'))
                clf
                %a = size(ebsd.mineralList);
               %a = a(2);
                %for x = 2:(a+1)
                   % % Skips first one because it is 'notIndexed'
                    try
                    plotIPDF(ebsd.orientations,yvector, 'smooth')
                    mtexColorbar
                    catch
                        close
                        f = [];
                        display('Specify Partition!');
                        return
                    end
               % end
                gb = 0;
                f.Name = 'IPF_Y_DensityFunction';

            elseif (strcmp(type, 'SubgrainDensity'))
                hold on
                plot(grains, grains.subBoundarySize ./ grains.area);
                mtexColorbar('title','Subgrain segment number per grain area [\mum^-^2]');
                setColorRange([0,10])
                mtexColorMap LaboTeX
                hold on;
                f.Name = 'SubgrainBoundaryDensityPerGrainArea';
                
            elseif (strcmp(type, 'KAM'))
                
                try
                    hold on
                    plot(ebsd, ebsd.prop.KAM ./degree);
                    caxis([0,2]);
                    mtexColorbar('title','Kernel Average Misorientation (KAM) [°]');
                    mtexColorMap LaboTeX %WhiteJetColorMap  
                    hold on;
                    f.Name = 'KAM';
                catch
                    obj.calculateKAM();
                    close;
                    plot(obj, type, gb, iq, partition)
                    return;
                end
            
            elseif (strcmp(type, 'GROD'))              
                try
                    hold on
                    plot(ebsd, ebsd.prop.GROD.angle./degree,'micronbar','on');
                    mtexColorbar('title','Misorientation angle to grain mean orientation [°]');
                    mtexColorMap LaboTeX %WhiteJetColorMap
                    setColorRange([0,30])

                    hold on;
                    f.Name = 'GROD';
                catch
                    obj.calculateGROD();
                    close;
                    plot(obj, type, gb, iq, partition)
                    return;
                end
 
            elseif (strcmp(type, 'GOS'))
                %GROD = ebsdOriginal.calcGROD(originalGrains);
                %GOS = grainMean(ebsdOriginal, GROD.angle);
                
                GOS = originalGrains.prop.GOS;
                
                %plot(originalGrains, GOS ./ originalGrains.size, 'noboundary')
                % setColorRange([0,0.002])
                hold on
                plot(originalGrains, GOS ./ degree, 'noboundary')
                mtexColorbar('title','Grain orientation spread [°]')
                setColorRange([0,20])
                mtexColorMap LaboTeX
                hold on;
                f.Name = 'GOS';
                
            elseif (strcmp(type, 'CombinedTwinGOS'))
                %https://mtex-toolbox.github.io//GrainMerge.html
                
                GOS = grains.prop.GOS;
                
                %plot(originalGrains, GOS ./ originalGrains.size, 'noboundary')
                % setColorRange([0,0.002])
                hold on
                plot(grains, GOS ./ degree, 'noboundary')
                setColorRange([0,20])
                mtexColorMap LaboTeX
                mtexColorbar('title','GOS [°]')
                hold on;
                f.Name = 'GOS';
                
            elseif (strcmp(type, 'GAM'))
                %EBSD ORIGINAL IS NOT CLEANED! USE DIFFERENT!
                %GAM = ebsdOriginal.grainMean(ebsdOriginal.gridify.KAM('threshold',2.5*degree,'order',1), originalGrains);                
                try
                    hold on
                    plot(grains, grains.prop.GAM ./degree, 'noboundary')
                    %plot(obj.MergedGrains, GAM ./degree, 'noboundary')
                    mtexColorbar('title','Grain Average Minsorientation [°]')
                    setColorRange([0,2])
                    mtexColorMap LaboTeX

                   % histogram(GAM./degree)
                    hold on;
                    f.Name = 'GAM';
                catch
                    obj.calculateKAM();
                    close;
                    plot(obj, type, gb, iq, partition)
                    return;
                end
                
                
             elseif (strcmp(type, 'GND'))
                
                if (isempty(obj.GND) == 1)
                    obj.calculateGND()                    
                end
                
                for i=1:(length(ebsd.gridify.mineralList)-1)
                    gnd = obj.GND{i};
                    
                    hold on;
                    try
                        plot(ebsd(ebsd.phase==i).gridify,gnd,'micronbar','off')
                        hold on;
                    end
                end
                
                %plot(ebsd.gridify,obj.GND,'micronbar','off')
                mtexColorMap LaboTeX
                mtexColorbar('title','Geometrically Necessary Dislocations [m^-^2]')

                set(gca,'ColorScale','log'); % this works only starting with Matlab 2018a
                set(gca,'CLim',[1e13 5e16]);
                
                hold on;
                f.Name = 'GND';
                
                % MORE INFO: https://mtex-toolbox.github.io/GND.html    
                
            elseif(strcmp(type, 'SF'))
                %Schmidfactor
                sS = slipSystem.fcc(ebsd.CS);
                sS = sS.symmetrise;

                sSLocal = grains.meanOrientation * sS;

                sigma = stressTensor.uniaxial(vector3d.Y);
                SF = sSLocal.SchmidFactor(sigma);
                [SFMax,active] = max(SF,[],2);

                % Move to plot function:
                plot(grains,SFMax,'noboundary')
                set(gca,'CLim',[0.3 0.5]);
                mtexColorbar('title','Schmid Factor [-]')
                hold on;
                f.Name = 'SF';
%                 mtexColorbar location southoutside
% 
%                 sSactive = grains.meanOrientation .* sS(active);
% 
%                 hold on
%                 % visualize the trace of the slip plane
%                 quiver(grains,sSactive.trace,'color','b')
% 
%                 % and the slip direction
%                 quiver(grains,sSactive.b,'color','r')
%                 hold off
                
            end
            
            
            %% Plotting of GBs
            if (gb >= 1)       
                % Plotting of GBs
                plot(grains.boundary,'lineWidth',1,'DisplayName','GBs');
                
                if (gb == 2 || gb == 3)
                    % Plotting of Twins
                    plot(obj.Twins,'lineColor','white','linewidth',1,'DisplayName','CSL 3');
                    %plot(obj.Twins,'lineColor','black','linewidth',2.5,'DisplayName','CSL 3');
                    %plot(obj.Twins,'lineColor','white','linewidth',2,'DisplayName','CSL 3');

                    if (~isempty(obj.TwinsDDRX))
                        plot(obj.TwinsDDRX, 'linecolor', 'black', 'linewidth', 1.5)
                        plot(obj.TwinsDDRX, 'linecolor', 'white', 'linewidth', 0.1)
                    end
                end
                
                if (gb == 3)
                   % Plotting of Subgrains
                   
                    % gets alls phases saved in EBSD
                    a = size(ebsd.mineralList);
                    a = a(2);

                    % runs through all phases, and plot their subgrains
                    % seperatly 
                    for x = 2:(a+1)
                        % Skips first one because it is 'notIndexed'
                        try
                            hold on;
                            phase = ebsd.mineralList{x};

                            alpha = grains(phase).innerBoundary.misorientation.angle / (5*degree);
                            
                            % Colour all SB-Boundaries in blue
                            plot(grains(phase).innerBoundary,'linewidth',0.5,'edgeAlpha',alpha,'linecolor','b');
                            
                            % Red = potential tilt boundaries; blue = for sure no tilt 
                            %plot(grains(phase).innerBoundary,angle(grains(phase).innerBoundary.direction,axS)./degree,'linewidth',0.5,'edgeAlpha',alpha)
                            %mtexColorMap blue2red
                            %mtexColorbar
                            
                        catch
                            % in case there is a phase in data where no indexed
                            % points are found
                            continue
                        end
                    end
                   
                 end
            end
            
            %% Micronbar settings
%             f=gcm;
%             mp=getappdata(f.currentAxes,'mapPlot');
%             mp.micronBar.backgroundAlpha=1;
%             mp.micronBar.backgroundColor = [1 1 1]
%             mp.micronBar.lineColor = [0 0 0];
%             mp.micronBar.txt.FontName = 'Arial';
%             mp.micronBar.txt.comiler = 'none';
%             %mp.micronBar.txt.FontWeight = 'normal';
     
            hold off;

        end
               
        function grainD = grainSize(obj, phase)
            %% Determines equivalent diameter of Grains
            % INPUT
            %   phase: phase which should be analysed
            %
            % OUTPUT
            %   histogram
            %   grainD: equivalent diameter
            
            % Binning for Histograms
            bin = 20;
            
            %% Grain Statistics
            % Just considers grains on the inside for grains statistics
            % using the excludeEdge function
            innerGrains = excludeEdge(obj);
            innerGrains = innerGrains(phase);
            
            % Determines equiv. diameter of grains
            grainD =  2*equivalentRadius(innerGrains);
            
            % Number of detected grains
            number = size(grainD);
            number = number(1);
            
            for (i = 1:number)
               if (grainD(i) < obj.minGrainDiameter)
                  grainD(i) = false; 
               end
            end
            grainD = abs(grainD(grainD~=0));
            
            % Average of equiv. diameter
            d = mean(grainD);
            % Standard deviation of equiv. diameter
            dStd = std(grainD);
            % Volume fraction of phase
            vol = sum(area(obj.Grains(phase)))/sum(area(obj.Grains));
            
            % Write determined values in Report
            obj.writeReport('-','-')
            obj.writeReport('Grain Statistic', phase)
            obj.writeReport('Number Grains', number)
            obj.writeReport('Equiv. diameter', d)
            obj.writeReport('Equiv. diameter std.', dStd)
            obj.writeReport('Median diameter', median(grainD))
            obj.writeReport('Volume fraction', vol)
            
            
            %% Create Histograms
            
            % Saves detected grain diameters ro GrainStatistics for
            % analysis in different software
            [tf, phaseIndex] = ismember(phase, obj.Ebsd.mineralList);
            obj.GrainStatistics{phaseIndex} = grainD;
                        
            % Create histogram:
            f = figure('Name', ['GS_' phase]);
            h = histogram(grainD,bin);
            title( join(['Grain Diameter Distribution: ',phase]) );
            hold on  
            xline(d, 'LineWidth', 2, 'Color', 'r');
            xline(median(grainD), 'LineWidth', 2, 'Color', 'y');
            xlabel('Grain diameter [\mum]');
            ylabel('Frequency [-]');
            hold off
        end
                
        function f = boundaryChar(obj, grains, ph1, ph2)
            %% Determines the frequency of boundaries between phases
            % Misorientation angle between the defined partitions
            %
            % INPUT
            %   grains: grains of interest
            %   ph1: First phase
            %   ph2: Second phase
            %
            % OUTPUT
            %   histogram f
            
            % Detect boundary angles
            AOboundary = grains.boundary(ph1,ph2)
            Angle = AOboundary.misorientation.angle;
            
            % Plot in Histogram
            f = figure('Name', ['BoundaryCharacter_' ph1 '_' ph2]);
            histogram(Angle./degree, 'Normalization','probability')
            xlabel([' Misorientation [°]'])
            ylabel('Relative Frequency')
        end
        
        function writeReport(obj, property, value) 
            %% Add comments to Report property in style (property: value)
            % in style "property: value"
            %
            % INPUT
            %   property: name of property
            %   value: value of property
            
           if (isnan(value))
               value = 'NaN';
           end
            
            property = string(property);
            value = string(value);
            
            obj.Report = obj.Report + newline + property + ': ' + value;
        end
        
        function save(analysis, dic)
            %% Saves all open figures, report and grain statistics
            % Figures are saved twice (1x with legend, 1x without)
            % Exports: FunctionsOutline_Combined.m + ebsd_analysis.m
            %
            % INPUT
            %   dic: dictionary it should save the files to
            %   analysis: usually called 'obj' in other functions, here
            %           different to be able to export/import it easier
            
            %% Option to safe fig files
            fig = false; %if true saves .fig (large files!)
                        
            %% Create new folder with date
            dic = fullfile(dic, [date '_Analysis']);
            mkdir (dic)
            
            dicTiff = fullfile(dic, 'TIFF_Figures');
            dicLabels = fullfile(dic, 'TIFF_Labels');
            
            mkdir (dicTiff)
            mkdir (dicLabels)
            
            if (fig == true)
                dicFig = fullfile(dic, 'FIG_Figures');
                mkdir (dicFig)
            end
            
            %% Name of files is based on property Name
            name = sprintf(analysis.Name, 's');
            name = urlencode(name);
            
            %% Save Report
            % Opens file
            fid = fopen( fullfile (dic, sprintf('%s_Report.txt', name)) ,'a');
            % Save Report to file
            fprintf(fid, '%s\n%s\n%s', analysis.Report);
            % Close file
            fclose(fid);
            
            %% Save all open Figures
            % get all figures
            figHandles = findall(0,'Type','figure');
            legHandles = findall(0,'Type','legend'); 
            
            % Loop through figures 2:end and saves them
            for i = 1:numel(figHandles)
               
                figName = [int2str(i) '_' name '_' figHandles(i).Name];
                               
               %% Save with Label
                try
                    legHandles(i).Visible = 'on'; %legHandles(i).Visible = 'on';
                end

                print(figHandles(i), fullfile(dicLabels, [figName '.tif']), '-dtiff', '-r300')
                
                 %% Save as .fig (can be edited later e.g. micronbar removed)
                 if (fig == true)
                     %Very large files! Maybe enough to save .mat file?
                     saveas(figHandles(i), fullfile(dicFig, [figName '.fig']));
                 end
                
                %% Save without Label
                try
                    legHandles(i).Visible = 'off';
                    micronHandle = getappdata(figHandles(i).CurrentAxes,'mapPlot');
                    micronHandle.micronBar.visible = 'off';
                end
                print(figHandles(i), fullfile(dicTiff, [figName '.tif']), '-dtiff', '-r300')
                
                close(figHandles(i))

            end
            
            %% Saves Grain Statistics
            % Checks if grainSize() has been applied
            if (isempty(analysis.GrainStatistics) == 0)
                
                % Checks for how many phases available and loops through
                for i = 1:numel(analysis.GrainStatistics)
                    grains = analysis.GrainStatistics{i};
                    
                    % checks if grain statistics available for that phase
                    if (isempty(grains) == 0)
                        phase = string(analysis.Ebsd.mineralList(i));
                        
                        % Saves grain statistics for that phase
                        fid = fopen( fullfile (dic, sprintf('%s_EquivGrainDiam4eters.txt', phase)) ,'a');
                        for k = 1:numel(grains)
                            fprintf(fid, '%s\n%s\n%s', grains(k));
                        end
                        
                        fclose(fid);
                    end
                end
                                
            end
            
           %% Export ebsd_execution.m + ebsd_analysis.m + workspace
           % Export FunctionsOutline_Combined
           copyfile('ebsd_execution.m', fullfile(dic, 'ebsd_execution.m'))
           
           % Export current ebsd_analysis version
           copyfile('ebsd_analysis.m', fullfile(dic, 'ebsd_analysis.m'))
           
%            %% Export workspace variables (too large!)
%            try
%                save( fullfile(dic, [name '_analysisVariable.mat']), 'analysis','-v7.3');
%            catch
%                display('.mat file could not be saved')
%            end
        end
        
        function innerGrains = excludeEdge(obj)
            %% Exclude all grains which are on the edge of Scan
            % OUTPUT
            %   innerGrains: all grains that do not lay on the edge of scan
            
            % Check if Twins have been detected
            if (isempty(obj.MergedGrains) == 1)
                grains = obj.Grains;
            else
                grains = obj.MergedGrains;
            end
            
            % Checks if Grains lay on outer boundary and removes them
            outerBoundary_id = any(grains.boundary.grainId==0,2);
            gId = grains.boundary(outerBoundary_id).grainId;
            gId(gId==0) = [];
            grains(gId) = [];
            innerGrains = grains;
        end
        
        function cleanEDS(obj, matrix, maxMatrixSize, angle, prec)
              %% Clean EDS
              %     Removes the wrongly assignd matrix grains on the edge
              %     of precipitates due to the larger interaction volume of
              %     EDS compared to EBSD
              %
              % INPUT
              %     matrix: name of matrix phase
              %     maxMatrixSize: limit of larges matrix grain which
              %                    should be considered
              %     angle: limit of max. angle difference between matrix
              %            and precipitate phase to be removed
              %     prec: precipitate phase
              %
              %     If matrix and prec are entered the other way around it
              %     removes secondary gamma prime!
              %     createGrains() must be used before!
              %     use before denoising!
              %
              % OUTPUT
              %     redefines obj.Grains, obj.Ebsd
              %     figure
              
             %% Select relevant boundaries
             % Selects all GBs between matrix and precipitate where matrix
             % is smaller than limit
             g = obj.Grains(obj.Grains.grainSize < maxMatrixSize);
             gB = g(matrix).boundary(matrix, prec);
              
             % Selects all GBs where orientation difference is smaller than angle
             gB = gB(gB.misorientation.angle < angle*degree)

             %% Change Phase of detected matrix grains
             % detect affected grains
             changedGrainIds = unique(gB.grainId(:,1));
             % determine phase number of entered precipitate phase
             [tf, precIndex] = ismember(prec, obj.Ebsd.mineralList);
             
             % Counting variable for how many pixel have been changed
             pixelPhaseChanged = 0;
             pixelTotalMap = size(obj.Ebsd);
             pixelTotalMap = pixelTotalMap(1);
             pixelIndexed = size(obj.Ebsd(obj.Ebsd.isIndexed == 1));
             pixelIndexed = pixelIndexed(1);
             
             % loop through all affected grains
             for (i = 1:size(changedGrainIds))
                 % id of grain
                 id = changedGrainIds(i,1);
                 % change phase of detected matrix to precipitate phase
                 obj.Ebsd(obj.Ebsd.grainId == id).phaseId = precIndex;
                 
                 % increase counting variable by amount of pixel changed
                 pixel = size(obj.Ebsd(obj.Ebsd.grainId == id));
                 pixel = pixel(1);
                 
                 pixelPhaseChanged = pixelPhaseChanged + pixel;
             end
             
             %% Re-calculate Grains after Merging again
             obj.createGrains(obj.minGrainAngle, obj.minGrainSize);
             %% add to Report 
             obj.writeReport('EDS Clean-UP performed with following parameters', '')
             obj.writeReport('Matrix phase', matrix)
             obj.writeReport('Precipitate phase', prec)
             obj.writeReport('Max. size of considered matrix grains', maxMatrixSize)
             obj.writeReport('Max. angle difference between phases', angle)
             obj.writeReport('Phase of Pixel changed', '');
             obj.writeReport('% of total map', pixelPhaseChanged/pixelTotalMap*100);
             obj.writeReport('% of indexed pixel', pixelPhaseChanged/pixelIndexed*100);
             obj.writeReport('-','-')

             %% PLOT CONTROL
             f = figure('Name', 'EDS_CleanUp')
             plot(obj.Grains)
             hold on
             plot(obj.Grains.boundary,'lineColor','black','linewidth',1,'DisplayName','GBs');
             plot(gB,'lineColor','red','linewidth',1.5,'DisplayName','Corrected EDS/EBSD');
             f.Name = 'EDS_CleanUp';
             hold off
        end
        
        function highlightHREX(obj, matrix, minMatrixSize, angle, prec)
            %% Similar to cleanEDS just max<>min Matrix size changed
            
            %% cleanEDS with all grains smaller than minMatrixSize
%             obj.cleanEDS(matrix, minMatrixSize, angle, prec)
            
            %% Select relevant boundaries
             % Selects all GBs between matrix and precipitate where matrix
             % is smaller than limit
             g = obj.Grains(obj.Grains.grainSize > minMatrixSize);
             gBTotal = g(matrix).boundary(matrix, prec);
              
             % Selects all GBs where orientation difference is smaller than angle
             gB = gBTotal(gBTotal.misorientation.angle < angle*degree)
             
             % add to Report 
             obj.writeReport('HREX Grains Highlited', '')
             obj.writeReport('Matrix phase', matrix)
             obj.writeReport('Precipitate phase', prec)
             obj.writeReport('Min. size of considered matrix grains', minMatrixSize)
             obj.writeReport('Max. angle difference between phases', angle)
             obj.writeReport('% matrix/prec GB segments detected', gB.length/gBTotal.length);
             obj.writeReport('-','-')
             % Add Report for deformation parameters of HREX grains
             
              %% PLOT CONTROL
             %f = figure('Name', 'HREX_Highlighted_IPF')
             %hold on;
             obj.plot('ipf', 2,0,0)
             hold on;
             plot(gB,'lineColor','green','linewidth',1.5,'DisplayName','HREX GBs');
             f = gcf;
             f.Name = 'HREX_Highlighted_IPF';
             hold off
             
             %f = figure('Name', 'HREX_Highlighted_phase')
             %hold on;
             obj.plot('phase', 2,0,0)
             hold on;
             plot(gB,'lineColor','green','linewidth',1.5,'DisplayName','HREX GBs');
             f = gcf;
             f.Name = 'HREX_Highlighted_phase';
             hold off
             
              %% Deformation analysis
              classifyDRX(obj, prec, matrix)

        end
                
        function classifyDRX(obj, particle, matrix)
            %% Classification of Rx grains into partitions
            % does not just analyse PSN but classifies Rx grains
            %
            % INPUT
            %   particle: PSN particles of interest (put '' if none)
            %   matrix: name of matrix phase that recrystallises
                        
            %% If RX Function has not been applied yet, apply it
            % Seperate Recrystallised/Deformed Grains
            if (isempty(obj.RxGrains) == true && isempty(obj.DeformedGrains) == true)
                obj.recrystallisedGrains();
            end
            
            %% Section in Analysis Report
            obj.writeReport('-','-')
            obj.writeReport('PSN Matrix', matrix)
            obj.writeReport('PSN Particle', particle)

            if(~strcmp(particle, ''))

                %% Select matrix grains which border particle
                selectedGrains = neighbouringGrains(obj, particle, matrix);
                
                if (~isempty(selectedGrains))
                    %% Select PSN-Rx Grains
                    threshold = obj.RxThreshold;
                    
                    % Save all particle grains to later add to
                    particleGrains = selectedGrains(particle);

                    %% Select matrix Rx grains which border particle + all particles
                    selectedGrains(particle).prop.GOS = 0; % so they are not picked up by next line
                    PSNGrains = selectedGrains(selectedGrains.prop.GOS./degree < threshold);  
                    PSNGrains(particle) = particleGrains;
                    obj.PSNGrains = PSNGrains(matrix);

                    %% Select deformed matrix grains which border particle + all particles
                    selectedGrains(particle).prop.GOS = threshold;
                    noPSNGrains = selectedGrains(selectedGrains.prop.GOS./degree >= threshold);
                    noPSNGrains(particle) = particleGrains; 

                    %% DRX Grains without PSN
                    DRXGrains = obj.RxGrains(matrix);
                    PSNGrainsMatrix = PSNGrains(matrix);

                    DRXGrains('id', PSNGrainsMatrix.id) = [];
                    obj.DRXGrains = DRXGrains;
                else
                    obj.PSNGrains = [];
                    DRXGrains = obj.RxGrains(matrix);
                    obj.DRXGrains = DRXGrains;
                end
            else
                    obj.PSNGrains = [];
                    DRXGrains = obj.RxGrains(matrix);
                    obj.DRXGrains = DRXGrains;
            end
            
            %% DRX Grains within grains (not on GBs!)
            %Some particles migght be wrongly indixed:
            insideRxGrains = obj.DRXGrains(obj.DeformedGrains.checkInside(obj.DRXGrains)); %DRX (not Rx) because PSN are already disgarded
            obj.DRXGrainsWithinGrain = insideRxGrains(matrix);
   
            %% DRX Grains on HAGBs
            obj.DRXGrainsOnGBs = obj.DRXGrains(matrix);
            
            %Delte GB-Grain Fraction out of Fraction:
            obj.DRXGrainsOnGBs('id', obj.DRXGrainsWithinGrain.id) = [];
            
            %% DRX Grains on prior-TWIN GBs
            if (~isempty(obj.DRXGrainsWithinGrain))
                obj.DRXGrainsOnPriorTwinGBs = insideRxGrains(insideRxGrains.prop.InclusionBeforeMerge == 0);
                obj.DRXGrainsWithinGrain = insideRxGrains(insideRxGrains.prop.InclusionBeforeMerge == 1);
            else
                obj.DRXGrainsOnPriorTwinGBs = obj.DRXGrainsWithinGrain;
                obj.DRXGrainsWithinGrain = obj.DRXGrainsWithinGrain;
            end
            
            %% Visualise different Fractions
            f = obj.plot('phase', 2, 0, 0);
            f.Name = (['Phase_PSN_' matrix '_around_' particle]);
            hold on
            plot(obj.DRXGrainsOnGBs(matrix).boundary,'lineWidth',1,'linecolor',1/255*[0 158 115],'DisplayName','DRX Grains (on GBs)');
            plot(obj.DRXGrainsOnPriorTwinGBs(matrix).boundary,'lineWidth',1,'linecolor',1/255*[230 159 0],'DisplayName','DRX Grains (on prior Twins)');
            plot(obj.DRXGrainsWithinGrain(matrix).boundary,'lineWidth',1,'linecolor',1/255*[240 228 66],'DisplayName','DRX Grains (not on GBs)');
            if (~isempty(obj.PSNGrains))
                plot(obj.PSNGrains(matrix).boundary,'lineWidth',1,'linecolor',1/255*[213 94 0],'DisplayName','PSN Grains');
            end
                %plot(obj.TwinsDDRX, 'linecolor', 1/255*[0 158 115], 'linewidth', 2)
            if (~isempty(obj.TwinsDDRX))
                plot(obj.TwinsDDRX, 'linecolor', 'white', 'linewidth', 0.1)
            end

            hold off
            
            %% Fraction of PSNGrains in all RxGrains
            if (~isempty(obj.PSNGrains))
                volFraction = sum(abs(PSNGrains(matrix).area))/sum(abs(obj.RxGrains(matrix).area));
                numFraction = length(PSNGrains(matrix))/length(obj.RxGrains(matrix));

                obj.writeReport('PSN volumne Fraction of all Rx grains', volFraction)
                obj.writeReport('PSN number Fraction of all Rx grains', numFraction)

                %% Average Size PSN GRains
                diameterPSN = mean(abs(PSNGrains(matrix).equivalentRadius))*2;
                stdDiameterPSN = std(abs(PSNGrains(matrix).equivalentRadius))*2;
                medianDiameterPSN = median(abs(PSNGrains(matrix).equivalentRadius))*2;

                gosPSN = mean(PSNGrains(matrix).prop.GOS./degree);
                stdGosPSN = std(PSNGrains(matrix).prop.GOS./degree);

                obj.writeReport('Number of PSN Grains', length(PSNGrains(matrix)))
            
                obj.writeReport('PSN average diameter [um]', diameterPSN)
                obj.writeReport('PSN diameter std', stdDiameterPSN)
                 obj.writeReport('PSN median diameter [ym]', medianDiameterPSN)
                obj.writeReport('PSN average GOS [°]', gosPSN)
                obj.writeReport('PSN GOS std', stdGosPSN)
            end
            obj.writeReport('Number of DRX Grains', length(DRXGrains(matrix)))
            obj.writeReport('DRX average diameter [um]', mean(abs(DRXGrains(matrix).equivalentRadius))*2)
            obj.writeReport('DRX diameter std', std(abs(DRXGrains(matrix).equivalentRadius))*2)
            obj.writeReport('DRX median diameter [um]', median(abs(DRXGrains(matrix).equivalentRadius))*2)
            obj.writeReport('DRX average GOS [°]', mean(DRXGrains(matrix).prop.GOS./degree))
            obj.writeReport('DRX GOS std', std(DRXGrains(matrix).prop.GOS./degree))
            
            obj.writeReport('Number of DRX (Necklace) Grains', length(obj.DRXGrainsOnGBs(matrix)))
            if (length(obj.DRXGrainsOnGBs(matrix)) > 0)
                obj.writeReport('DRX (Necklace) average diameter [um]', mean(abs(obj.DRXGrainsOnGBs(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (Necklace) diameter std', std(abs(obj.DRXGrainsOnGBs(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (Necklace) median diameter [um]', median(abs(obj.DRXGrainsOnGBs(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (Necklace) average GOS [°]', mean(obj.DRXGrainsOnGBs(matrix).prop.GOS./degree))
                obj.writeReport('DRX (Necklace) GOS std', std(obj.DRXGrainsOnGBs(matrix).prop.GOS./degree))
                obj.writeReport('DRX (Necklace) vol. fraction of DRX', sum(abs(obj.DRXGrainsOnGBs(matrix).area))/sum(abs(obj.DRXGrains(matrix).area)) )
            end
                
            obj.writeReport('Number of DRX (not on GBs) Grains', length(obj.DRXGrainsWithinGrain(matrix)))
            if (length(obj.DRXGrainsWithinGrain(matrix)) > 0)
                obj.writeReport('DRX (not on GBs) average diameter [um]', mean(abs(obj.DRXGrainsWithinGrain(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (not on GBs) diameter std', std(abs(obj.DRXGrainsWithinGrain(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (not on GBs) median diameter [um]', median(abs(obj.DRXGrainsWithinGrain(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (not on GBs) average GOS [°]', mean(obj.DRXGrainsWithinGrain(matrix).prop.GOS./degree))
                obj.writeReport('DRX (not on GBs) GOS std', std(obj.DRXGrainsWithinGrain(matrix).prop.GOS./degree))
                obj.writeReport('DRX (not on GBs) vol. fraction of DRX', sum(abs(obj.DRXGrainsWithinGrain(matrix).area))/sum(abs(obj.DRXGrains(matrix).area)) )
            end
            
            obj.writeReport('Number of DRX (Twin GBs) Grains', length(obj.DRXGrainsOnPriorTwinGBs(matrix)))
            if (length(obj.DRXGrainsOnPriorTwinGBs(matrix)) > 0)
                obj.writeReport('DRX (Twin GBs) average diameter [um]', mean(abs(obj.DRXGrainsOnPriorTwinGBs(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (Twin GBs) diameter std', std(abs(obj.DRXGrainsOnPriorTwinGBs(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (Twin GBs) median diameter [um]', median(abs(obj.DRXGrainsOnPriorTwinGBs(matrix).equivalentRadius))*2)
                obj.writeReport('DRX (Twin GBs) average GOS [°]', mean(obj.DRXGrainsOnPriorTwinGBs(matrix).prop.GOS./degree))
                obj.writeReport('DRX (Twin GBs) GOS std', std(obj.DRXGrainsOnPriorTwinGBs(matrix).prop.GOS./degree))
                obj.writeReport('DRX (Twin GBs) vol. fraction of DRX', sum(abs(obj.DRXGrainsOnPriorTwinGBs(matrix).area))/sum(abs(obj.DRXGrains(matrix).area)) )
            end
            
            %% Histogram: Boundary character between particle & matrix grains             
            %figureBoundaryChar = obj.boundaryChar(PSNGrains(matrix), particle, matrix)
            if (isempty(obj.PSNGrains))
               obj.writeReport('No PSN detected', '')
               return
            end
            AOboundaryPSN = PSNGrains(matrix).boundary(particle,matrix);
            AnglePSN = AOboundaryPSN.misorientation.angle ./degree;
            
            AngleNoPSN = [];
            if (length(noPSNGrains(matrix)) > 0)
                AOboundaryNoPSN = noPSNGrains(matrix).boundary(particle,matrix);
                AngleNoPSN = AOboundaryNoPSN.misorientation.angle ./degree;
            end
            
            % Plot in Histogram
            step = 5;
            fBoundaryChar = figure('Name', ['PSN_MisorAng_of_' matrix '_around_' particle]);
            fBoundaryChar = histogram(AnglePSN, 0:step:62, 'Normalization','probability');
            hold on
            fBoundaryChar = histogram(AngleNoPSN, 0:step:62, 'Normalization','probability'); 
            legend('PSN Grains', 'Deformed Grains')
            xlabel([matrix ' misorientation to ' particle ' [°]'])
            ylabel('Relative Frequency')
            hold off
            
            f = figure('Name', ['PSN_MisorAxis_of_' matrix '_around_' particle]);
            plot(AOboundaryPSN.misorientation.axis, 'fundamentalRegion','MarkerSize', 1, 'DisplayName', 'PSN Grains')
            hold on
            plot(AOboundaryNoPSN.misorientation.axis, 'fundamentalRegion', 'MarkerSize',1, 'DisplayName', 'Deformed Grains')
            f.Name = ['PSN_MisorAxis_of_' matrix '_around_' particle];
            legend show
            
            %% Fraction of particle/matrix boundary length, which is recrystallised
            boundaryPSN = PSNGrains(matrix).boundary(particle, matrix); %(matrix) important, so that just border with matrix considered, not entire TiC
            boundaryNoPSN = noPSNGrains(matrix).boundary(particle, matrix);
            
            boundaryFraction = sum(boundaryPSN.segLength)/( sum(boundaryPSN.segLength)+sum(boundaryNoPSN.segLength) );
            
            obj.writeReport('Boundary length fraction of PSN Grains around particles', boundaryFraction)
            
            %obj.Grains.prop.GOS;
         
            %% TiC-GROESSE vs. Segment Fraction Histogram
            maxParticleD = max(particleGrains(particle).equivalentRadius)*2;
            minParticleD = min(particleGrains(particle).equivalentRadius)*2;
            
            for i=1:length(particleGrains(particle))
                particleGrainsBinning = particleGrains(i);
                
                %Segment Length as Parameter
                segLengthPSN = sum(boundaryPSN(particleGrainsBinning).segLength);
                segLengthNoPSN = sum(boundaryNoPSN(particleGrainsBinning).segLength);
                
                particleDiameterVsPSN(i,1) = abs(particleGrainsBinning.equivalentRadius)*2;
                particleDiameterVsPSN(i,2) = segLengthPSN/(segLengthPSN+segLengthNoPSN);
                
                %GOS of particle
                particleGOS = particleGrainsBinning.prop.GOS ./degree;
                particleDiameterVsPSN(i,3) = particleGOS;
                
                %% average GOS of matrix grains arouns (psn)
                PSNGrainsBinning = PSNGrains(matrix);
                neighbors = unique(boundaryPSN(particleGrainsBinning).grainId);
                               
                nGOS = [];
                for j=1:length(neighbors)
                   
                   n = neighbors(j);                    
                   g = PSNGrainsBinning(PSNGrainsBinning.id == n);
                   
                   nGOS = [nGOS; g.prop.GOS ./degree];                 
                   
                end
                
                particleDiameterVsPSN(i,4) = mean(nGOS);
                particleDiameterVsPSN(i,5) = std(nGOS);
                particleDiameterVsPSN(i,6) = length(neighbors);
                
                
                %% average GOS of matrix grains arouns (NO PSN
                NoPSNGrainsBinning = noPSNGrains(matrix);
                neighborsNoPSN = unique(boundaryNoPSN(particleGrainsBinning).grainId);
                               
                nGOS = [];
                for j=1:length(neighborsNoPSN)
                   
                   n = neighborsNoPSN(j);                    
                   g = NoPSNGrainsBinning(NoPSNGrainsBinning.id == n);
                   
                   nGOS = [nGOS; g.prop.GOS ./degree];                 
                   
                end
                
                particleDiameterVsPSN(i,7) = mean(nGOS);
                particleDiameterVsPSN(i,8) = std(nGOS);
                particleDiameterVsPSN(i,9) = length(neighborsNoPSN);
                
            end
            particleDiameterVsPSN = sortrows(particleDiameterVsPSN, 1);
            
            %% Histogram: Size of Particles which lead to PSN  
            % Plot in Histogram
            step = 1;
            
            maxRxGrain = max(obj.RxGrains(matrix).equivalentRadius);
            minRxGrain = min(particleGrains(particle).equivalentRadius);
            
            fParticleSize = figure('Name', [particle '_Diameter_(PSN)']);
            fParticleSize = histogram(particleDiameterVsPSN(particleDiameterVsPSN(:,6)>0), 0:step:round(maxParticleD)+1, 'Normalization','probability');
            hold on
            fParticleSize = histogram(particleDiameterVsPSN(particleDiameterVsPSN(:,6)==0), 0:step:round(maxParticleD)+1, 'Normalization','probability');
            legend(['PSN around ' particle], ['No PSN around ' particle])
            xlabel([particle ' Diameter [\mum]'])
            ylabel('Relative Frequency')
            hold off
            
            %% Absolute number of Particles which lead to PSN: 
            obj.writeReport(['Total number of ' particle], length(particleDiameterVsPSN(:,1)))
            obj.writeReport(['Number of ' particle ' which lead to PSN'], length(particleDiameterVsPSN(particleDiameterVsPSN(:,6)>0)))
            
            %% Histogram: Size of PSN vs Rx Matrix Grains
            % Plot in Histogram
            step = 2;      
             
            maxRxD = max(abs(obj.RxGrains(matrix).equivalentRadius))*2;
            fRxSize = figure('Name', [matrix '_Diameter_(PSNvsDRX)']);
            fRxSize = histogram(abs(PSNGrains(matrix).equivalentRadius)*2, 0:step:round(maxRxD)+1, 'Normalization','probability');
            hold on
            fRxSize = histogram(abs(DRXGrains.equivalentRadius)*2, 0:step:round(maxRxD)+1, 'Normalization','probability');
            %fRxSize = histogram(abs(obj.RxGrains(matrix).equivalentRadius)*2, 0:step:round(maxRxD)+1, 'Normalization','probability')
            legend('PSN', 'DRX grains')
            xlabel(['Rx ' matrix ' Diameter [\mum]'])
            ylabel('Relative Frequency')
            hold off
                      
            
            %% Histogram: GOS 
            step = 0.15;
            
            fFreqGOS = figure('Name', [matrix '_Rx_GOS_Distribution']);
            fFreqGOS = histogram(PSNGrains(matrix).prop.GOS./degree, 0:step:obj.RxThreshold, 'Normalization','probability');
            hold on
            fFreqGOS = histogram(DRXGrains(matrix).prop.GOS./degree, 0:step:obj.RxThreshold, 'Normalization','probability');
            legend('PSN', 'DRXGrains')
            xlabel('GOS [°]')
            ylabel('Relative Frequency')
            hold off
            
        end
        
        function selectedGrains = neighbouringGrains(obj, particle, matrix);
            %% Find matrix grains that border particle
            %
            % INPUT
            %   particle: name of particle
            %   matrix: name of matrix
            %
            % OUTPUT
            %   selctedGrains: matrix grains that border particle
            
            %% Check if Twins have been merged
            if (isempty(obj.MergedGrains) == 1)
                grains = obj.Grains;
            else
                grains = obj.MergedGrains;
            end       
                
            %% Select Matrix grains which borders particle
            AOboundary = grains.boundary(particle,matrix);
            selectedGrains = grains(unique(AOboundary.grainId));
            
        end
     
        function lineProfile(obj)
            %% Draw misorientation line profile
            % Get Corrdinates from Cursor Data (rightclick -> cursor Data)
            
            %f = obj.plot('ipf', 3,0,0)
            f = gcf;
            hold on
            [x,y] = ginput(2);
            startCoordinates = [x(1) y(1)];
            endCoordinates = [x(2) y(2)];
            
            % If cursor data used from here on:
            lineSec = [startCoordinates; endCoordinates];
            name = ['x-y_' num2str(round(x(1))) '-' num2str(round(y(1))) '_' num2str(round(x(2))) '-' num2str(round(y(2))) ];
            
            %% Length of Segment
            vector = [(endCoordinates(1)-startCoordinates(1)); (endCoordinates(2)-startCoordinates(2))];
            lineLength = sqrt(vector(1)^2 + vector(2)^2);
            unitVector = vector/lineLength;
            
            obj.LineProfileCount = obj.LineProfileCount+1;
            %%
            
            %f = obj.plot('ipf', 3,0,0)
            %f.Name = ['LineProfile_ ' obj.LineProfileCount];
            hold on;
            line(lineSec(:,1),lineSec(:,2),'linewidth',2, 'color', 'white') %make to arrow! thickness?
            text(x(1),y(1), num2str(obj.LineProfileCount), 'color', 'white')
            
            obj.LineProfileCount = obj.LineProfileCount+1;
            text(x(2),y(2), num2str(obj.LineProfileCount), 'color', 'white')
            hold off;
            
            % Treats every phase as Ni-superalloy! Just okay with cubic
            % phases!
            ebsd = obj.Ebsd;
            ebsd.phase = 1;
            
            ebsd_line = spatialProfile(ebsd, lineSec);
                              
            f = figure('Name', ['LineProfile_ ' name ]);
           
            % Distance as x axis
            step = 1/(length(ebsd_line.y)-1)*lineLength;
            distance = transpose([0.0000:step:lineLength]);
            
            %plot(ebsd_line('Jadeite').y,angle(ebsd_line.orientations,ebsd_line('Jadeite').orientations)/degree) 
            
            % misorientation angle to the first orientation on the line
            plot(distance,...
              angle(ebsd_line(1).orientations,ebsd_line.orientations)/degree)

            % misorientation gradient
              hold on
              plot(0.5*(distance(1:end-1)+distance(2:end)),...
                angle(ebsd_line(1:end-1).orientations,ebsd_line(2:end).orientations)/degree)
            
              text(distance(1),0, num2str(obj.LineProfileCount-1), 'color', 'black')
              text(distance(end),0, num2str(obj.LineProfileCount), 'color', 'black')
              hold off

             xlabel('Distance [\mum]'); ylabel('Misorientation [°]')

             legend('Point to origin','Point to point');
             
                          
%             %% Check for GB segment misorientation (Function Test)
%             grains_line = obj.MergedGrains(unique(ebsd_line.grainId));
%             boundary = grains_line.boundary('Ni-Superlaloy','TiC');
%             misorientation = boundary.misorientation.angle ./degree

        end
        
        function traceSlipPlane(obj)
            %% Draw slip plane traces and calculate Schmid factors
            
            f = gcf;
            hold on
            [x,y] = ginput(1);
           
            ebsd = obj.EbsdOriginalGrains(x,y);
            grain = obj.Grains(ebsd.grainId);
            
            cs = ebsd.CS;
            linestyles = {'-' '--' ':' '-.'};
            %colors = {'darkblue' 'fuchsia' 'purple' 'coral' };
            
            h=symmetrise(Miller(1,1,1,cs),'unique','antipodal');
            n111 = grain.meanOrientation*h;
            nn111.antipodal=1;
            
            %% Calculate Schmid Factors
            
            sS = slipSystem.fcc(cs);
              sS = sS.symmetrise('antipodal');
              sSLocal = grain.meanOrientation * sS;
              
              % compute Schmid factor
                sigma = stressTensor.uniaxial(vector3d.Y);
                SF = sSLocal.SchmidFactor(sigma);
                
                
%                 % take the maxium allong the rows
                 [SFMax,active] = max(SF,[],2);

                  sSactive = grain.meanOrientation .* sS;
                  sStrace = sSactive.trace;
                  
%                   sSmax = grain.meanOrientation .* sS(active);
% %                 hold on
% %                 % visualize the trace of the slip plane
% %                 quiver(grain,sSactive.trace,'color','g')
%                 % and the slip direction
%                 hold on
%                 quiver(grain,sSmax.b,'color','r')
                
               
            
            %% Plot Data
            
            for i = 1:length(h)
                hold on
                v = cross(n111(:,i),-zvector);  
                
                %quiver(grain,cross(n111(:,i),-zvector),'linestyle',linestyles{i},'color',colors{i},'legend', 'off')
                factors = '';
                text = '';
                sum = 0;
                    
                for j = 1:length(sStrace)  
                    %dot(sStrace(j), normalize(v))
                    if (round(dot(sStrace(j), normalize(v)),5) == double(1)) %round in case there are small deviations
                        sf = round(abs(SF(j)), 2);
                        
                        factors = [factors  num2str(sf) '+'];
                        sum = sum + sf;
                    end
                end
                
                text = [factors '=' num2str(sum)];
                
                quiver(grain, v,'linestyle',linestyles{i},'color','darkblue','DisplayName',text)
            end
            hold off
              
        end
      
        function plotCrystal(obj)
            %% Plot crystal orientation at curser location
            
            f = gcf;
            hold on
            [x,y] = ginput(1);
            
            ebsd = obj.Ebsd(x,y);
            cS = crystalShape.cube(ebsd.CS);
            scaling = 5; % scale the crystal shape to have a nice size

            plot(x, y, -50, ebsd.orientations * cS * scaling);

        end
        
        function poleFigure(obj, n, region ,h, k, l)
            %% Plot pole figure
            %
            % INPUT
            %   n: number of overlayed orientations
            %   region: 'pt' for curser click position
            %   h,k,l: hkl of pole figure
            %
            % OUTPUT
            %   figure
            %   misorientation: misorientation matrix
            %   minAngle: smallest hkl pole deviation matrix
            
            f = gcf;
            %COLORS = {'darkblue' 'fuchsia' 'purple' 'coral'}; %more?
            COLORS = [];
            EBSD = [];
            ORI = [];
            CS = [];
            
            for i=1:n
                hold on
                if (strcmp(region, 'area')) 
                    [ebsd, rec] = selectInteractive(obj.Ebsd);
                    plot(rec, 'DisplayName', 'selected')
                elseif(strcmp(region, 'poly'))
                    poly = selectPolygon
                    ebsd = ebsd(inpolygon(ebsd,poly))
                else
                    [x,y] = ginput(1); 
                    ebsd = obj.Ebsd(x,y);
                    
                    cS = crystalShape.cube(ebsd.CS);
                    scaling = 5; % scale the crystal shape to have a nice size
                    
                    ipfKey = ipfHSVKey(ebsd);
                    COLORS(i,:) = ipfKey.orientation2color(ebsd.orientations);
                    
                    plot(x, y, -50, ebsd.orientations * cS * scaling, 'FaceColor', COLORS(i,:), 'FaceAlpha', 0.6); %change
                 end

                ori = ebsd.orientations;
                cs = ebsd.CS;
                
                EBSD = [EBSD ebsd];
                CS = [CS cs];
                ORI = [ORI ori];
            end
            
            g = figure;
            
            for i=1:n
                hold on
                plotPDF(ORI(i),  Miller(h,k,l, CS(i)), 'MarkerColor', COLORS(i,:), 'DisplayName', int2str(i), 'upper','projection','stereo', 'grid', 'grid_res', 5*degree)
                
                for j=(i+1):n
                    misorientation(i,j) = angle(ORI(i), ORI(j)) /degree;
                    
                    %https://mtex-toolbox.github.io/MisorientationTheory.html
                    mori = inv(ORI(i)) * ORI(j); %2 into 1 reference frame;
                    a = angle(mori.symmetrise * Miller(h,k,l, CS(2)), Miller(h,k,l, CS(1))) /degree;
                    MinAngle(i,j) = min(a);
                end
            end
            g.Name = ['PoleFigure_ ' int2str(h) '' int2str(k) '' int2str(l)];
            legend({},'location','NorthWest','FontSize',13);
            
            try
                misorientation
                MinAngle
            end
%             %angle(ORI(1)*Miller(h,k,l, CS(1)), ORI(2)*Miller(h,k,l, CS(1))) /degree
%             mori = inv(ORI(1)) * ORI(2); %2 into 1 reference frame;
%              
%             a = angle(mori.symmetrise * Miller(h,k,l, CS(2)), Miller(h,k,l, CS(1))) /degree;
%             minAngle = min(a)

        end
                    
        function calculateKAM(obj)
            %% Calculate kernel average misorientation (KAM)
            % Order 1 with 2.5 degrees 
            % Also calculates grain average misorientation (GAM)
            %
            % OUTPUT
            %   obj.KAM, obj.Ebsd.prop.KAM
            %   obj.Grains.prop.GAM, obj.MergedGrains.prop.GAM 
            
                ebsd = obj.Ebsd;
                kam = ebsd.KAM('threshold',2.5*degree,'order',1);
                
                obj.KAM = kam;
                obj.Ebsd.prop.KAM = kam;
                
                if (isempty(obj.MergedGrains) == 1)
                   obj.Grains.prop.GAM = obj.Ebsd.grainMean(obj.Ebsd.prop.KAM, obj.MergedGrains);
                else
                   obj.MergedGrains.prop.GAM = obj.Ebsd.grainMean(obj.Ebsd.prop.KAM, obj.MergedGrains);
                end
                
        end
        
        function calculateGROD(obj)
            %% Calculates grain orientation deviation (GROD)
            %
            % OUTPUT
            %   obj.Ebsd.prop.GROD
            
                obj.EbsdOriginalGrains.prop.GROD = obj.EbsdOriginalGrains.calcGROD(obj.Grains);
                obj.Ebsd.prop.GROD = obj.EbsdOriginalGrains.prop.GROD;
        end
        
        function calculateGND(obj) 
                %% GND calculation based on doi.org/10.1016/j.scriptamat.2008.01.050
                %Pantleon, Resolving the geometrically necessary dislocation content by conventional electron backscattering diffraction, Scripta Materialia, 2008
                
                %try
                    totalEbsd = obj.Ebsd('indexed').gridify; 
                    totalEbsd.prop.GND = nan(length(totalEbsd),1);      

                    for i=1:(length(totalEbsd.mineralList)-1)
                        
                        ebsd = totalEbsd(totalEbsd.phase==i).gridify;
                        %gnd = [];
                        %gnd = ones(length(ebsd), 1);

                        % compute the curvature tensor
                        kappa = ebsd.curvature;

                        %Kroener: curvature is directly related to dislocation
                        %density
                        alpha = kappa.dislocationDensity;

                        %Crystallographic Dislocations
                        dS = dislocationSystem.fcc(ebsd.CS);

                        % size of the unit cell
                        a = norm(ebsd.CS.aAxis);

                        % in bcc and fcc the norm of the burgers vector is sqrt(3)/2 * a
                        [norm(dS(1).b), norm(dS(end).b), sqrt(3)/2 * a]

                        %gnd = clalcGND(ebsd, dS)

                        %% Energy of dislocations
                        % b = Burger vector, nu = Posson ratio, R, r

                        nu = 0.3; %Poisson Ratio
                        E = 1; %Energy

                        % energy of the edge dislocations
                        dS(dS.isEdge).u = E;

                        % energy of the screw dislocations
                        dS(dS.isScrew).u = E - nu;

                        % Roate from crystal reference frame to specimen reference
                        % frame

                        dSRot = ebsd.orientations * dS;

                        %% Fittind Dislocations to the incomplete dislocation dneisty tensor
                        %Fit 
                        [rho,factor] = fitDislocationSystems(kappa,dSRot);

                        % the restored dislocation density tensors
                        alpha = sum(dSRot.tensor .* rho,2);

                        % we have to set the unit manualy since it is not stored in rho
                        alpha.opt.unit = '1/um';

                        kappa = alpha.curvature;

                        %The unit of the densities h in our example is 1/um * 1/au 
                        %where 1/um comes from the unit of the curvature tensor an 1/au 
                        %from the unit of the Burgers vector. In order to transform h to SI
                        %units, i.e., 1/m^2 we have to multiply it with 10^16 (= factor). 

                        gnd{i} = factor*sum(abs(rho .* dSRot.u),2);

                        %totalGND(totalEbsd.phase==i) = gnd;
                    end

                    obj.GND = gnd;
%                 catch
%                     display('GND: memory not sufficient.');
%                 end
        end
        
        function serratedGBs(obj, phase)
            %% Colour all serrated GBs in red!
            % smoothes GBs and selectes serrates above a thereshold change
            % in segment length change
            % Other Ideas: Difference between sum(segLength) and linear
            % connection between triple points; fourrier transform;
            % curvature of grains
            % or statisics: based on perimeter/equivialent Diameter
            
            %% Check if Twins have been merged
            if (isempty(obj.MergedGrains) == 1)
                grains = obj.Grains('indexed');
            else
                grains = obj.MergedGrains('indexed');
            end
            
                        
            boundaries = grains.boundary();
            
            %histogram(boundaries.segLength);
            
            sFactor = 1000;
            threshold = 0.1;
            
            grainsSmoothed = grains.smooth(sFactor);
            
            figure
            plot(grains.boundary)
            hold on
            plot(grainsSmoothed.boundary, 'linecolor', 'red');
            
            figure
            plot(grains)
            hold on;
            serratedGbId = {};
            for (i=1:max(grains.id))
                
                g = grains(grains.id(i));
                gSmoothed = grainsSmoothed(grains.id(i));
                
                neighbours = unique(g.boundary.grainId);
                

                for (k=1:length(neighbours))
                    
                    if(grains.id(i) == neighbours(k))
                       continue 
                    end
                    
                    gbOriginal = g.boundary(g.boundary.hasGrain(neighbours(k)));
                    gbSmoothed = gSmoothed.boundary(gSmoothed.boundary.hasGrain(neighbours(k)));
                    
                    diff = (sum(gbSmoothed.segLength) - sum(gbOriginal.segLength)) / sum(gbOriginal.segLength);
                
                    if (abs(diff) > threshold)
                        serratedGbId(end+1) = {[g.id neighbours(k) diff]};
                        
                        plot(g.boundary(g.boundary.hasGrain(neighbours(k))), 'linecolor', 'red')
                                                
                        hold on;
                    end
                end
                
            end
            
            
            %% AUTOMATE (So far MANUAL!)
            %General GB Distribution            
            AOboundary = grains.boundary(phase,phase)
            Angle = AOboundary.misorientation.angle;
            
            % Plot in Histogram
            f = figure('Name', ['BoundaryCharacter_' phase '_' phase]);
            histogram(Angle./degree,0:5:65)
            hold on;
            %Between Grain 32 &33
            b = grains(32).boundary(grains(32).boundary.hasGrain(33));
            Angle = b.misorientation.angle;
            histogram(Angle./degree, 0:5:65)
            
            legend('all GBs', 'selected serrated GB');

            
            
        end

    end
end

