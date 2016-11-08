function [ rsqLag ] = plotLag(mergedPacketCellArray, betaRawCarryOver, baselineSizeCarryOver, dropboxAnalysisDir)

% variables for plotting purposes
subAnalysisDirectory='lagPlots';
if isfield(mergedPacketCellArray{1}{1}.metaData, 'simulationStyle');
    if strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'arousal');
        simFileName = '_simulatedArousal';
    elseif strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'carryOver')
        simFileName = '_simulatedCarryOver';
    elseif strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'attention')
        simFileName = '_simulatedAttention';
        
    end
else
    simFileName = '';
end
subjectKey = [];
projectKey = [];
for ss = 1:size(mergedPacketCellArray,2);
    subjectKey{ss} = mergedPacketCellArray{ss}{1}.metaData.subjectName;
    projectKey{ss} = mergedPacketCellArray{ss}{1}.metaData.projectName;
end


% set how far we want to be able to slide the window
lagRange = -7000:200:15000;


% Pre-allocate for results storage
for ss = 1:size(mergedPacketCellArray,2);
    rsqLag{ss} = [];
    
end

plotFig = figure;
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
   ss
    for lagIndex = lagRange;
        
        for stimuliTypes = 1:5;
            baselineSize{stimuliTypes} = [];
            beta{stimuliTypes} = [];
        end
        for rr = 1:(size(mergedPacketCellArray{ss},2));
            for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1)
                % use splitOffInstancePacket to figure out stimOnset time
                % define the split params
                splitParams.instanceIndex=[]; % will hold the instance index
                splitParams.splitDurationMsecs=13000; % Grab 13 second windows
                splitParams.normalizationWindowMsecs=100; % define the size of the norm window
                splitParams.instanceIndex = ii;
                % grab the packet
                singlePacket=splitOffAnInstancePacket(mergedPacketCellArray{ss}{rr},splitParams);
                % grab stimOnset time
                stimOnset = singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs;

                if (stimOnset-lagIndex) > 0 && (stimOnset-lagIndex) < 448000;
                    for stimuliTypes = 1:5;
                        contrast = mergedPacketCellArray{ss}{rr}.stimulus.metaData.stimTypes(ii);
                        if contrast == stimuliTypes;
                            
                            baselineSize{stimuliTypes} = [baselineSize{stimuliTypes}, mergedPacketCellArray{ss}{rr}.response.values(stimOnset-lagIndex)];
                            beta{stimuliTypes} = [beta{stimuliTypes}, betaRawCarryOver(ss,rr,ii)];
                        end
                    end
                    
                end
            end
        end
        for stimuliTypes = 1:5;
            baselineSize{stimuliTypes} = [ones(1,length(baselineSize{stimuliTypes})); baselineSize{stimuliTypes}];
            baselineSize{stimuliTypes} = baselineSize{stimuliTypes}';
            beta{stimuliTypes} = beta{stimuliTypes}';
            [B,BINT,R,RINT,STATS] = regress(beta{stimuliTypes},baselineSize{stimuliTypes});
            rsq(stimuliTypes) = STATS(1);
            rsqLag{ss}(1,(find(lagRange==lagIndex))) = lagIndex;
            rsqLag{ss}(2,(find(lagRange==lagIndex))) = mean(rsq);
        end
    end
    
    if strcmp(projectKey{ss},'MelanopsinMRMaxLMSCRF');
        color = 'b';
        
    end
    if strcmp(projectKey{ss},'MelanopsinMRMaxMelCRF');
        color = 'r';
        
    end
    
    hold on
    plot(rsqLag{ss}(1,:),rsqLag{ss}(2,:), color);
    ylabel('R2');
    xlabel('Trial Lag');
    
end





outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, ['lagPlot', simFileName, '.png']), 'png');
close(plotFig);


end
