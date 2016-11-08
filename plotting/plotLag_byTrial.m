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


% set how far we want to be able to slide the window
lagRange = -10:10;

% Pre-allocate for results storage
for ss = 1:size(mergedPacketCellArray,2);
    rsqLag{ss} = [];
    
end

plotFig = figure;
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
    for lagIndex = lagRange;
        for stimuliTypes = 1:5;
            baselineSize{stimuliTypes} = [];
            beta{stimuliTypes} = [];
        end
        for rr = 1:(size(mergedPacketCellArray{ss},2));
            for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1)
                if (ii-lagIndex) > 0 && (ii-lagIndex) < (size(mergedPacketCellArray{ss}{rr}.stimulus.values,1))
                    for stimuliTypes = 1:5;
                        contrast = mergedPacketCellArray{ss}{rr}.stimulus.metaData.stimTypes(ii);
                        if contrast == stimuliTypes;
                            
                            baselineSize{stimuliTypes} = [baselineSize{stimuliTypes}, baselineSizeCarryOver(ss,rr,(ii-lagIndex))];
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
            rsqLag{ss}(1,(lagIndex+max(lagRange)+1)) = lagIndex;
            rsqLag{ss}(2,(lagIndex+max(lagRange)+1)) = mean(rsq);
        end
    end
    
    
    hold on
    plot(rsqLag{ss}(1,:),rsqLag{ss}(2,:));
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
