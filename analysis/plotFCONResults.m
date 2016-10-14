function [theResult] = plotFCONResults(mergedPacketCellArray)

effectiveContrast_sessxrunxevent = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/pupilMelanopsinMRIAnalysis/fitFCONModelToIndividualResponses/effectiveContrast_sessxrunxevent.mat');
myResultsVariable = effectiveContrast_sessxrunxevent.myResultsVariable;

for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
     
    for contrast = 1:6;
        fcon{ss,contrast}=[];
    end
end

for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
     
    for rr = 1:size(mergedPacketCellArray{ss},2); % loop over runs
        for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1); % loop over instances
            contrast = mergedPacketCellArray{ss}{rr}.stimulus.metaData.stimTypes(ii);
            fcon{ss,contrast} = [fcon{ss,contrast} myResultsVariable(ss,rr,ii)];
        end
    end
end

fconMean =zeros(8,6)
fconSEM = zeros(8,6)
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
    
    for contrast = 1:6;
        fconMean(ss,contrast) = mean(fcon{ss,contrast})
        fconSEM(ss,contrast) = std(fcon{ss,contrast})/sqrt(length(fcon{ss,contrast}))
    end
end

    
                
plot(fconMean)
    





% do each session separately
% plot mean and SEM across contrast level