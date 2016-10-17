function [theResult] = plotFCONResults(mergedPacketCellArray, dropboxAnalysisDir)


subAnalysisDirectory='fitFCONModelToIndividualResponses';

effectiveContrast_sessxrunxevent = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/pupilMelanopsinMRIAnalysis/fitFCONModelToIndividualResponses/effectiveContrast_sessxrunxevent.mat');
myResultsVariable = effectiveContrast_sessxrunxevent.myResultsVariable;
contrastList = ['25%' '50%' '100%' '200%' '400%'];
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};


for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
     
    for contrast = 1:6;
        fcon{ss,contrast}=[];
        fcon_log{ss,contrast}=[];
    end
end

for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
     
    for rr = 1:size(mergedPacketCellArray{ss},2); % loop over runs
        for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1); % loop over instances
            contrast = mergedPacketCellArray{ss}{rr}.stimulus.metaData.stimTypes(ii);
            fcon_log{ss,contrast} = [fcon_log{ss,contrast} myResultsVariable(ss,rr,ii)];
            fcon{ss,contrast} = [fcon{ss,contrast} 10^(myResultsVariable(ss,rr,ii))*100];
        end
    end
end



fconMean =zeros(8,6);
fconSEM = zeros(8,6);
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
    
    for contrast = 1:6;
        fconMean(ss,contrast) = mean(fcon{ss,contrast});
        fconSEM(ss,contrast) = std(fcon{ss,contrast})/sqrt(length(fcon{ss,contrast}));
    end
end

fconMean
    
plotFig = figure;
hold on
set(gca,'Xtick',[25 50 100 200 400],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
ylim([0 600])
xlim([0 450])
xlabel('Presented Contrast')
ylabel('Effective Contrast')   
unity=plot([0:600],[0:600], ':');



for ss = 1:size(mergedPacketCellArray,2);
    
    % set color to differentiate stimulus type; LMS in blue, melanopsin in
    % red
    if strcmp(projectKey{ss},'LMS');
        color = 'b';
        
    end
    if strcmp(projectKey{ss},'Melanopsin');
        color = 'r';
        
    end

    % perform the plotting
    plot([25 50 100 200 400],fconMean(ss,1:5), 'o', 'Color', color)
  
    
    
    % add SEM
    h=errorbar([25 50 100 200 400],fconMean(ss,1:5),fconSEM(ss,1:5), color); set(h,'linestyle','-')
    
    % save representative from each stimulus type for legend
    if ss == 4;
        lms=errorbar([25 50 100 200 400],fconMean(ss,1:5),fconSEM(ss,1:5), color); set(h,'linestyle','-');
    end
    if ss == 8;
        mel=errorbar([25 50 100 200 400],fconMean(ss,1:5),fconSEM(ss,1:5), color); set(h,'linestyle','-');
    end
    
    
end

legend([unity lms mel], 'Unity', 'LMS', 'Melanopsin', 'Location', 'northwest');
outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, 'effectiveContrastByPresentedContrast.png'), 'png');
close(plotFig);
 

end