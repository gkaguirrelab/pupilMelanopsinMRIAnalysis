function [theResult] = plotFCONResults(mergedPacketCellArray, dropboxAnalysisDir, myResultsVariable)


subAnalysisDirectory='fitFCONModelToIndividualResponses';

%effectiveContrast_sessxrunxevent = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/pupilMelanopsinMRIAnalysis/fitFCONModelToIndividualResponses/effectiveContrast_sessxrunxevent.mat');
%myResultsVariable = effectiveContrast_sessxrunxevent.myResultsVariable;

% describe subjects for labeling purposes
contrastList = ['25%' '50%' '100%' '200%' '400%'];
subjectKey = {'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1';	'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1'};
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};
subjectsLMS = [];
subjectsMel = [];
for x = 1:length(projectKey);
    if strcmp(projectKey(x), 'LMS');
        subjectsLMS = [subjectsLMS x];
    end
    if strcmp(projectKey(x), 'Melanopsin');
        subjectsMel = [subjectsMel x];
    end
end


%% Average effective contrast for each contrast level for all subjects

% Make room for the combined dataset
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects    
    for contrast = 1:6;
        fcon{ss,contrast}=[];
        fcon_log{ss,contrast}=[];
    end
end

% Actually make the combined dataset. Result will be an 8x6 (subject x 
% contrast) cell array. Each cell will contain the effective contrast
% values for all stimulus shown of a given contrast for a given subject
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects  
    for rr = 1:size(mergedPacketCellArray{ss},2); % loop over runs
        for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1); % loop over instances
            contrast = mergedPacketCellArray{ss}{rr}.stimulus.metaData.stimTypes(ii);
            fcon_log{ss,contrast} = [fcon_log{ss,contrast} myResultsVariable(ss,rr,ii)];
            fcon{ss,contrast} = [fcon{ss,contrast} 10^(myResultsVariable(ss,rr,ii))*100];
        end
    end
end


% make average per subject for a given contrast level, also calculate the
% SEM while we're at it
fconMean = zeros(8,6);
fconSEM = zeros(8,6);
for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
    
    for contrast = 1:6;
        fconMean(ss,contrast) = mean(fcon{ss,contrast});
        fconSEM(ss,contrast) = std(fcon{ss,contrast})/sqrt(length(fcon{ss,contrast}));
    end
end






% plot with individual subjects
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



 

%% Make average across subjects for each stimulus type

% preallocate space for combined across subjects data set. Note stimulus =
% 1 corresponds to LMS, stimulus = 2 corresponds to melanopsin
for contrast = 1:6;
    for stimulus = 1:2;
        fconAcrossSubjects{stimulus, contrast} = [];
        fconAcrossSubjects_log{stimulus, contrast} = [];
    end
end

% make the combined data set. Note stimulus =
% 1 corresponds to LMS, stimulus = 2 corresponds to melanopsin
for contrast = 1:6;
    for ss = 1:size(mergedPacketCellArray,2); % loop over subjects
        if strcmp(projectKey{ss}, 'LMS');
            fconAcrossSubjects{1, contrast} = [fconAcrossSubjects{1, contrast} fcon{ss,contrast}];
            fconAcrossSubjects_log{1, contrast} = [fconAcrossSubjects_log{1, contrast} fcon_log{ss,contrast}];
        end
        if strcmp(projectKey{ss}, 'Melanopsin');
            fconAcrossSubjects{2, contrast} = [fconAcrossSubjects{2, contrast} fcon{ss,contrast}];
            fconAcrossSubjects_log{2, contrast} = [fconAcrossSubjects_log{2, contrast} fcon_log{ss,contrast}];
        end
    end
end

fconAcrossSubjectsMean = zeros(2,6);
fconAcrossSubjectsSEM = zeros(2,6);
for stimulus = 1:2 % loop over stimulus types
    
    for contrast = 1:6;
        fconAcrossSubjectsMean(stimulus,contrast) = mean(fconAcrossSubjects{stimulus,contrast});
        fconAcrossSubjectsSEM(stimulus,contrast) = std(fconAcrossSubjects{stimulus,contrast})/sqrt(length(fconAcrossSubjects{stimulus,contrast}));
    end
end



    
plotFig = figure;
title('Mean')
hold on
set(gca,'Xtick',[25 50 100 200 400],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
ylim([0 600])
xlim([0 450])
xlabel('Presented Contrast')
ylabel('Effective Contrast')   
unity=plot([0:600],[0:600], ':');



for stimulus = 1:2
    
    % set color to differentiate stimulus type; LMS in blue, melanopsin in
    % red
    if stimulus == 1;
        color = 'b';
        
    end
    if stimulus == 2;
        color = 'r';
        
    end

    % perform the plotting
    plot([25 50 100 200 400],fconAcrossSubjectsMean(stimulus,1:5), 'o', 'Color', color)
  
    
    
    % add SEM
    h=errorbar([25 50 100 200 400],fconAcrossSubjectsMean(stimulus,1:5),fconAcrossSubjectsSEM(stimulus,1:5), color); set(h,'linestyle','-')
    
    % save representative from each stimulus type for legend
    if stimulus == 1;
        lms=errorbar([25 50 100 200 400],fconAcrossSubjectsMean(stimulus,1:5),fconAcrossSubjectsSEM(stimulus,1:5), color); set(h,'linestyle','-');
    end
    if stimulus == 2;
        mel=errorbar([25 50 100 200 400],fconAcrossSubjectsMean(stimulus,1:5),fconAcrossSubjectsSEM(stimulus,1:5), color); set(h,'linestyle','-');
    end
    
    
end

legend([unity lms mel], 'Unity', 'LMS', 'Melanopsin', 'Location', 'northwest');
outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, 'effectiveContrastByPresentedContrast_acrossSubjects.png'), 'png');
close(plotFig);

% plot by median, not average
fconAcrossSubjectsMedian = zeros(2,6);
fconAcrossSubjectsSEM = zeros(2,6);
for stimulus = 1:2 % loop over stimulus types
    
    for contrast = 1:6;
        fconAcrossSubjectsMedian(stimulus,contrast) = median(fconAcrossSubjects{stimulus,contrast});
        fconAcrossSubjectsSEM(stimulus,contrast) = std(fconAcrossSubjects{stimulus,contrast})/sqrt(length(fconAcrossSubjects{stimulus,contrast}));
    end
end



    
plotFig = figure;
hold on
title('Median')
set(gca,'Xtick',[25 50 100 200 400],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
ylim([0 600])
xlim([0 450])
xlabel('Presented Contrast')
ylabel('Effective Contrast')   
unity=plot([0:600],[0:600], ':');



for stimulus = 1:2
    
    % set color to differentiate stimulus type; LMS in blue, melanopsin in
    % red
    if stimulus == 1;
        color = 'b';
        
    end
    if stimulus == 2;
        color = 'r';
        
    end

    % perform the plotting
    plot([25 50 100 200 400],fconAcrossSubjectsMedian(stimulus,1:5), 'o', 'Color', color)
  
    
    
    % add SEM
    h=errorbar([25 50 100 200 400],fconAcrossSubjectsMedian(stimulus,1:5),fconAcrossSubjectsSEM(stimulus,1:5), color); set(h,'linestyle','-')
    
    % save representative from each stimulus type for legend
    if stimulus == 1;
        lms=errorbar([25 50 100 200 400],fconAcrossSubjectsMedian(stimulus,1:5),fconAcrossSubjectsSEM(stimulus,1:5), color); set(h,'linestyle','-');
    end
    if stimulus == 2;
        mel=errorbar([25 50 100 200 400],fconAcrossSubjectsMedian(stimulus,1:5),fconAcrossSubjectsSEM(stimulus,1:5), color); set(h,'linestyle','-');
    end
    
    
end

legend([unity lms mel], 'Unity', 'LMS', 'Melanopsin', 'Location', 'northwest');
outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, 'effectiveContrastByPresentedContrast_acrossSubjects_median.png'), 'png');
close(plotFig);

end