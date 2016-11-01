function [ theResponseCellArray ] = pupilPMEL_makeCarryOverMatrix(mergedPacketCellArray, resultsVariable, dropboxAnalysisDir)

subAnalysisDirectory='carryOver';


% have to do something a little bit different to handle the
% effectiveContrast variable
if strcmp(inputname(2), 'effectiveContrast_sessxrunxevent');
    effectiveContrast_sessxrunxevent = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/pupilMelanopsinMRIAnalysis/fitFCONModelToIndividualResponses/effectiveContrast_sessxrunxevent.mat');
    resultsVariable = effectiveContrast_sessxrunxevent.myResultsVariable;
end

nSubjects=size(resultsVariable,1);
%% Assemble the carry-over matrix

stimTypes = 1:6; % stimuli are 1 for 25% up to 5 for 400%, 6 is attention

% Pre-allocate results variable
for ss = 1:nSubjects;
    for priorStimulus = 1:6;
        for currentStimulus = 1:6;
            theResponseCellArray{ss}{priorStimulus,currentStimulus} = [];
        end
    end
end
    

for ss = 1:nSubjects; % loop over subjects
    nRuns = sum(resultsVariable(ss,:,1)~=0); % determine how many runs for that subject
    for rr = 1:nRuns; % loop over runs
        stimulusList = [];
        nInstances = length(resultsVariable(ss,rr,:));
        for ii = 1:nInstances;
            currentStimulus = mergedPacketCellArray{1,ss}{rr}.stimulus.metaData.stimTypes(ii);
            stimulusList = [stimulusList, currentStimulus];
            if ii == 1;
            else
                priorStimulus = stimulusList(ii-1);
                theResponseCellArray{ss}{priorStimulus,currentStimulus}= ...
                    [theResponseCellArray{ss}{priorStimulus,currentStimulus}, ...
                    resultsVariable(ss,rr,ii)];
            end
        end
    end
end

% Create mean response value for a given current stimulus per subject
for ss = 1:nSubjects;
    meanResultsVariable{ss} = zeros(5,5);
end
   

for ss = 1:nSubjects;
    for previousStimulus = 1:5;
        for currentStimulus = 1:5;
            meanResponseCellArray{ss}(previousStimulus,currentStimulus) = mean(theResponseCellArray{ss}{previousStimulus,currentStimulus});
        end
    end
end


%% Do some plotting

% Get some information about what we're plotting
subjectKey = [];
projectKey = [];
for ss = 1:size(mergedPacketCellArray,2);
    subjectKey{ss} = mergedPacketCellArray{ss}{1}.metaData.subjectName;
    projectKey{ss} = mergedPacketCellArray{ss}{1}.metaData.projectName;
end
colorList = ['r', 'b', 'g', 'y', 'm'];
contrastList = {'25%', '50%', '100%', '200%', '400%'};
if strcmp(inputname(2),'betaPercentCarryOver');
    yAxis = 'Percent Change (%)';
    save = 'percentChange';
elseif strcmp(inputname(2), 'betaRawCarryOvery');
    yAxis = 'Raw Pupil Diameter (mm)';
    save = 'raw';
elseif strcmp(inputname(2), 'baselineSizeCarryOver');
    yAxis = 'Baseline Pupil Diameter (mm)';
    save = 'baselineSize';
elseif strcmp(inputname(2), 'effectiveContrast_sessxrunxevent');
    yAxis = 'Effective Contrast (%)';
    save = 'effectiveContrast';
end


% label whether a session is LMS or Mel
subjectsLMS = [];
subjectsMel = [];
for x = 1:length(projectKey);
    if strcmp(projectKey(x), 'MelanopsinMRMaxLMSCRF');
        subjectsLMS = [subjectsLMS x];
    end
    if strcmp(projectKey(x), 'MelanopsinMRMaxMelCRF');
        subjectsMel = [subjectsMel x];
    end
end

% store the type of data for file saving purposes
if isfield(mergedPacketCellArray{1}{1}.metaData, 'simulationStyle');
    if strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'arousal');
        simFileName = '_simulatedArousal';
    elseif strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'carryOver')
        simFileName = '_simulatedCarryOver';
    
    end
else
        simFileName = '';
end

% Plot the LMS subjects
plotFig = figure;
set(gca,'Xtick',[1 2 3 4 5],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
xlabel('Previous Stimuli')
xlim([0 6])
ylabel(yAxis) 

for s = subjectsLMS;
    hold on
    if length(subjectsLMS)>1
        subplot(2, length(subjectsLMS)/2, s)
        set(gca,'Xtick',[1 2 3 4 5],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
        xlabel('Previous Stimuli')
        xlim([0 6])
        ylabel(yAxis)
        hold on
        
    end
    for priorStimulus = 1:5;
        for currentStimulus = 1:5;
            for ii = 1:length(theResponseCellArray{s}{priorStimulus,currentStimulus});
                x = priorStimulus;
                y = theResponseCellArray{1,s}{priorStimulus,currentStimulus}(ii);
                plot(x,y,'o', 'Color', colorList(currentStimulus));
            end
        end
    % plot mean values
    plot(1:5,meanResponseCellArray{1,s}(priorStimulus,1:5), 'Color', colorList(priorStimulus))
    end
    
    
end

outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, ['LMS_carryOver_', save, simFileName, '.png']), 'png');
close(plotFig);

% Plot the Melanopsin subjects
plotFig = figure;
set(gca,'Xtick',[1 2 3 4 5],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
xlabel('Previous Stimuli')
xlim([0 6])
ylabel(yAxis) 

for s = subjectsMel;
    hold on
    if length(subjectsMel)>1
        subplot(2, length(subjectsMel)/2, find(subjectsMel==s))
        set(gca,'Xtick',[1 2 3 4 5],'XTickLabel',{'25%', '50%', '100%', '200%', '400%'})
        xlabel('Previous Stimuli')
        xlim([0 6])
        ylabel(yAxis)
        hold on
        
    end
    for priorStimulus = 1:5;
        for currentStimulus = 1:5;
            for ii = 1:length(theResponseCellArray{s}{priorStimulus,currentStimulus});
                x = priorStimulus;
                y = theResponseCellArray{1,s}{priorStimulus,currentStimulus}(ii);
                plot(x,y,'o', 'Color', colorList(currentStimulus));
            end
        end
    % plot mean values
    plot(1:5,meanResponseCellArray{1,s}(priorStimulus,1:5), 'Color', colorList(priorStimulus))
    end
    
    
end

outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, ['Mel_carryOver_', save, simFileName, '.png']), 'png');
close(plotFig);
