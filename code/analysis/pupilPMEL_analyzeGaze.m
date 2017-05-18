function [ avgPackets ] = pupilPMEL_analyzeGaze( mergedPacketCellArray, analysisDir)
% function [ avgPackets ] = pupilPMEL_analyzeGaze( mergedPacketCellArray, analysisDir)
%
%  Analysis of the horizontal (x) position of the eye during presentation
%  of the stimulus and during the ISI.
%
%  We use the x-position for two reasons: 1) it is the better tracked
%  parameter in the LiveTrack data; 2) the eye has greater variation in
%  position in the horizontal plane.
%
%  The standard deviation of gaze position is measured for each trial
%  during the stimulus and during the ensuing ISI period. The means across
%  trials, runs, and subjects are ultimately obtained, and the difference
%  in the mean std of gaze position is examined between the stimulus and
%  ISI time periods.


% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=1; % No normalization for the gaze and blink vectors

nStimTypes=6;

%% Loop through the sessions
nSessions=size(mergedPacketCellArray,2);

% pre-allocate the array
% I assume that there are not more than 12 runs for any direction/subject,
% and that there are not more than 10 of a given stimulus instance per run.
stdGazeX_StimulusAllInstances=nan(nSessions,nStimTypes,12,10);
stdGazeX_ISIAllInstances=nan(nSessions,nStimTypes,12,10);

for ss = 1:nSessions
    
    nRuns=size(mergedPacketCellArray{ss},2);
    
    for rr = 1:nRuns
        
        % grab a packet that corresponds to a run for a given subject
        theRunPacket=mergedPacketCellArray{ss}{rr};
        theRunPacket_GazeX=theRunPacket;
        theRunPacket_GazeX.response.values=theRunPacket.response.gazeX;
        
        % How many stimtypes are there?
        theUniqueStimTypes = unique(theRunPacket.stimulus.metaData.stimTypes);
        stimTypeCounter=zeros(1,length(theUniqueStimTypes));
        
        for ii = 1:size(theRunPacket.stimulus.values,1)
            
            % Put the instance index into the splitParams
            splitParams.instanceIndex = ii;
            
            % grab the packet for this particular instance
            singlePacket=splitOffAnInstancePacket(theRunPacket_GazeX,splitParams);
            
            % identify the stim type for this instance
            thisStimType = theRunPacket.stimulus.metaData.stimTypes(ii);
            uniqueStimIdx = find(theUniqueStimTypes==thisStimType);
            
            % Obtain a measure of the std of gaze position in the time
            % window of the stimulus and 6 seconds later.
            stimTypeCounter(uniqueStimIdx) = stimTypeCounter(uniqueStimIdx)+1;
            stdGazeX_StimulusAllInstances(ss, uniqueStimIdx, rr, stimTypeCounter(uniqueStimIdx))= ...
                std(singlePacket.response.values(1:3000));
            stdGazeX_ISIAllInstances(ss, uniqueStimIdx, rr, stimTypeCounter(uniqueStimIdx))= ...
                std(singlePacket.response.values(3001:13000));
        end % loop over instances
    end % loop over runs
end % loop over sessions

stdGazeX_Stimulus_means = nanmean(nanmean(stdGazeX_StimulusAllInstances,4),3);
stdGazeX_ISI_means = nanmean(nanmean(stdGazeX_ISIAllInstances,4),3);
stdGazeMeasure = stdGazeX_ISI_means - stdGazeX_Stimulus_means;


%% Plot CRFs of the mean gaze effect
contrastLabels={'25','50','100','200','400'};

experimentSet{1}=[1 2 3 4];
experimentSet{2}=[5 6 7 8];
experimentNames={'MaxLMSCRF','MaxMelCRF'};
figHandle = figure;
set(gcf, 'PaperSize', [8.5 11]);
subPlotHandle=subplot(1,1,1);

for ee = 1:length(experimentNames)
    
    switch experimentNames{ee}
        case 'MaxLMSCRF'
            lineColorBase=[.25 .25 .25];
        case 'MaxMelCRF'
            lineColorBase=[0 0 1];
    end
    meanGazeEffect = squeeze(mean(stdGazeMeasure(experimentSet{ee},:),1));
    semGazeEffect = squeeze(std(stdGazeMeasure(experimentSet{ee},:),1))/sqrt(length(experimentSet{ee}));
    
    plotTitle=['Subject mean gaze X-position effect (std ISI vs. stim) ±SEM'];
    % plot the response and SEM across subjects
    pupilPMEL_PlotCRF( subPlotHandle, 1:5, meanGazeEffect(1:5), semGazeEffect(1:5), ...
        'xTickLabels',contrastLabels,...
        'xlim',[0 6],...
        'ylim',[0 0.1],...
        'lineColor',lineColorBase,...
        'markerColor',lineColorBase,...
        'errorColor',lineColorBase,...
        'plotTitle',plotTitle,...
        'xLabel','contrast',...
        'yLabel','gaze effect [diff std]');
    hold on
end % loop over experiments
hold off
% Save the plots
plotFileName=fullfile(analysisDir, 'figures', ['GazeEffect.pdf']);
pupilPMEL_suptitle(figHandle,['Gaze X position during ISI vs. stimulus']);
set(figHandle,'Renderer','painters');
print(figHandle, plotFileName, '-dpdf', '-fillpage');
close(figHandle);

end % function