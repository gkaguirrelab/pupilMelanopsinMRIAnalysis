function [ ] = pupilPMEL_plotAveragePupilResponses( mergedPacketCellArray, twoComponentFitToData, dropboxAnalysisDir )
% Create plots of average evoked pupil responses within and across subjects

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, 'lowFreqClean', true, 'aggregateMethod', 'mean');

%% Plot the within-subject data
NSessionsMerged=size(mergedPacketCellArray,2);
NStimTypes=5;

% Assign a name for this sub-analysis directory for saving plots and fits
subAnalysisDirectory='fitTPUPModelToAverageResponse';


for ee = 1:2 % loop over experiments, LMS and Mel
    figHandle = figure;
    set(gcf, 'PaperSize', [8.5 11]);
    for ss = 1:NSessionsMerged/2 % loop over subjects
        
        thisPacketIndex=(ee-1)*(NSessionsMerged/2)+ss;
        % Identify the stimulus type
        stimType=mergedPacketCellArray{6}{1,1}.stimulus.metaData.stimLabels;
        stimType=stimType{1};
        lineColorBase=[.25 .25 .25];
        switch stimType
            case 'MelanopsinMRMaxLMS_PulseMaxLMS_3s_CRF16sSegment_025Pct'
                lineColorBase=[.25 .25 .25];
            case 'MelanopsinMRMaxMel_PulseMaxMel_3s_CRF16sSegment_025Pct'
                lineColorBase=[0 0 1];
            case 'MelanopsinMR_SplatterControlPulse_3s_CRF16sSegment_025Pct'
                lineColorBase=[1 0 0];
        end
        
        for mm = 1:NStimTypes
            if mm==1
                dataOnly=false;
            else
                dataOnly=true;
            end            
            subplotHandle=subplot(NSessionsMerged/2,NStimTypes,(ss-1)*NStimTypes+mm);
            plotTitle=[ mergedPacketCellArray{thisPacketIndex}{1}.metaData.projectName, strrep(mergedPacketCellArray{thisPacketIndex}{1}.metaData.subjectName, '_', '\_')];
            pupilPMEL_PlotEvokedResponse( subplotHandle, avgPackets{thisPacketIndex, mm}.response.timebase, 100*avgPackets{thisPacketIndex, mm}.response.values, 100*avgPackets{thisPacketIndex, mm}.response.sem, 'dataOnly', dataOnly, 'ylim', [-50 50], 'xAxisAspect', 1, 'yAxisAspect', 2, 'lineColor', lineColorBase, 'plotTitle', plotTitle);
            hold on
            pupilPMEL_PlotEvokedResponse( subplotHandle, avgPackets{thisPacketIndex, mm}.response.timebase, twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.values, [], 'dataOnly', dataOnly, 'ylim', [-50 50], 'xAxisAspect', 1, 'yAxisAspect', 2, 'lineColor', [1 0 0],'plotTitle', plotTitle);
        end % Loop over stim types
        
    end % Loop over subjects
    
    %             % Save the plot. If the saving directory doesn't exist, create it.
    %         outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory,mergedPacketCellArray{ss}{1}.metaData.subjectName);
    %         if ~exist(outDir, 'dir')
    %             mkdir(outDir);
    %         end
    %         saveas(figHandle, fullfile(outDir, [mergedPacketCellArray{ss}{1}.metaData.projectName '.pdf']), 'pdf');
    %         close(figHandle);
    %
    %
    %     %% Plot the across-subject data LMS
    %     sessionTypes={'LMS','Mel'};
    %     sessionSets=[1 2 3 4; 5 6 7 8];
    %     sessionTypeColor(1,:,:)=[0.5 0.5 0.5; 0.625 0.375 0.375; 0.75 0.25 0.25; 0.875 0.125 0.125; 1 0 0];
    %     sessionTypeColor(2,:,:)=[0.5 0.5 0.5; 0.375 0.375 0.625; 0.25 0.25 0.75; 0.125 0.125 0.875; 0 0 1; ];
    %
    %     for pp=1:length(sessionTypes)
    %         figHandle = figure;
    %         plot([avgPackets{sessionSets(pp,1), 1}.response.timebase(1) avgPackets{sessionSets(pp,1), 1}.response.timebase(end)], [0 0], '-k');
    %         hold on
    %         acrossSubjectResponseMatrix=[];
    %         for cc=1:size(avgPackets,2)-1
    %             for ss=1:length(sessionSets(pp,:))
    %                 acrossSubjectResponseMatrix(ss,:)=avgPackets{sessionSets(pp,ss),cc}.response.values;
    %             end % loop over subjects in the session set
    %
    %             % Calculate the mean response and plot it
    %             meanAcrossSubjectResponse=nanmean(acrossSubjectResponseMatrix);
    %             semAcrossSubjectResponse=nanstd(acrossSubjectResponseMatrix)/sqrt(length(sessionSets(pp,:)));
    %             plot(avgPackets{sessionSets(pp,ss),cc}.response.timebase, 100*meanAcrossSubjectResponse,...
    %                 'LineWidth',3,'Color',sessionTypeColor(pp,cc,:));
    %         end % loop over stimulus types
    %
    %         % Add the stimulus profile to the plot
    %         plot(avgPackets{sessionSets(pp,ss),cc}.stimulus.timebase, 10*avgPackets{sessionSets(pp,ss),cc}.stimulus.values,'Color',[0 1 0]);
    %
    %         % Clean up the plot
    %         ylim(100*[-0.4 0.2]);
    %         pbaspect([1 1 1]);
    %         xlabel('Time [msecs]');
    %         ylabel('Amplitude [%]');
    %         adjustPlot(figHandle);
    %         title(['Average subject response - ' sessionTypes{pp}]);
    %         hold off
    %
    %         % Save the plot.
    %         outFileName = fullfile(dropboxAnalysisDir,subAnalysisDirectory,['AcrossSubject_' sessionTypes{pp} '.pdf']);
    %         saveas(figHandle, outFileName, 'pdf');
    %         close(figHandle);
    
end % Loop over session sets (Mel and LMS)


end

