function [ ] = pupilPMEL_plotAveragePupilResponses( mergedPacketCellArray, twoComponentFitToData, dropboxAnalysisDir )
% Create plots of average evoked pupil responses within and across subjects

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, 'lowFreqClean', true, 'aggregateMethod', 'mean');

%% Plot the within-subject data
NSessionsMerged=size(mergedPacketCellArray,2);
NStimTypes=5;

for ee = 1:2 % loop over experiments, LMS and Mel
    figHandle = figure;
    set(gcf, 'PaperSize', [8.5 11]);
    for ss = 1:NSessionsMerged/2 % loop over subjects
        
        thisPacketIndex=(ee-1)*(NSessionsMerged/2)+ss;
        % Identify the stimulus direction (LMS or Mel)
        experimentName=mergedPacketCellArray{thisPacketIndex}{1}.metaData.projectName;
        lineColorBase=[.25 .25 .25];
        switch experimentName
            case 'MaxLMSCRF'
                lineColorBase=[.25 .25 .25];
            case 'MaxMelCRF'
                lineColorBase=[0 0 1];
        end
        
        for mm = 1:NStimTypes
            if mm==1
                dataOnly=false;
            else
                dataOnly=true;
            end
            subplotHandle=subplot(NSessionsMerged/2,NStimTypes,(ss-1)*NStimTypes+mm);
            plotTitle=mergedPacketCellArray{thisPacketIndex}{1}.metaData.subjectName;
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                100*avgPackets{thisPacketIndex, mm}.response.values, ...
                100*avgPackets{thisPacketIndex, mm}.response.sem, ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 10], ...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', lineColorBase, 'plotTitle', plotTitle);
            hold on
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.values, ...
                [], ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 10], ...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', [1 0 0],'plotTitle', plotTitle);
        end % loop over stim types
        
    end % loop over subjects
    
    % Save the plots
    plotFileName=fullfile(dropboxAnalysisDir, 'Figures', [experimentName '_PupilTrialMeanResponses.pdf']);
    pupilPMEL_suptitle(figHandle,[experimentName ' - Pupil responses']);
    set(gca,'FontSize',6);
    set(figHandle,'Renderer','painters');
    print(figHandle, plotFileName, '-dpdf', '-fillpage');
    close(figHandle);
    
end % loop over session sets (Mel and LMS)

end % function