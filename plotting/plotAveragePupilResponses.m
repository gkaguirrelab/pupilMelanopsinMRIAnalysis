function [ ] = plotAveragePupilResponses( mergedPacketCellArray, twoComponentFitToData, dropboxAnalysisDir )
% Create plots of average evoked pupil responses within and across subjects

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag );

%% Plot the within-subject data
NSessionsMerged=size(mergedPacketCellArray,2);
NStimTypes=6;

% Assign a name for this sub-analysis directory for saving plots and fits
subAnalysisDirectory='fitTPUPModelToAverageResponse';

for ss = 1:NSessionsMerged
    plotFig = figure;
    for mm = 1:NStimTypes
        plot([avgPackets{ss, mm}.response.timebase(1) avgPackets{ss, mm}.response.timebase(end)], [0 0], '-k'); hold on;
        % plot a model fit if it is available
        if isfield(twoComponentFitToData{ss,mm}, 'modelResponseStruct')
            plot(twoComponentFitToData{ss, mm}.modelResponseStruct.timebase, 100*twoComponentFitToData{ss,mm}.modelResponseStruct.values,'--k');
        end
        plot(avgPackets{ss, mm}.response.timebase, 100*avgPackets{ss, mm}.response.values);
        xlim([avgPackets{ss, mm}.response.timebase(1) avgPackets{ss, mm}.response.timebase(end)]);
    end
    ylim(100*[-0.5 0.5]);
    pbaspect([1 1 1]);
    xlabel('Time [msecs]');
    ylabel('Amplitude [%]');
    adjustPlot(plotFig);
    title({ mergedPacketCellArray{ss}{1}.metaData.projectName, strrep(mergedPacketCellArray{ss}{1}.metaData.subjectName, '_', '\_')});
    
    % Save the plot. If the saving directory doesn't exist, create it.
    outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory,mergedPacketCellArray{ss}{1}.metaData.subjectName);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveas(plotFig, fullfile(outDir, [mergedPacketCellArray{ss}{1}.metaData.projectName '.pdf']), 'pdf');
    close(plotFig);
end

%% Plot the across-subject data LMS
sessionTypes={'LMS','Mel'};
sessionSets=[1 2 3 4; 5 6 7 8];
sessionTypeColor(1,:,:)=[0.5 0.5 0.5; 0.625 0.375 0.375; 0.75 0.25 0.25; 0.875 0.125 0.125; 1 0 0];
sessionTypeColor(2,:,:)=[0.5 0.5 0.5; 0.375 0.375 0.625; 0.25 0.25 0.75; 0.125 0.125 0.875; 0 0 1; ];

for pp=1:length(sessionTypes)
    plotFig = figure;
    plot([avgPackets{sessionSets(pp,1), 1}.response.timebase(1) avgPackets{sessionSets(pp,1), 1}.response.timebase(end)], [0 0], '-k');
    hold on
    acrossSubjectResponseMatrix=[];
    for cc=1:size(avgPackets,2)-1
        for ss=1:length(sessionSets(pp,:))
            acrossSubjectResponseMatrix(ss,:)=avgPackets{sessionSets(pp,ss),cc}.response.values;
        end % loop over subjects in the session set
        
        % Calculate the mean response and plot it
        meanAcrossSubjectResponse=nanmean(acrossSubjectResponseMatrix);
        semAcrossSubjectResponse=nanstd(acrossSubjectResponseMatrix)/sqrt(length(sessionSets(pp,:)));
        plot(avgPackets{sessionSets(pp,ss),cc}.response.timebase, 100*meanAcrossSubjectResponse,...
            'LineWidth',3,'Color',sessionTypeColor(pp,cc,:));
    end % loop over stimulus types

    % Add the stimulus profile to the plot
    plot(avgPackets{sessionSets(pp,ss),cc}.stimulus.timebase, 10*avgPackets{sessionSets(pp,ss),cc}.stimulus.values,'Color',[0 1 0]);
    
    % Clean up the plot
    ylim(100*[-0.4 0.2]);
    pbaspect([1 1 1]);
    xlabel('Time [msecs]');
    ylabel('Amplitude [%]');
    adjustPlot(plotFig);
    title(['Average subject response - ' sessionTypes{pp}]);
    hold off

    % Save the plot.
    outFileName = fullfile(dropboxAnalysisDir,subAnalysisDirectory,['AcrossSubject_' sessionTypes{pp} '.pdf']);
    saveas(plotFig, outFileName, 'pdf');
    close(plotFig);

end % Loop over session sets (Mel and LMS)


end

