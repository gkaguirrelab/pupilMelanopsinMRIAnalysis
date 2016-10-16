function [twoComponentFitToData] = fitTPUPModelToAverageResponse(mergedPacketCellArray, dropboxAnalysisDir)

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag );

%% SETUP

% Assign a name for this sub-analysis directory for saving plots and fits
subAnalysisDirectory='fitTPUPModelToAverageResponse';

% Identify the number of (sessions x subjects) and the number of stimulus
% types from the dimensionality of the avgPackets cell array

NSessionsMerged=size(avgPackets,1);
NStimTypes=size(avgPackets,2);

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;

% Announce what we are about to do
fprintf('>> Fitting two-component model to pupil data (TPUP)\n');

% Pre-allocate a structure to hold the fit results
twoComponentFitToData{NSessionsMerged,NStimTypes}.paramsFit=[];
twoComponentFitToData{NSessionsMerged,NStimTypes}.fVal=[];
twoComponentFitToData{NSessionsMerged,NStimTypes}.fitResponse=[];

% Loop over subjects and stimulus types
% Skipping the attention task for now
for ss = 1:NSessionsMerged
    for mm = 1:NStimTypes-1
        % Update the user
        fprintf('* Subject, stimulus <strong>%g</strong> , <strong>%g</strong>', ss, mm);
        fprintf('\n');
        
        % Construct the model object
        temporalFit = tfeTPUP('verbosity','none');
        
        % Grab a single packet
        singlePacket=avgPackets{ss, mm};
        
        % Conduct the fit
        [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo, ...
            'paramLockMatrix',paramLockMatrix);
        
        % Store the fitResponse
        twoComponentFitToData{ss,mm}.paramsFit=paramsFit;
        twoComponentFitToData{ss,mm}.fVal=fVal;
        twoComponentFitToData{ss,mm}.modelResponseStruct=modelResponseStruct;
        
        % Clear the object
        delete(temporalFit);
        
    end % loop over subjects
end % loop over stimuli


%% Plot the data

NSessionsMerged=size(mergedPacketCellArray,2);
NStimTypes=6;

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


end % function