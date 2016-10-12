%% pupilMelanopsinMRIAnalysis.m
%
% Program to run data analysis on pupil data collected within the
% `MelanopsinMR` project.

% harry added line 84,85 to get pupil units back into mm

% Housekeeping
clearvars; close all; clc;
warning on;

% Discover user name and set Dropbox path
[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxAnalysisDir = ...
    fullfile('/Users', userName, ...
    '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/pupilMelanopsinMRIAnalysis/packetCache');


% Define packetCacheBehavior. Options include:
%    'make' - load and process stim/response files, save the packets
%    'load' - load the packets from the passed hash name
packetCacheBehavior='load';
packetCellArrayTag='maxMelLSM_CRF_Pupil';
packetCellArrayHash='0460faac0fa14f9837e5b7a6775f6258';

%% Create or load the packetCellArray
switch packetCacheBehavior
    
    case 'make'  % If we are not to load the mergedPacketCellArray, then we must generate it
        
        % Make the packetCellArray
        [ mergedPacketCellArray ] = pupilPMEL_makeMergedPacketCellArray( userName );
        
        % calculate the hex MD5 hash for the packetCellArray
        packetCellArrayHash = DataHash(mergedPacketCellArray);
        
        % Set path to the packetCache and save it using the MD5 hash name
        packetCacheFileName=fullfile(dropboxAnalysisDir, [packetCellArrayTag '_' packetCellArrayHash '.mat']);
        save(packetCacheFileName,'mergedPacketCellArray','-v7.3');
        fprintf(['Saved the packetCellArray with hash ID ' packetCellArrayHash '\n']);
        
    case 'load'  % load a cached packetCellArray
        
        fprintf('>> Loading cached packetCellArray\n');
        packetCacheFileName=fullfile(dropboxAnalysisDir, [packetCellArrayTag '_' packetCellArrayHash '.mat']);
        load(packetCacheFileName);
        
    otherwise
        
        error('Please define a legal packetCacheBehavior');
end

%% Obtain the responses
fprintf('>> Making average responses\n');
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray );

%% Fit model to avg packets
%twoComponentFitToData = fitTPUPModelToAverageResponse(avgPackets);

%% Fit IAMP model to individual events
fprintf('>> Fitting IAMP model to individual responses\n');
[~] = fitIAMPModelToIndividualResponse(mergedPacketCellArray);

%% Fit FCON model to individual events
[myResultsVariable] = fitFCONModelToIndividualResponses(mergedPacketCellArray, twoComponentFitToData);


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
    outDir = fullfile(dropboxAnalysisDir, mergedPacketCellArray{ss}{1}.metaData.projectName, mergedPacketCellArray{ss}{1}.metaData.subjectName);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveas(plotFig, fullfile(outDir, [mergedPacketCellArray{ss}{1}.metaData.projectName '.png']), 'png');
    saveas(plotFig, fullfile(outDir, [mergedPacketCellArray{ss}{1}.metaData.projectName '.pdf']), 'pdf');
    close(plotFig);
end
