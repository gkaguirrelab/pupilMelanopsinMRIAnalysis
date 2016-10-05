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
packetCacheBehavior='make';
packetCellArrayTag='maxMelLSM_CRF_Pupil';
packetCellArrayHash='';

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




%% Assemble the different stimTypes
for ss = 1:NSessionsMerged
    % Make an 'accumulator'
    NStimTypes = 6;
    for mm = 1:NStimTypes
        accumStimTypesResp{ss, mm} = [];
        accumStimTypesStim{ss, mm} = [];
    end
    
    NRunsTotal = length(mergedPacketCellArray{ss});
    for ii = 1:NRunsTotal
        % Find the stimulus onsets so that we can align the data to it. We
        % do that by finding a [0 1] edge from a difference operator.
        tmp = sum(mergedPacketCellArray{ss}{ii}.stimulus.values);
        tmp2 = diff(tmp);
        tmp2(tmp2 < 0) = 0;
        tmp2(tmp2 > 0) = 1;
        stimOnsets = strfind(tmp2, [0 1]);
        
        % Get the number of segments from the stimulus onsets
        NSegments = length(stimOnsets);
        t = (0:extractionDurInd)/1000;
        for jj = 1:length(stimOnsets)
            if (stimOnsets(jj)+extractionDurInd) <= length(mergedPacketCellArray{ss}{ii}.response.values)
                idxToExtract = stimOnsets(jj):(stimOnsets(jj)+extractionDurInd);
            else
                idxToExtract = stimOnsets(jj):length(mergedPacketCellArray{ss}{ii}.response.values);
            end
            thisPacket.packetType = 'pupil';
            thisPacket.stimulusFile = mergedPacketCellArray{ss}{ii}.metaData.stimulusFile;
            thisPacket.responseFile = mergedPacketCellArray{ss}{ii}.metaData.responseFile;
            thisPacket.respValues = mergedPacketCellArray{ss}{ii}.response.values(idxToExtract);
            % Normalize the pupil data
            thisPacket.respValues = (thisPacket.respValues - nanmean(thisPacket.respValues(1:normalizationDurInd)))./nanmean(thisPacket.respValues(1:normalizationDurInd));
            thisPacket.respTimeBase = mergedPacketCellArray{ss}{ii}.response.timebase(idxToExtract);
            thisPacket.respTimeBase = thisPacket.respTimeBase-thisPacket.respTimeBase(1);
            thisPacket.stimValues = mergedPacketCellArray{ss}{ii}.stimulus.values(jj, idxToExtract);
            thisPacket.stimTimeBase = mergedPacketCellArray{ss}{ii}.stimulus.timebase(idxToExtract);
            thisPacket.stimTimeBase =  thisPacket.stimTimeBase - thisPacket.stimTimeBase(1);
            thisPacket.stimMetaData.stimTypes = mergedPacketCellArray{ss}{ii}.stimulus.metaData.stimTypes(jj);
            thisPacket.stimMetaData.stimLabels = params.stimMetaData.stimLabels;
            
            % Could make packets here for each event, but not doing it...
            % Just making the 'accumStimTypes' variable for now.
            accumStimTypesResp{ss, thisPacket.stimMetaData.stimTypes} = [accumStimTypesResp{ss, thisPacket.stimMetaData.stimTypes} ; thisPacket.respValues];
            accumStimTypesStim{ss, thisPacket.stimMetaData.stimTypes} = [accumStimTypesStim{ss, thisPacket.stimMetaData.stimTypes} ; thisPacket.stimValues];
        end
    end
end

%% Make average packets per subject
for ss = 1:NSessionsMerged
    for mm = 1:NStimTypes
        thisPacket = [];
        thisPacket.packetType = 'pupil';
        thisPacket.sessionDir = '';
        thisPacket.stimulusFile = mergedPacketCellArray{ss}{1}.metaData.stimulusFile;
        thisPacket.responseFile = mergedPacketCellArray{ss}{1}.metaData.responseFile;
        thisPacket.respValues =  nanmean(accumStimTypesResp{ss, mm});
        thisPacket.respTimeBase = mergedPacketCellArray{ss}{1}.response.timebase(idxToExtract);
        thisPacket.respTimeBase = thisPacket.respTimeBase-thisPacket.respTimeBase(1);
        thisPacket.stimValues = max(accumStimTypesStim{ss, mm});
        thisPacket.stimTimeBase = mergedPacketCellArray{ss}{1}.stimulus.timebase(idxToExtract);
        thisPacket.stimTimeBase = thisPacket.stimTimeBase-thisPacket.stimTimeBase(1);
        thisPacket.stimMetaData.stimTypes = params.stimMetaData.stimTypes(jj);
        thisPacket.stimMetaData.stimLabels = params.stimMetaData.stimLabels;
        avgPackets{ss, mm} = makePacket(thisPacket);
    end
end

%% Fit model to avg packets
twoComponentFitToData = fitTPUPModelToAverageResponse(avgPackets);

%% call out to fitIAMPModelToIndividualResponse here
[~] = fitIAMPModelToIndividualResponse(mergedPacketCellArray);

%% Plot the data
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
