%% MelanopsinMR_PupilAnalysis.m
%
% Program to run data analysis on pupil data collected within the
% `MelanopsinMR` project.

% Housekeeping
clearvars; close all; clc;

% Discover user name and set Dropbox path
[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data');
saveDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis');

% Define the session directories
sessDirs = {...
    'MelanopsinMRMaxLMSCRF/HERO_asb1/060816'...
    'MelanopsinMRMaxLMSCRF/HERO_aso1/060116' ...
    'MelanopsinMRMaxLMSCRF/HERO_gka1/060616' ...
    'MelanopsinMRMaxLMSCRF/HERO_mxs1/062416' ...
    'MelanopsinMRMaxLMSCRF/HERO_mxs1/062816' ...
    'MelanopsinMRMaxMelCRF/HERO_asb1/060716' ...
    'MelanopsinMRMaxMelCRF/HERO_aso1/053116' ...
    'MelanopsinMRMaxMelCRF/HERO_gka1/060216' ...
    'MelanopsinMRMaxMelCRF/HERO_mxs1/060916' ...
    'MelanopsinMRMaxMelCRF/HERO_mxs1/061016'
    };

%sessDirs = {'MelanopsinMRMaxLMSCRF/HERO_gka1/060616'};
% Define which sessions we'd like to merge
whichSessionsToMerge = {[1], [2], [3], [4 5], [6], [7], [8], [9 10]};
%whichSessionsToMerge = {[1]};

% Set some parameters we need
normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;

for ss = 1:length(sessDirs)
    % Clear some information
    Data_Per_Segment = [];
    uniqueCombos = [];
    allDataTmp = [];
    allIndicesTmp = [];
    
    % Extract some information about this session
    tmp = strsplit(sessDirs{ss}, '/');
    params.sessionType = tmp{1};
    params.sessionObserver = tmp{2};
    params.sessionDate = tmp{3};
    
    % Display some useful information
    fprintf('>> Processing <strong>%s</strong> | <strong>%s</strong> | <strong>%s</strong>\n', params.sessionType, params.sessionObserver, params.sessionDate);
    
    % Determine some parameters
    switch params.sessionDate
        case {'053116' '060116' '060216'}
            params.acquisitionFreq      = 30;
        otherwise
            params.acquisitionFreq      = 60;
    end
    params.LiveTrackSamplingRate        = 60; % Hz
    params.ResamplingFineFreq           = 1000; % 1 msec
    params.BlinkWindowSample            = -50:50; % Samples surrounding the blink event
    params.TRDurSecs                    = 0.8;
    
    % Make the packets
    params.packetType       = 'pupil';
    params.sessionDir       = fullfile(dropboxDir, sessDirs{ss});
    NRuns = length(listdir(fullfile(params.sessionDir, 'MatFiles', '*.mat'), 'files'));
    
    % Iterate over runs
    for ii = 1:NRuns;
        fprintf('\t* Run <strong>%g</strong> / <strong>%g</strong>\n', ii, NRuns);
        % Set up some parameters
        params.runNum           = ii;
        params.stimulusFile     = fullfile(params.sessionDir, 'MatFiles', [params.sessionObserver '-' params.sessionType '-' num2str(ii, '%02.f') '.mat']);
        params.responseFile     = fullfile(params.sessionDir, 'EyeTrackingFiles', [params.sessionObserver '-' params.sessionType '-' num2str(ii, '%02.f') '.mat']);
        [params.respValues params.respTimeBase] = loadPupilDataForPackets(params);
        [params.stimValues params.stimTimeBase params.stimMetaData] = makeStimStruct(params);
        packets{ss, ii} = makePacket(params);
    end
    fprintf('\n');
end

%% Merge sessions
NSessionsMerged = length(whichSessionsToMerge);
for mm = 1:NSessionsMerged
    mergeIdx = whichSessionsToMerge{mm};
    mergedPacket = {packets{mergeIdx, :}};
    mergedPacket = mergedPacket(~cellfun('isempty', mergedPacket));
    mergedPackets{mm} = mergedPacket;
end

%% Assemble the different stimTypes
for ss = 1:NSessionsMerged
    % Make an 'accumulator'
    NStimTypes = 6;
    for mm = 1:NStimTypes
        accumStimTypesResp{ss, mm} = [];
        accumStimTypesStim{ss, mm} = [];
    end
    
    NRunsTotal = length(mergedPackets{ss});
    for ii = 1:NRunsTotal
        % Find the stimulus onsets so that we can align the data to it. We
        % do that by finding a [0 1] edge from a difference operator.
        tmp = sum(mergedPackets{ss}{ii}.stimulus.values);
        tmp2 = diff(tmp);
        tmp2(tmp2 < 0) = 0;
        tmp2(tmp2 > 0) = 1;
        stimOnsets = strfind(tmp2, [0 1]);
        
        % Get the number of segments from the stimulus onsets
        NSegments = length(stimOnsets);
        t = (0:extractionDurInd)/1000;
        for jj = 1:length(stimOnsets)
            if (stimOnsets(jj)+extractionDurInd) <= length(mergedPackets{ss}{ii}.response.values)
                idxToExtract = stimOnsets(jj):(stimOnsets(jj)+extractionDurInd);
            else
                idxToExtract = stimOnsets(jj):length(mergedPackets{ss}{ii}.response.values);
            end
            thisPacket.packetType = 'pupil';
            thisPacket.stimulusFile = mergedPackets{ss}{ii}.metaData.stimulusFile;
            thisPacket.responseFile = mergedPackets{ss}{ii}.metaData.responseFile;
            thisPacket.respValues = mergedPackets{ss}{ii}.response.values(idxToExtract);
            % Normalize the pupil data
            thisPacket.respValues = (thisPacket.respValues - nanmean(thisPacket.respValues(1:normalizationDurInd)))./nanmean(thisPacket.respValues(1:normalizationDurInd));
            thisPacket.respTimeBase = mergedPackets{ss}{ii}.response.timebase(idxToExtract);
            thisPacket.respTimeBase = thisPacket.respTimeBase-thisPacket.respTimeBase(1);
            thisPacket.stimValues = mergedPackets{ss}{ii}.stimulus.values(jj, idxToExtract);
            thisPacket.stimTimeBase = mergedPackets{ss}{ii}.stimulus.timebase(idxToExtract);
            thisPacket.stimTimeBase =  thisPacket.stimTimeBase - thisPacket.stimTimeBase(1);
            thisPacket.stimMetaData.stimTypes = mergedPackets{ss}{ii}.stimulus.metaData.stimTypes(jj);
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
        thisPacket.stimulusFile = mergedPackets{ss}{1}.metaData.stimulusFile;
        thisPacket.responseFile = mergedPackets{ss}{1}.metaData.responseFile;
        thisPacket.respValues =  nanmean(accumStimTypesResp{ss, mm});
        thisPacket.respTimeBase = mergedPackets{ss}{1}.response.timebase(idxToExtract);
        thisPacket.respTimeBase = thisPacket.respTimeBase-thisPacket.respTimeBase(1);
        thisPacket.stimValues = max(accumStimTypesStim{ss, mm});
        thisPacket.stimTimeBase = mergedPackets{ss}{1}.stimulus.timebase(idxToExtract);
        thisPacket.stimTimeBase = thisPacket.stimTimeBase-thisPacket.stimTimeBase(1);
        thisPacket.stimMetaData.stimTypes = params.stimMetaData.stimTypes(jj);
        thisPacket.stimMetaData.stimLabels = params.stimMetaData.stimLabels;
        avgPackets{ss, mm} = makePacket(thisPacket);
    end
end

%% Fit model to avg packets
% Insert code here. Average packets are now in avgPackets{ss, mm}, where ss
% is subject (with merged session), and mm is stimulus type

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instamce.
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
        temporalFit = tmriTwoComponentPupilResponse();
        
        % Grab a single packet
        singlePacket=avgPackets{ss, mm};
        
        % For now, we need to define an empty HRF field for the
        % temporalFittingEnginge. This should go away soon.
        singlePacket.HRF.values=[];
        singlePacket.HRF.timebase=[];
        
        % report fitting progress
        fprintf('iterations:');
        
        % Conduct the fit
        [paramsFit,fVal,fitResponse] = temporalFit.fitResponse(singlePacket, 'DefaultParamsInfo', defaultParamsInfo, ...
            'paramLockMatrix',paramLockMatrix);
        
        % Store the fitResponse
        twoComponentFitToData{ss,mm}.paramsFit=paramsFit;
        twoComponentFitToData{ss,mm}.fVal=fVal;
        twoComponentFitToData{ss,mm}.fitResponse=fitResponse;
        
        % Report the fit error value:
        fprintf('\n\t> Fit error value: %g', fVal);
        fprintf('\n');
        
        % Clear the object
        delete(temporalFit); %clear temporalFit;
        
    end % loop over subjects
end % loop over stimuli

%% Plot the data
for ss = 1:NSessionsMerged
    plotFig = figure;
    for mm = 1:NStimTypes
        plot([avgPackets{ss, mm}.response.timebase(1) avgPackets{ss, mm}.response.timebase(end)], [0 0], '-k'); hold on;
        % plot a model fit if it is available
         if ~isempty(twoComponentFitToData{ss,mm}) & ~isempty(twoComponentFitToData{ss,mm}.fitResponse)
            plot(avgPackets{ss, mm}.response.timebase, 100*twoComponentFitToData{ss,mm}.fitResponse,'--k');
        end
        plot(avgPackets{ss, mm}.response.timebase, 100*avgPackets{ss, mm}.response.values);
        xlim([avgPackets{ss, mm}.response.timebase(1) avgPackets{ss, mm}.response.timebase(end)]);
    end
    ylim(100*[-0.5 0.5]);
    pbaspect([1 1 1]);
    xlabel('Time [msecs]');
    ylabel('Amplitude [%]');
    adjustPlot(plotFig);
    title({ mergedPackets{ss}{1}.metaData.projectName, strrep(mergedPackets{ss}{1}.metaData.subjectName, '_', '\_')});
    
    % Save the plot. If the saving directory doesn't exist, create it.
    outDir = fullfile(saveDir, mergedPackets{ss}{1}.metaData.projectName, mergedPackets{ss}{1}.metaData.subjectName);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveas(plotFig, fullfile(outDir, [mergedPackets{ss}{1}.metaData.projectName '.png']), 'png');
    close(plotFig);
end