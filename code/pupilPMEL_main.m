function pupilPMEL_main(ppsRawDataDir, ppsPupilPacketsDir, analysisDir)

fprintf('> Running %s\n', mfilename);

%
% Program to run data analysis on pupil data collected within the
% `MelanopsinMR` project.

% Housekeeping
close all; clc;
warning on;

% Define output folders for figures and tables
outFigureDir = fullfile(analysisDir, 'figures');
if ~isdir(outFigureDir)
    mkdir(outFigureDir);
end
outTableDir = fullfile(analysisDir, 'tables');
if ~isdir(outTableDir)
    mkdir(outTableDir);
end

% Define cache and analysis behavior.

packetCacheBehavior='load';
packetCacheTag='maxMelLSM_CRF_Pupil';
packetCacheHash='f3e507e694a98311485df5adadfc373a';

fitTPUPCacheBehavior='load';
fitTPUPCacheTag='TPUPModelFits';
fitTPUPCacheHash='b2bd2f43419d0ef28f985fe91a8f0957';

refitLockedTPUPCacheBehavior='load';
refitLockedTPUPCacheTag='refitTPUPLockedTimeFits';
refitLockedTPUPCacheHash='e4e6677382753e3ae1d7c628a7ebad4b';

makePupilPlots='make';
analyzeBlinksBehavior='make';

%% Create or load the packetCellArray
switch packetCacheBehavior
    case 'make'  % If we are not to load the mergedPacketCellArray, then we must generate it
        % Make the packetCellArray
        [ mergedPacketCellArray ] = pupilPMEL_makeMergedPacketCellArray(ppsRawDataDir);
        % calculate the hex MD5 hash for the packetCellArray
        packetCacheHash = DataHash(mergedPacketCellArray);
        % Set path to the packetCache and save it using the MD5 hash name
        packetCacheFileName=fullfile(ppsPupilPacketsDir, 'packetCache', [packetCacheTag '_' packetCacheHash '.mat']);
        save(packetCacheFileName,'mergedPacketCellArray','-v7.3');
        fprintf(['Saved the ' packetCacheTag ' with hash ID ' packetCacheHash '\n']);
    case 'load'  % load a cached packetCellArray
        fprintf('>> Loading cached packetCellArray\n');
        packetCacheFileName=fullfile(ppsPupilPacketsDir, 'packetCache', [packetCacheTag '_' packetCacheHash '.mat']);
        load(packetCacheFileName);
    otherwise
        error('Please define a legal packetCacheBehavior');
end

%% Fit TPUP model to avg packets
switch fitTPUPCacheBehavior    
    case 'make'
        [ twoComponentFitToData ] = pupilPMEL_fitTPUPModelToAverageResponse(...
            mergedPacketCellArray, ...
            ppsPupilPacketsDir);
        % calculate the hex MD5 hash for the twoComponentFitToData
        fitTPUPCacheHash = DataHash(twoComponentFitToData);        
        % Set path to the packetCache and save it using the MD5 hash name
        fitTPUPCacheFileName=fullfile(ppsPupilPacketsDir, 'analysisCache', [fitTPUPCacheTag '_' fitTPUPCacheHash '.mat']);
        save(fitTPUPCacheFileName,'twoComponentFitToData','-v7.3');
        fprintf(['Saved the ' fitTPUPCacheTag ' with hash ID ' fitTPUPCacheHash '\n']);        
    case 'load'  % load a cached twoComponentFitToData        
        fprintf(['>> Loading cached ' fitTPUPCacheTag ' \n']);
        fitTPUPCacheFileName=fullfile(ppsPupilPacketsDir, 'analysisCache', [fitTPUPCacheTag '_' fitTPUPCacheHash '.mat']);
        load(fitTPUPCacheFileName);        
    otherwise        
        error('Please define a legal fitTPUPCacheBehavior');
end

%% Fit TPUP model to avg packets with pinned time parameters
switch refitLockedTPUPCacheBehavior    
    case 'make'
        [ lockedTimeTPUPFitToData ] = pupilPMEL_refitLockedTPUPToAverageResponse(...
            mergedPacketCellArray, ...
            ppsPupilPacketsDir, ...
            twoComponentFitToData);
        % calculate the hex MD5 hash for the twoComponentFitToData
        refitLockedTPUPCacheHash = DataHash(lockedTimeTPUPFitToData);        
        % Set path to the packetCache and save it using the MD5 hash name
        refitLockedTPUPCacheFileName=fullfile(ppsPupilPacketsDir, 'analysisCache', [refitLockedTPUPCacheTag '_' refitLockedTPUPCacheHash '.mat']);
        save(refitLockedTPUPCacheFileName,'lockedTimeTPUPFitToData','-v7.3');
        fprintf(['Saved the ' refitLockedTPUPCacheTag ' with hash ID ' refitLockedTPUPCacheHash '\n']);        
    case 'load'  % load a cached refit TPUP with locked temporal components        
        fprintf(['>> Loading cached ' refitLockedTPUPCacheTag ' \n']);
        refitLockedTPUPCacheFileName=fullfile(ppsPupilPacketsDir, 'analysisCache', [refitLockedTPUPCacheTag '_' refitLockedTPUPCacheHash '.mat']);
        load(refitLockedTPUPCacheFileName);        
    otherwise        
        error('Please define a legal refitLockedTPUPCacheBehavior');
end

%% Plot the average pupil responses
switch makePupilPlots
    case 'make'
        pupilPMEL_plotAveragePupilResponses(mergedPacketCellArray, twoComponentFitToData, lockedTimeTPUPFitToData, analysisDir)
end

% Analyze blinks and gaze
switch analyzeBlinksBehavior
    case 'make'
        pupilPMEL_analyzeGaze(mergedPacketCellArray, analysisDir)
end