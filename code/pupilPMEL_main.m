%% pupilMelanopsinMRIAnalysis.m
%
% Program to run data analysis on pupil data collected within the
% `MelanopsinMR` project.


% Housekeeping
clearvars; close all; clc;
warning on;

% Discover user name and set Dropbox path
[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxAnalysisDir = ...
    fullfile('/Users', userName, ...
    '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/pupilMelanopsinMRIAnalysis');

% Define cache behavior.

packetCacheBehavior='load';
packetCacheTag='maxMelLSM_CRF_Pupil';
packetCacheHash='82a4db868b14431c9c4cc65795a2c58d';

fitTPUPCacheBehavior='load';
fitTPUPCacheTag='TPUPModelFits';
fitTPUPCacheHash='1b052ee837c7576303b4f88053491300';


%% Create or load the packetCellArray
switch packetCacheBehavior
    case 'make'  % If we are not to load the mergedPacketCellArray, then we must generate it
        % Make the packetCellArray
        [ mergedPacketCellArray ] = pupilPMEL_makeMergedPacketCellArray( userName );
        % calculate the hex MD5 hash for the packetCellArray
        packetCacheHash = DataHash(mergedPacketCellArray);
        % Set path to the packetCache and save it using the MD5 hash name
        packetCacheFileName=fullfile(dropboxAnalysisDir, 'packetCache', [packetCacheTag '_' packetCacheHash '.mat']);
        save(packetCacheFileName,'mergedPacketCellArray','-v7.3');
        fprintf(['Saved the ' packetCacheTag ' with hash ID ' packetCacheHash '\n']);
    case 'load'  % load a cached packetCellArray
        fprintf('>> Loading cached packetCellArray\n');
        packetCacheFileName=fullfile(dropboxAnalysisDir, 'packetCache', [packetCacheTag '_' packetCacheHash '.mat']);
        load(packetCacheFileName);
    otherwise
        error('Please define a legal packetCacheBehavior');
end

%% Fit TPUP model to avg packets
switch fitTPUPCacheBehavior    
    case 'make'
        [ twoComponentFitToData ] = fitTPUPModelToAverageResponse(...
            mergedPacketCellArray, ...
            dropboxAnalysisDir);
        % calculate the hex MD5 hash for the twoComponentFitToData
        fitTPUPCacheHash = DataHash(twoComponentFitToData);        
        % Set path to the packetCache and save it using the MD5 hash name
        fitTPUPCacheFileName=fullfile(dropboxAnalysisDir, 'analysisCache', [fitTPUPCacheTag '_' fitTPUPCacheHash '.mat']);
        save(fitTPUPCacheFileName,'twoComponentFitToData','-v7.3');
        fprintf(['Saved the ' fitTPUPCacheTag ' with hash ID ' fitTPUPCacheHash '\n']);        
    case 'load'  % load a cached twoComponentFitToData        
        fprintf(['>> Loading cached ' fitTPUPCacheTag ' \n']);
        fitTPUPCacheFileName=fullfile(dropboxAnalysisDir, 'analysisCache', [fitTPUPCacheTag '_' fitTPUPCacheHash '.mat']);
        load(fitTPUPCacheFileName);        
    otherwise        
        error('Please define a legal packetCacheBehavior');
end


%% Plot the average pupil responses
pupilPMEL_plotAveragePupilResponses( mergedPacketCellArray, twoComponentFitToData, dropboxAnalysisDir )

% %% Fit FCON model to individual events
% [myResultsVariable] = ...
%     fitFCONModelToIndividualResponses(mergedPacketCellArray, ...
%     twoComponentFitToData, ...
%     dropboxAnalysisDir);
%
% % plot FCON results
% plotFCONResults(mergedPacketCellArray, dropboxAnalysisDir, myResultsVariable);
%
