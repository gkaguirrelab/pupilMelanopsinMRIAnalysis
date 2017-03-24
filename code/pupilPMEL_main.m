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
    '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/pupilMelanopsinMRIAnalysis');

% Define packetCacheBehavior. Options include:
%    'make' - load and process stim/response files, save the packets
%    'load' - load the packets from the passed hash name
packetCacheBehavior='load';
packetCellArrayTag='maxMelLSM_CRF_Pupil';
packetCellArrayHash='d702e5898732600707fb32a2302d4772';

%% Create or load the packetCellArray
switch packetCacheBehavior
    
    case 'make'  % If we are not to load the mergedPacketCellArray, then we must generate it
        
        % Make the packetCellArray
        [ mergedPacketCellArray ] = pupilPMEL_makeMergedPacketCellArray( userName );
        
        % calculate the hex MD5 hash for the packetCellArray
        packetCellArrayHash = DataHash(mergedPacketCellArray);
        
        % Set path to the packetCache and save it using the MD5 hash name
        packetCacheFileName=fullfile(dropboxAnalysisDir, 'packetCache', [packetCellArrayTag '_' packetCellArrayHash '.mat']);
        save(packetCacheFileName,'mergedPacketCellArray','-v7.3');
        fprintf(['Saved the packetCellArray with hash ID ' packetCellArrayHash '\n']);
        
    case 'load'  % load a cached packetCellArray
        
        fprintf('>> Loading cached packetCellArray\n');
        packetCacheFileName=fullfile(dropboxAnalysisDir, 'packetCache', [packetCellArrayTag '_' packetCellArrayHash '.mat']);
        load(packetCacheFileName);
        
    otherwise
        
        error('Please define a legal packetCacheBehavior');
end


%% Fit TPUP model to avg packets
twoComponentFitToData = ...
    fitTPUPModelToAverageResponse(...
    mergedPacketCellArray, ...
    dropboxAnalysisDir);

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

