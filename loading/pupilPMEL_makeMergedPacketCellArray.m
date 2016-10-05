function [ mergedPacketCellArray, params ] = pupilPMEL_makeMergedPacketCellArray( userName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define the dropbox path
dropboxDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data');

% Define the session directories
sessDirs = {...
    'MelanopsinMRMaxLMSCRF/HERO_asb1/060816' ...
    'MelanopsinMRMaxLMSCRF/HERO_aso1/060116' ...
    'MelanopsinMRMaxLMSCRF/HERO_gka1/060616' ...
    'MelanopsinMRMaxLMSCRF/HERO_mxs1/062816' ...
    'MelanopsinMRMaxMelCRF/HERO_asb1/060716' ...
    'MelanopsinMRMaxMelCRF/HERO_aso1/053116' ...
    'MelanopsinMRMaxMelCRF/HERO_gka1/060216' ...
    'MelanopsinMRMaxMelCRF/HERO_mxs1/060916' ...
    'MelanopsinMRMaxMelCRF/HERO_mxs1/061016'
    };

sessDirs = {'MelanopsinMRMaxLMSCRF/HERO_gka1/060616'};
% Define which sessions we'd like to merge
whichSessionsToMerge = {[1], [2], [3], [4], [5], [6], [7], [8 9]};
whichSessionsToMerge = {[1]};

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
        params.lowFreqFilterHz  = 0.025; % The boundary of the low-pass filter that will be used to estimate the component of slow pupil modulation in the signal
        
        % new line here
        calibrationFileTemp  = cellstr(ls(fullfile(params.sessionDir, 'EyeTrackingFiles/ScaleCalibration', '*.mat')));
        params.calibrationFile = calibrationFileTemp{1,1};
        [params.respValues params.respTimeBase params.lowFreqComponentValues] = loadPupilDataForPackets(params);
        [params.stimValues params.stimTimeBase params.stimMetaData] = pupilMelanopsinMRImakeStimStruct(params);
        packets{ss, ii} = makePacket(params);
        packets{ss, ii}.response.metaData.lowFreqComponent=params.lowFreqComponentValues;
    end
    fprintf('\n');
end

%% Merge sessions
NSessionsMerged = length(whichSessionsToMerge);
for mm = 1:NSessionsMerged
    mergeIdx = whichSessionsToMerge{mm};
    mergedPacket = {packets{mergeIdx, :}};
    mergedPacket = mergedPacket(~cellfun('isempty', mergedPacket));
    mergedPacketCellArray{mm} = mergedPacket;
end


end

