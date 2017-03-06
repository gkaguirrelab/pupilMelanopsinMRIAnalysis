function [ mergedPacketCellArray, params ] = pupilPMEL_makeMergedPacketCellArray( userName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define the dropbox path
dropboxDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data');

% Define the session directories
sessDirs = {...
    'MelanopsinMR_fMRI/MaxLMSCRF/HERO_asb1/060816' ...
    'MelanopsinMR_fMRI/MaxLMSCRF/HERO_aso1/060116' ...
    'MelanopsinMR_fMRI/MaxLMSCRF/HERO_gka1/060616' ...
    'MelanopsinMR_fMRI/MaxLMSCRF/HERO_mxs1/062816' ...
    'MelanopsinMR_fMRI/MaxMelCRF/HERO_asb1/060716' ...
    'MelanopsinMR_fMRI/MaxMelCRF/HERO_aso1/053116' ...
    'MelanopsinMR_fMRI/MaxMelCRF/HERO_gka1/060216' ...
    'MelanopsinMR_fMRI/MaxMelCRF/HERO_mxs1/060916' ...
    'MelanopsinMR_fMRI/MaxMelCRF/HERO_mxs1/061016'
    };
sessionTypes = {'MelanopsinMRMaxLMSCRF' ...
    'MelanopsinMRMaxLMSCRF' ...
    'MelanopsinMRMaxLMSCRF' ...
    'MelanopsinMRMaxLMSCRF' ...
    'MelanopsinMRMaxMelCRF' ...
    'MelanopsinMRMaxMelCRF' ...
    'MelanopsinMRMaxMelCRF' ...
    'MelanopsinMRMaxMelCRF' ...
    'MelanopsinMRMaxMelCRF'};

% Define which sessions we'd like to merge
whichSessionsToMerge = {[1], [2], [3], [4], [5], [6], [7], [8 9]};


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
    params.sessionType = sessionTypes{ss};
    params.sessionObserver = tmp{3};
    params.sessionDate = tmp{4};
    
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
        [pupilSize params.respTimeBase] = loadEyeTrackingDataForPackets(params, 'pupil');
        [gazeX params.respTimeBase] = loadEyeTrackingDataForPackets(params, 'gazeX');
        [gazeY params.respTimeBase] = loadEyeTrackingDataForPackets(params, 'gazeY');
        params.respValues = gazeX;
        [params.stimValues params.stimTimeBase params.stimMetaData] = pupilMelanopsinMRImakeStimStruct(params);
        packets{ss, ii} = makePacket(params);
        packets{ss, ii}.response.values = [];
        packets{ss, ii}.response.gazeX = gazeX;
        packets{ss, ii}.response.gazeY = gazeY;
        packets{ss, ii}.response.pupilSize = pupilSize;
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

