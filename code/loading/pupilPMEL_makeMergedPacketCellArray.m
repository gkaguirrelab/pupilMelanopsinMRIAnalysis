function [ mergedPacketCellArray, params ] = pupilPMEL_makeMergedPacketCellArray(ppsDataPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define the dropbox path
baseDir = fullfile(ppsDataPath);

% Define the session directories
sessDirs = {...
    'MaxLMSCRF/HERO_asb1/060816' ...
    'MaxLMSCRF/HERO_aso1/060116' ...
    'MaxLMSCRF/HERO_gka1/060616' ...
    'MaxLMSCRF/HERO_mxs1/062816' ...
    'MaxMelCRF/HERO_asb1/060716' ...
    'MaxMelCRF/HERO_aso1/053116' ...
    'MaxMelCRF/HERO_gka1/060216' ...
    'MaxMelCRF/HERO_mxs1/060916' ...
    'MaxMelCRF/HERO_mxs1/061016'
    };

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
    params.sessionType = ['MelanopsinMR' tmp{2}];
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
    params.sessionDir       = fullfile(baseDir, sessDirs{ss});
    NRuns = length(listdir(fullfile(params.sessionDir, 'MatFiles', '*.mat'), 'files'));
    
    % Iterate over runs
    for ii = 1:NRuns
        fprintf('\t* Run <strong>%g</strong> / <strong>%g</strong>\n', ii, NRuns);
        % Set up some parameters
        params.runNum           = ii;
        params.stimulusFile     = fullfile(params.sessionDir, 'MatFiles', [params.sessionObserver '-' params.sessionType '-' num2str(ii, '%02.f') '.mat']);
        params.responseFile     = fullfile(params.sessionDir, 'EyeTrackingFiles', [params.sessionObserver '-' params.sessionType '-' num2str(ii, '%02.f') '.mat']);
        params.lowFreqFilterHz  = 0.025; % The boundary of the low-pass filter that will be used to estimate the component of slow pupil modulation in the signal
        
        % new line here
        calibrationFileTemp  = cellstr(ls(fullfile(params.sessionDir, 'EyeTrackingFiles/ScaleCalibration', '*.mat')));
        params.calibrationFile = calibrationFileTemp{1,1};
        [params.respValues params.respTimeBase params.gazeX params.gazeY params.blinks params.lowFreqComponentValues] = pupilPMEL_loadPupilDataForPackets(params);
        [params.stimValues params.stimTimeBase params.stimMetaData] = pupilPMEL_makeStimStruct(params);
        packets{ss, ii} = makePupilPacket(params);
        packets{ss, ii}.response.metaData.lowFreqComponent=params.lowFreqComponentValues;
        packets{ss, ii}.response.gazeX=params.gazeX;
        packets{ss, ii}.response.gazeY=params.gazeY;
        packets{ss, ii}.response.blinks=params.blinks;
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


end % function


function packet = makePupilPacket(params)


%   Outputs a 'packet' structure with stimulus, response, metaData, and
%   (optionally) HRF information
%
%   Usage:
%   packet = makePacket(params)
%
%   params structure:
%   params.packetType       - 'bold' or 'pupil'
%   params.sessionDir       - session directory, full path
%   params.stimulusFile     - full path to stimulus file
%   params.stimValues       - 1 x N vector of stimulus values
%   params.stimTimeBase     - 1 x N vector of stimulus times (msec)
%   params.stimMetaData     - Any stimulus meta data
%   params.responseFile     - full path to response file
%   params.respValues       - 1 x N vector of response values
%   params.respTimeBase     - 1 x N vector of response times (msec)
%
%   If strcmp(params.packetType,'bold')
%
%   params.hrfFile          - full path to HRF file
%
%   Output fields in packets:
%
%   stimulus.values         - M x N matrix modeling M stimulus events
%   stimulus.timebase       - 1 x N vector of stimulus times (msec)
%   stimulus.metaData       - structure with info about the stimulus
%
%   response.values         - 1 x N vector of response values
%   response.timebase       - 1 x N vector of response times (msec)
%   response.metaData       - structure with info about the response
%
%   metaData.projectName    - project name (e.g. 'MelanopsinMR');
%   metaData.subjectName    - subject name (e.g. 'HERO_asb1');
%   metaData.sessionDate    - session date (e.g. '041416');
%   metaData.stimulusFile   - fullfile(sessionDir,'MatFiles',matFiles{i});
%   metaData.responseFile   - fullfile(sessionDir,boldDirs{i},[func '.nii.gz']);
%
%   If packetType == 'bold', also outputs:
%
%   kernel.values           - 1 x N vector of response values
%   kernel.timebase         - 1 x N vector of response times (msec)
%   kernel.metaData         - structure with info about the HRF
%
%   Otherwise, the fields in kernel are the empty matrix ([]).
%
%   Written by Andrew S Bock Aug 2016

%% Metadata
[subjectStr,sessionDate]            = fileparts(params.sessionDir);
[projectStr,subjectName]            = fileparts(subjectStr);
[~,projectName]                     = fileparts(projectStr);
metaData.projectName                = projectName;
metaData.subjectName                = subjectName;
metaData.sessionDate                = sessionDate;
metaData.stimulusFile               = params.stimulusFile;
metaData.responseFile               = params.responseFile;
%% Stimulus
stimulus.values                     = params.stimValues;
stimulus.timebase                   = params.stimTimeBase;
stimulus.metaData.filename          = params.stimulusFile;
stimulus.metaData                   = params.stimMetaData;
%% Response
response.values                     = params.respValues;
response.timebase                   = params.respTimeBase;
response.metaData.filename          = params.responseFile;
%% Kernel
switch params.packetType
    case 'bold'
        % HRF (if applicable)
        tmp                         = load(params.hrfFile);
        kernel.values               = tmp.HRF.mean;
        kernel.timebase             = 0:length(kernel.values)-1;
        kernel.metaData             = tmp.HRF.metaData;
    otherwise
        kernel.values               = [];
        kernel.timebase             = [];
        kernel.metaData             = [];
end
%% Save the packets
packet.stimulus                     = stimulus;
packet.response                     = response;
packet.metaData                     = metaData;
packet.kernel                       = kernel;

end
