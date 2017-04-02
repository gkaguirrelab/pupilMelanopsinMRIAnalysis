function [Data_LiveTrack_PupilDiameter_FineMasterTime, TimeVectorFine, Data_LiveTrack_gazeX_FineMasterTime, Data_LiveTrack_gazeY_FineMasterTime, Data_LiveTrack_blinks_FineMasterTime, Data_LiveTrack_PupilDiameter_FineMasterTime_LowFreqComponent] = pupilPMEL_loadPupilDataForPackets(params)
% values = loadPupilDataForPackets(input_dir, stimulus, params)

%% STIMULUS VALUES
Data_Stimulus = load(params.stimulusFile);

%% Calibration Factor: based on calibration runs, the conversion factor from camera units to mm width
calibration = load(params.calibrationFile);
calibrationFactor = calibration.params.cameraUnitsToMmWidthMean;

% Extract how many runs are needed
NTRsExpected = sum(Data_Stimulus.params.trialDuration)/params.TRDurSecs;

TimeVectorFine = 0:1:sum(Data_Stimulus.params.trialDuration)*1000-1;

%%% EYE TRACKING DATA %%%
% Load the eye tracking data
Data_LiveTrack = load(params.responseFile);

% Find cases in which the TTL pulse signal was split over the two
% samples, and remove the second sample.
Data_LiveTrack_TTLPulses_raw = [Data_LiveTrack.params.Report.Digital_IO1];
tmpIdx = strfind(Data_LiveTrack_TTLPulses_raw, [1 1]);
Data_LiveTrack_TTLPulses_raw(tmpIdx) = 0;

Data_LiveTrack_TTLPulses = [];
Data_LiveTrack_PupilDiameter = [];
Data_LiveTrack_gazeX = [];
Data_LiveTrack_gazeY = [];
Data_LiveTrack_IsTracked = [];
% We reconstruct the data set collected at 30/60 Hz.
for rr = 1:length(Data_LiveTrack.params.Report)
    % Depending on how we set up the acquisiton frequency, we have to
    % do different things to extract the data
    switch params.acquisitionFreq
        case 60
            Data_LiveTrack_TTLPulses = [Data_LiveTrack_TTLPulses Data_LiveTrack_TTLPulses_raw(rr) 0];
            
            % We use the pupil width as the index of pupil diameter
            Data_LiveTrack_PupilDiameter = [Data_LiveTrack_PupilDiameter Data_LiveTrack.params.Report(rr).PupilWidth_Ch01/calibrationFactor ...
                Data_LiveTrack.params.Report(rr).PupilWidth_Ch02/calibrationFactor];
            
            % Obtain the gaze position
            Data_LiveTrack_gazeX = [Data_LiveTrack_gazeX ...
                (Data_LiveTrack.params.Report(rr).PupilCameraX_Ch01 - Data_LiveTrack.params.Report(rr).Glint1CameraX_Ch01) /calibrationFactor ...
                (Data_LiveTrack.params.Report(rr).PupilCameraX_Ch02 - Data_LiveTrack.params.Report(rr).Glint1CameraX_Ch02) /calibrationFactor];
            Data_LiveTrack_gazeY = [Data_LiveTrack_gazeY ...
                (Data_LiveTrack.params.Report(rr).PupilCameraY_Ch01 - Data_LiveTrack.params.Report(rr).Glint1CameraY_Ch01) /calibrationFactor ...
                (Data_LiveTrack.params.Report(rr).PupilCameraY_Ch02 - Data_LiveTrack.params.Report(rr).Glint1CameraY_Ch02) /calibrationFactor];
                        
            % Special case
            if strcmp(params.sessionDate, '060616') && strcmp(params.sessionObserver, 'HERO_gka1')
                Data_LiveTrack_IsTracked = [Data_LiveTrack_IsTracked Data_LiveTrack.params.Report(rr).PupilTracked ...
                    Data_LiveTrack.params.Report(rr).S2Tracked];
            else
                Data_LiveTrack_IsTracked = [Data_LiveTrack_IsTracked Data_LiveTrack.params.Report(rr).PupilTracked_Ch01 ...
                    Data_LiveTrack.params.Report(rr).PupilTracked_Ch02];
            end
        case 30
            Data_LiveTrack_TTLPulses = [Data_LiveTrack_TTLPulses Data_LiveTrack_TTLPulses_raw(rr)];
            Data_LiveTrack_PupilDiameter = [Data_LiveTrack_PupilDiameter Data_LiveTrack.params.Report(rr).LeftPupilWidth/calibrationFactor];
            Data_LiveTrack_gazeX = [Data_LiveTrack_gazeX (Data_LiveTrack.params.Report(rr).LeftPupilCameraX - Data_LiveTrack.params.Report(rr).LeftGlint1CameraX) /calibrationFactor];
            Data_LiveTrack_gazeY = [Data_LiveTrack_gazeY (Data_LiveTrack.params.Report(rr).LeftPupilCameraY - Data_LiveTrack.params.Report(rr).LeftGlint1CameraY) /calibrationFactor];
            Data_LiveTrack_IsTracked = [Data_LiveTrack_IsTracked Data_LiveTrack.params.Report(rr).PupilTracked];
            
    end
end

% Now, we reconstruct the time vector of the data.
TTLPulseIndices = find(Data_LiveTrack_TTLPulses); FirstTTLPulse = TTLPulseIndices(1);
TimeVectorLinear = zeros(1, size(Data_LiveTrack_TTLPulses, 2));
TimeVectorLinear(TTLPulseIndices) = (1:NTRsExpected)-1;

% Replace zeros with NaN
tmpX = 1:length(TimeVectorLinear);
TimeVectorLinear(TimeVectorLinear == 0) = NaN;
TimeVectorLinear(isnan(TimeVectorLinear)) = interp1(tmpX(~isnan(TimeVectorLinear)), ...
    TimeVectorLinear(~isnan(TimeVectorLinear)), tmpX(isnan(TimeVectorLinear)), 'linear', 'extrap');

% Resample the timing to 1 msecs sampling
Data_LiveTrack_PupilDiameter_FineMasterTime = interp1(TimeVectorLinear*params.TRDurSecs*1000, ...
    Data_LiveTrack_PupilDiameter, TimeVectorFine);
Data_LiveTrack_gazeX_FineMasterTime = interp1(TimeVectorLinear*params.TRDurSecs*1000, ...
    Data_LiveTrack_gazeX, TimeVectorFine);
Data_LiveTrack_gazeY_FineMasterTime = interp1(TimeVectorLinear*params.TRDurSecs*1000, ...
    Data_LiveTrack_gazeY, TimeVectorFine);

% Identify when the eye was tracked
Data_LiveTrack_IsTracked_FineMasterTime = interp1(TimeVectorLinear*params.TRDurSecs*1000, ...
    Data_LiveTrack_IsTracked, TimeVectorFine, 'nearest'); % Use NN interpolation for the binary tracking state

% Extract the stimulus timing
keyPressWhich = [];
keyPressWhen = [];
for rr = 1:length(Data_Stimulus.params.responseStruct.events)
    keyPressWhich = [keyPressWhich Data_Stimulus.params.responseStruct.events(rr).buffer.keyCode];
    keyPressWhen = [keyPressWhen Data_Stimulus.params.responseStruct.events(rr).buffer.when];
end

% Extract only the ts
keyPressWhen = keyPressWhen(keyPressWhich == 18);

% Tack the first t also in this vector
stimulus.metaData_TTL = [Data_Stimulus.params.responseStruct.tBlockStart keyPressWhen];

% Subtract the absolute time of the first t
stimulus.metaData_TTL_t0 = stimulus.metaData_TTL(1);
stimulus.metaData_TTL = stimulus.metaData_TTL-stimulus.metaData_TTL_t0;

% Check that we have as many TRs as we expect
fprintf('\t> Expecting <strong>%g</strong> TRs - Found <strong>%g</strong> (LiveTrack) and <strong>%g</strong> (OneLight record).\n', ...
    NTRsExpected, sum(Data_LiveTrack_TTLPulses), length(stimulus.metaData_TTL));
if (NTRsExpected == sum(Data_LiveTrack_TTLPulses)) || (NTRsExpected == length(stimulus.metaData_TTL))
    fprintf('\t\t>> Expected number of TRs matches actual number.\n');
else
    error('\t\t>> Mismatch between expected and actual number of TRs received.');
end

% We now have four variables of interest
%   Data_LiveTrack_IsTracked_FineMasterTime <- Binary array indicating tracking state
%   Data_LiveTrack_PupilDiameter_FineMasterTime <- Pupil diameter
%   TimeVectorFine <- Time vector in TR time

% Remove blinks from the pupil data
Data_LiveTrack_BlinkIdx = [];
for rr = 1:length(params.BlinkWindowSample)
    Data_LiveTrack_BlinkIdx = [Data_LiveTrack_BlinkIdx find(~Data_LiveTrack_IsTracked_FineMasterTime)+params.BlinkWindowSample(rr)];
end
% Remove any blinks from before the first sample
Data_LiveTrack_BlinkIdx(Data_LiveTrack_BlinkIdx < 1) = []; %% ALSO CUT THE BLINKING OFF AT THE END

% Remove any blinks after the last sample
Data_LiveTrack_BlinkIdx(Data_LiveTrack_BlinkIdx > length(Data_LiveTrack_IsTracked_FineMasterTime)) = []; %% ALSO CUT THE BLINKING OFF AT THE END
Data_LiveTrack_BlinkIdx = unique(Data_LiveTrack_BlinkIdx);

% Apply the blink index to the data
Data_LiveTrack_PupilDiameter_FineMasterTime(Data_LiveTrack_BlinkIdx) = NaN;
Data_LiveTrack_gazeX_FineMasterTime(Data_LiveTrack_BlinkIdx) = NaN;
Data_LiveTrack_gazeY_FineMasterTime(Data_LiveTrack_BlinkIdx) = NaN;

% Createa a vector to report back the time of the blinks
Data_LiveTrack_blinks_FineMasterTime=Data_LiveTrack_PupilDiameter_FineMasterTime*0;
Data_LiveTrack_blinks_FineMasterTime(Data_LiveTrack_BlinkIdx)=1;

% Interpolate the elements
Data_LiveTrack_PupilDiameter_FineMasterTime(isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)) = interp1(TimeVectorFine(~isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)), Data_LiveTrack_PupilDiameter_FineMasterTime(~isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)), TimeVectorFine(isnan(Data_LiveTrack_PupilDiameter_FineMasterTime)));
Data_LiveTrack_gazeX_FineMasterTime(isnan(Data_LiveTrack_gazeX_FineMasterTime)) = interp1(TimeVectorFine(~isnan(Data_LiveTrack_gazeX_FineMasterTime)), Data_LiveTrack_gazeX_FineMasterTime(~isnan(Data_LiveTrack_gazeX_FineMasterTime)), TimeVectorFine(isnan(Data_LiveTrack_gazeX_FineMasterTime)));
Data_LiveTrack_gazeY_FineMasterTime(isnan(Data_LiveTrack_gazeY_FineMasterTime)) = interp1(TimeVectorFine(~isnan(Data_LiveTrack_gazeY_FineMasterTime)), Data_LiveTrack_gazeY_FineMasterTime(~isnan(Data_LiveTrack_gazeY_FineMasterTime)), TimeVectorFine(isnan(Data_LiveTrack_gazeY_FineMasterTime)));

% Get the DC

%%
% For now, we want to look at relationship between baseline pupil size and
% amplitude on data that has not already been filtered. The only change
% is to add a final line that defines values without manipulating by the
% filtered component

% Low-pass filter the pupil data
% Set up filter properties
nHarmonicsToFilter = round(...    % Take the closeset integer number of frequencies
                 ((max(TimeVectorFine)+1)/1000.)/ ... % Total duration of the data in seconds (assuming a deltaT of 1 msec)
                 (1/params.lowFreqFilterHz));  % the boundary low frequency

for ii = 1:nHarmonicsToFilter
    X(2*ii-1,:) = sin(linspace(0, 2*pi*ii, size(TimeVectorFine, 2)));
    X(2*ii,:) = cos(linspace(0, 2*pi*ii, size(TimeVectorFine, 2)));
end
X(nHarmonicsToFilter*2+1,:) = ones(1,size(TimeVectorFine, 2)); % create an intercept term

% Filter it
[b, bint, r] = regress(Data_LiveTrack_PupilDiameter_FineMasterTime',X');

% Obtain the low-frequency component
Data_LiveTrack_PupilDiameter_FineMasterTime_LowFreqComponent = (X'*b)';

end % function