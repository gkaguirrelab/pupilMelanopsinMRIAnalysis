function [theResult] = fitIAMPModelToIndividualResponse(accumStimTypesResponse,accumStimTypesStimulus, avgPackets)

% Loop across stimulus types (contrast level)

% Build a single packet container
% the appropriate stim type
% the kernel (which is the average response across trials for that stim
% type, available in average packets)

% Loop over stimulus instances
% place the response into the packet
% perform the model fit
% record the amplitude
% record the intial value of the response

% Loop structure
%  sessions
%    stimTypes (contrast levels plus attention task)
%       [make a packet that has a given stimulus profile and kernel]
%        instances / events







% wanted to have a modified accumStimTypesResp that had data in it's native
% units (diameter, mm) and that began each trial at 0
for ss = 1:length(mergedPackets)
    % update process
    fprintf('>> Processing <strong>%s</strong> | <strong>%s</strong> | <strong>%s</strong>\n', params.sessionType, params.sessionObserver, params.sessionDate);
    % create average run from all contrasts levels to be used as a model
    
    
    % create matrices to store results
    accumStimTypesRespRaw{ss, 1} = [];
    accumStimTypesRespRaw{ss, 2} = [];
    accumStimTypesRespRaw{ss, 3} = [];
    accumStimTypesRespRaw{ss, 4} = [];
    accumStimTypesRespRaw{ss, 5} = [];
    accumStimTypesRespRaw{ss, 6} = [];
    
    
    NRunsTotal = length(mergedPackets{ss});
    
    % make a model for each contrast level. avgPackets doesn't quite get
    % what we want for two reasons: it's shape will be based on percent
    % amplitude change, and a given run's starting point is not 0
    
    
    
    % loop over runs
    for ii = 1:NRunsTotal
        fprintf('\t* Run <strong>%g</strong> / <strong>%g</strong>\n', ii, NRunsTotal);
        % Find the stimulus onsets so that we can align the data to it. We
        % do that by finding a [0 1] edge from a difference operator.
        tmp = sum(mergedPackets{ss}{ii}.stimulus.values);
        tmp2 = diff(tmp);
        tmp2(tmp2 < 0) = 0;
        tmp2(tmp2 > 0) = 1;
        stimOnsets = strfind(tmp2, [0 1]);
        
        % Get the number of segments from the stimulus onsets
        NSegments = length(stimOnsets);
        
        t = (0:extractionDurInd)/1000; % tells how long each stimulus should be, based on information specified in MelanopsinMR_PupilAnalysis
        
        % loop over individual stimuli
        for jj = 1:length(stimOnsets)
            
            % figure out temporal boundaries of pupil response, adjusting
            %  if stimulus happens near the end of the run
            if (stimOnsets(jj)+extractionDurInd) <= length(mergedPackets{ss}{ii}.response.values);
                idxToExtract = stimOnsets(jj):(stimOnsets(jj)+extractionDurInd);
            else
                idxToExtract = stimOnsets(jj):length(mergedPackets{ss}{ii}.response.values);
            end
            thisPacket.packetType = 'pupil';
            thisPacket.stimulusFile = mergedPackets{ss}{ii}.metaData.stimulusFile;
            thisPacket.responseFile = mergedPackets{ss}{ii}.metaData.responseFile;
            
            % based on timing information, actually grab pupil response
            thisPacket.respValues = mergedPackets{ss}{ii}.response.values(idxToExtract);
            
            % for storing data, identify which type of contrast was
            % associated with the stimulus
            contrast = mergedPackets{ss}{ii}.stimulus.metaData.stimTypes(jj);
            
            % make new dataBin for data that has following properties: data
            % from raw data in pupil diameter (mm) (not percent change),
            % and each response starts from 0
            
            
            
            % this would be the main line to change in the
            % pupilMelanopsinMRIAnalysis if we wanted this to stick
            accumStimTypesRespRaw{ss,contrast} = [accumStimTypesRespRaw{ss,contrast}; (thisPacket.respValues - mergedPackets{ss}{ii}.response.values(stimOnsets(jj)))];
            
        end % end loop over stimuli per run
    end % end loop over runs
end % end loop over subjects/sessions


% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];
IAMPFitToData = []

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;

% Construct the model object
temporalFit = tfeIAMP('verbosity','high');

% Announce what we are about to do
fprintf('>> Fitting individual amplitude model to data (IAMP)\n');


for ss = 1:size(accumStimTypesRespRaw,1) % Looping over subjects
    betaPerAmplitude{ss,1} = []
            betaPerAmplitude{ss,2} = []
            betaPerAmplitude{ss,3} = []
            betaPerAmplitude{ss,4} = []
            betaPerAmplitude{ss,5} = []
    for cc = 1:(size(accumStimTypesRespRaw,2)-1) % Looping over stimulus types, skipping attention tasks
        % Update the user
        fprintf('* Subject, Contrast <strong>%g</strong> of <strong>%g</strong>', ss, 5);
        fprintf('\n');
        
        % allocate room for results
       
        
        % Update the packet to have current kernel and stimulus profile
        % Note that the kernel should be scaled to have a unit excursion
        
        singlePacket.stimulus.values = zeros(1,13000)% blip to be convolved with kernel; fixed per subject per contrast
        singlePacket.stimulus.values(1,1) = 1
        singlePacket.stimulus.timebase = [0:12999]
        singlePacket.kernel.values = nanmean(accumStimTypesRespRaw{ss,cc})
        singlePacket.kernel.values = singlePacket.kernel.values/(abs(min(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
        singlePacket.response.values = [] % to be filled in for each individual trial
        singlePacket.response.timebase = [0:12999]
        singlePacket.metaData = mergedPackets{1,ss}{1,1}.metaData % steals from one of the runs of the same subject
        
        
        % Loop over instances / events
        for ii = 1:size(accumStimTypesRespRaw{ss,cc},1)
            
            % Update the packet to have the response values for this event
            singlePacket.response.values= accumStimTypesRespRaw{ss,cc}(ii,:)
            
            % Conduct the fit
            [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix);
            
            % Store the fitResponse
            % Save the paramsFit. You'll discover that this is a structure,
            % with one field being "amplitude". Save the value associated with
            % that.
            IAMPFitToData{ss,cc}{1,ii}.paramsFit=paramsFit;
            IAMPFitToData{ss,cc}{1,ii}.fVal=fVal;
            IAMPFitToData{ss,cc}{1,ii}.modelResponseStruct=modelResponseStruct;
            betaPerAmplitude{ss,cc} = [betaPerAmplitude{ss,cc} paramsFit.paramMainMatrix];
          
            
        end % loop over events
        
    end % loop over stimulus types
end % loop over sessions

% Clear the object
delete(temporalFit);

end % function