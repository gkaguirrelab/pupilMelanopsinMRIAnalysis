%% The code below is a direct copy-and-paste from pupilPMEL_makeSimulatedData.m 
% I've noticed this error across a couple different versions of my
% simulated data, but this specific simulation has the property that I know
% which instance is affected.

lowFreq = 'carryOver'
%lowFreq = 'attention'

% Here just stealing the code to load the packetCache to make it easy
% Discover user name and set Dropbox path
[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxAnalysisDir = ...
    fullfile('/Users', userName, ...
    '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/pupilMelanopsinMRIAnalysis');
packetCacheBehavior='load';
packetCellArrayTag='maxMelLSM_CRF_Pupil';
packetCellArrayHash='d702e5898732600707fb32a2302d4772';
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
% make space for final simulated data
simulatedMergedPacketCellArray = [];
for rr = 1:(size(mergedPacketCellArray{1},2)) % loop over runs
    % the final simulated data will be based off all runs for session 1
    thePacket= mergedPacketCellArray{1}{rr};
    thePacket.response.values = 5*ones(1,length(thePacket.stimulus.timebase));
    
    % code used for translating stimulus type into contrast
    contrastList = {0.25; 0.50; 1; 2; 4; 1};
    
    % make basic structures for stimuli and response vectors 
    thePacket.stimulus.vector = zeros(1,length(thePacket.stimulus.timebase));
    retinalResponse = zeros(1,length(thePacket.stimulus.timebase));
    neuronalResponse = zeros(1,length(thePacket.stimulus.timebase));
    

    
    
    % pre-allocate space for data
    simulatedMergedPacketCellArray{1}{rr} = [];
    stimOnsetList = [];
    for ii = 1:size(thePacket.stimulus.values,1); % loop over instances
        
        % Find the stimulus onsets so that we can align the data to it. We
        % do that by finding a [0 1] edge from a difference operator.
        tmp = diff(thePacket.stimulus.values(ii,:));
        tmp(tmp < 0) = 0;
        tmp(tmp > 0) = 1;
        
        % Check if the very first value is 1, in which case the stim onset is
        % at the initial value
        if tmp(1)==1
            stimOnset=1;
        else
            stimOnset = strfind(tmp, [0 1]);
        end
        contrast = thePacket.stimulus.metaData.stimTypes(ii);
        stimOnsetList = [stimOnsetList, stimOnset];
        
        
        % Make stimulus vector. In this case, just grab an existing stimulus
        % vector
        % from some random subject. contrastList converts from 1-6 stimuli
        % designation to 0.25-0.5-1-etc.
        thePacket.stimulus.vector(stimOnset) = contrastList{contrast};
        
        % Make retinalResponse vector. Simply encode stimulus contrast on
        % a log scale
        retinalResponse(stimOnset) = log10(thePacket.stimulus.vector(stimOnset))+1;
        
        
        
        % In this mode, neuronal response is proportional to the size at the
        % beginning of stimulus presentation. That is neuronalResponse is a
        % percent change related to baseline size, where the percent is
        % dictated by the stimulus type. Pupil size can determined by output
        % of pupilMotor below (if there is a carryover effect and we don't 
        % simply return back to baseline) or an independent low frequency 
        % component, but for the first stimulus, we're starting
        % with a baseline size of 5mm
        
        % if the low frequency component is attention, this corresponds to
        % a model in which attention can manipulate both baseline pupil
        % size and percentage change (outputed by the neuronalResponse
        % function). More specifically, higher attention states are
        % associated with increased baseline pupil size and increased
        % percentage change
        if strcmp(lowFreq, 'attention');
            if ii == 1;
                pre = 5;
                neuronalResponseChangeFactor = 1;
                attentionComponent = 0;
            else
                neuronalResponseChangeFactor = 1;
                %attentionComponent = randi([0,1])*2 - 1; % randomly flip between high and low attentional states between trials
                attentionComponent = sin(stimOnset*2*pi/44800); % baseline size fluctuates slowly between high and low attentional states
                pre = 5 + attentionComponent;
            end
        end
        
        if strcmp(lowFreq, 'arousal');
            if ii == 1;
                pre = 5;
                neuronalResponseChangeFactor = 1;
                attentionComponent = 0;
            else
                lowFreqComponent = -1 + (1+1).*rand(1,1);
                pre = delta+lowFreqComponent;
                neuronalResponseChangeFactor = 1;
                attentionComponent = 0;
            end
        end
        
        
        if strcmp(lowFreq, 'carryOver');
            if ii == 1;
                pre = 5;
                neuronalResponseChangeFactor = 1;
                attentionComponent = 0;
            else
                pre = delta;
                attentionComponent = 0;
            end
        end
        
         % One method of giving our data have a positive relationship
         % between baseline size and percent signal change is to have the
         % reintal response also be sensitive to light
        %retinalResponse(stimOnset) = (log10(thePacket.stimulus.vector(stimOnset))+1) * pre/10;
        
        % added a carry-over effect, where the neuronalResponse is changed
        % by the stimulus that came before. neuronalResponseChangeFactor is
        % defined below so that it applies to the stimuli that comes next.
        neuronalResponse(stimOnset) = retinalResponse(stimOnset)/10 * 1 / neuronalResponseChangeFactor * (attentionComponent/5+1);

        % make the pupilMotor response. The input for pupilMotor is
        % the neuronalResponse vector and  "pre" pupil size and the output
        % is a"stimulated" pupil size (how small it gets in response to the trial) and
        % "change" pupil size, the corresponds to a different starting 
        % position for the next run. Post pupil size is an area in which we
        % can add some of this memory effect -> like if in response to a
        % 400% contrast, the pupil stays a little extra dilated. Post pupil
        % response is also an area in which we can add a low frequency
        % component. Note that we're starting with a baseline size of 5 mm.
        
        % for now, when shown a stimulus, the pupil contrast by 2 times the
        % neuronalResponse
        stimulated = (-2 * neuronalResponse(stimOnset)*pre) + pre;
        
        % pupilMotor functions differently depending on which lowFreq
        % component is specified
        if strcmp(lowFreq, 'attention');
            delta = pre;
        end
        if strcmp(lowFreq, 'arousal');
            delta = pre;
        end
        
        
        if strcmp(lowFreq, 'carryOver');
            % stay a little bit more constricted after 400% and 200%
            % contrast, a little more dilated after 25% and 50% contrast.
            % no change after 100% contrast
            
            % to make baseline size totally determined by previous contrast
            % level
            delta = 5- 2*log10(thePacket.stimulus.vector(stimOnset));
            
            % to make baseline size move from where it's at to a new
            % location based on previous contrast level. note that this
            % means that baseline size is predicted less well by previous
            % contrast (since in this case, previous contrast is just one
            % of the long history of stimuli that have shaped the baseline
            % size)
            %delta = pre  - 2*log10(thePacket.stimulus.vector(stimOnset));
            
            
            neuronalResponseChangeFactor = ((log10(thePacket.stimulus.vector(stimOnset)))/5) + 1;
        end
        
        % turn pre, stimulated, and post into a reasonable response vector.
        % The response vector always comes back to the same baseline size
        % as before the stimulus so model fitting isn't ambiguous
        %
        thePacket.response.values(stimOnset:stimOnset+105)=pre;
        thePacket.response.values(stimOnset+105:stimOnset+3000) = stimulated;
        thePacket.response.values(stimOnset+3000:length(thePacket.stimulus.timebase)) = pre;
        
        
        
    end
    
  
    
    % save in proper format
    simulatedMergedPacketCellArray{1}{rr} = thePacket;
    simulatedMergedPacketCellArray{1}{rr}.metaData.simulationStyle = lowFreq;
end


%% Now we actually have the representative simulation data
%subject 1, run 4, instance 3 is one of the instances affected.

ss = 1;
rr = 4;
ii = 3;


%% Specific parameters for the fitting process, taken directly from fitIAMPModelToIndividualResponse.m
% apply the relevant parameters, directly from fitIAMP
splitParams.instanceIndex=ii; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normalizationWindowMsecs=100; % define the size of the norm window
splitParams.normFlag=3; % percent change
% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];
IAMPFitToData = [];
% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;
temporalFit = tfeIAMP('verbosity','none');
normFlag = 3;




%grab the avgPacket which serves as the model
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( simulatedMergedPacketCellArray, normFlag, 'LowFreqClean', false, 'aggregateMethod', 'median');


% assemble the instance packet
theRunPacket=simulatedMergedPacketCellArray{1,ss}{rr};
singlePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
contrast = theRunPacket.stimulus.metaData.stimTypes(ii);
% create stimulus profile -> has to be a blip with this
% configuration of IAMP (it convolves the stimulus profile
% with the kernel)
singlePacket.stimulus.values = zeros(1,13000);  % blip to be convolved with kernel; fixed per subject per contrast
singlePacket.stimulus.values(1,1) = 1;
% create kernel, which is the average run for a contrast
% level per subject
singlePacket.kernel.values = avgPackets{ss,contrast}.response.values;
singlePacket.kernel.values = singlePacket.kernel.values/(abs(min(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
if strcmp(contrast, 6);
    singlePacket.kernel.values = singlePacket.kernel.values/(abs(max(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
end
singlePacket.kernel.timebase = [0:12999];

% do the fit, first without linearRgression flag
[paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix);

% plot to show what's wrong
plot(singlePacket.response.values)
hold on
plot(modelResponseStruct.values)
title('Without LinearRegression Flag')
hold off

% now do the fit with linearRegression flag on
[paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix, 'searchMethod','linearRegression');

% plot to show how things are changed
figure
plot(singlePacket.response.values)
hold on
plot(modelResponseStruct.values)
title('With LinearRegression Flag')
hold off
