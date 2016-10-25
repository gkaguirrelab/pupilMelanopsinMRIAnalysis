function [simulatedMergedPacketCellArray] = pupilPMEL_makeSimulatedData()

% Specify how the low frequency component is created. Two options: (1) low 
% frequency component relates to arousal, and this low frequency component
% has no carryover effect. (2) low frequency component relates to a
% carryover effect, where seeing a 400% contrast leaves the pupil a little
% bit more constricted

%lowFreq = 'arousal';
lowFreq = 'carryOver'



%% Load template data, in this case the cached version of our main data set

% Here just stealing the code to load the packetCache to make it easy
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


%% Actually make the simulated data


% make space for final simulated data
simulatedMergedPacketCellArray = [];
for rr = 1:(size(mergedPacketCellArray{1},2)) % loop over runs
    % the final simulated data will be based off all runs for session 1
    thePacket= mergedPacketCellArray{1}{rr};
    thePacket.response.values = [];
    
    % code used for translating stimulus type into contrast
    contrastList = {0.25; 0.50; 1; 2; 4; 8};
    
    % make basic structures for stimuli and response vectors 
    thePacket.stimulus.vector = zeros(1,length(thePacket.stimulus.timebase));
    retinalResponse = zeros(1,length(thePacket.stimulus.timebase));
    neuronalResponse = zeros(1,length(thePacket.stimulus.timebase));
    
    
    % create pupil response function, which is to convolved with the
    % neuronalResponse. For now, this pupil response function is just a step to
    % 2 for 3 s. With a slight delay introduced so when later functions
    % norm the time series for a given trial we can still do so
    % appropriately.
    pupilResponseFunction = zeros(1,length(thePacket.stimulus.timebase));
    pupilResponseFunction(101:3000) = -2;
    
    
    % pre-allocate space for data
    simulatedMergedPacketCellArray{1}{rr} = [];
    
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
        
        % Make stimulus vector. In this case, just grab an existing stimulus
        % vector
        % from some random subject. contrastList converts from 1-6 stimuli
        % designation to 0.25-0.5-1-etc.
        thePacket.stimulus.vector(stimOnset) = contrastList{contrast};
        
        % Make neuronalResponse vector. Simply encode stimulus contrast on
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
        
        if strcmp(lowFreq, 'arousal');
            if ii == 1;
                pre = 5;
            else
                lowFreqComponent = -1 + (1+1).*rand(1,1);
                pre = post+lowFreqComponent;
            end
        end
        
        
        if strcmp(lowFreq, 'carryOver');
            if ii == 1;
                pre = 5;
            else
                pre = post;
            end
        end
        
         % an attempt to introduce "circularity"
        retinalResponse(stimOnset) = log10(thePacket.stimulus.vector(stimOnset))+1 * pre/10;
        
        neuronalResponse(stimOnset) = retinalResponse(stimOnset)/10 * pre;

        % make the pupilMotor response. The input for pupilMotor is
        % the neuronalResponse vector and  "pre" pupil size and the output
        % is a"stimulated" pupil size (how small it gets in response to the trial) and
        % "post" pupil size. Post pupil size is an area in which we
        % can add some of this memory effect -> like if in response to a
        % 400% contrast, the pupil stays a little extra dilated. Post pupil
        % response is also an area in which we can add a low frequency
        % component. Note that we're starting with a baseline size of 5 mm.
        
        % for now, when shown a stimulus, the pupil contrast by 2 times the
        % neuronalResponse
        stimulated = -2 * neuronalResponse(stimOnset);
        
        % pupilMotor functions differently depending on which lowFreq
        % component is specified
        if strcmp(lowFreq, 'arousal');
            post = pre;
        end
        
        
        if strcmp(lowFreq, 'carryOver');
            % stay a little bit more constricted after 400% and 200%
            % contrast, a little more dilated after 25% and 50% contrast.
            % no change after 100% contrast
            post = pre  + - pre/2*log10(thePacket.stimulus.vector(stimOnset));
        end
        
        % turn pre, stimulated, and post into a reasonable response vector.
        % The response vector always comes back to the same baseline size
        % as before the stimulus so model fitting isn't ambiguous
        %
        thePacket.response.values(stimOnset:stimOnset+101)=pre;
        thePacket.response.values(stimOnset+101:stimOnset+3000) = stimulated;
        thePacket.response.values(stimOnset+3000:length(thePacket.stimulus.timebase)) = pre;
        
        
        
    end
    
  
    
    % save in proper format
    simulatedMergedPacketCellArray{1}{rr} = thePacket;
    simulatedMergedPacketCellArray{1}{rr}.metaData.simulationStyle = lowFreq;
end

end

