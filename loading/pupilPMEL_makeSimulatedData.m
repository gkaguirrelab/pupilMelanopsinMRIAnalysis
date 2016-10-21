function [simulatedMergedPacketCellArray] = pupilPMEL_makeSimulatedData()


% Make stimulus vector. In this case, just grab an existing stimulus
% vector
% from some random subject.

% Also make the neuronalResponse vector, which is log10(stimulus)+1 (+1 so
% that 100% contrast and below gives positive signal)


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

simulatedMergedPacketCellArray = [];
for rr = 1:(size(mergedPacketCellArray{1},2))
    thePacket= mergedPacketCellArray{1}{rr};
    
    
    contrastList = {0.25; 0.50; 1; 2; 4; 8};
    thePacket.stimulus.vector = zeros(1,length(thePacket.stimulus.timebase));
    neuronalResponse = zeros(1,length(thePacket.stimulus.timebase));
    pupilResponse = zeros(1,length(thePacket.stimulus.timebase));
    
    % create pupil response function, which is to convolved with the
    % neuronalResponse. For now, this pupil response function is just a step to
    % 2 for 3 s
    pupilResponseFunction = zeros(1,length(thePacket.stimulus.timebase));
    pupilResponseFunction(101:3000) = -2;
    
    % Create pupil size vector. For now, just a baseline size of 5 mm
    pupilSize = thePacket.response.metaData.lowFreqComponent;
    
    
    simulatedMergedPacketCellArray{1}{rr} = [];
    for ii = 1:size(thePacket.stimulus.values,1);
        
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
        
        thePacket.stimulus.vector(stimOnset) = contrastList{contrast};
        neuronalResponse(stimOnset) = log10(thePacket.stimulus.vector(stimOnset))+1;
        
        % In this mode, pupil response is proportional to the size at the
        % beginning of stimulus presentation
        pupilResponse(stimOnset) = neuronalResponse(stimOnset)/10 * pupilSize(stimOnset);
        
        
    end
    
    pupilMotor = conv(pupilResponse,pupilResponseFunction);
    pupilMotor = pupilMotor(1:448000);
    thePacket.response.values = pupilSize + pupilMotor;
    
    simulatedMergedPacketCellArray{1}{rr} = thePacket;
end

end

