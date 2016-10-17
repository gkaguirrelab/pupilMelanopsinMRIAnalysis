function [twoComponentFitToData] = fitTPUPModelToAverageResponse(mergedPacketCellArray, dropboxAnalysisDir)

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag );

%% SETUP

% Identify the number of (sessions x subjects) and the number of stimulus
% types from the dimensionality of the avgPackets cell array

NSessionsMerged=size(avgPackets,1);
NStimTypes=size(avgPackets,2);

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;

% Announce what we are about to do
fprintf('>> Fitting two-component model to pupil data (TPUP)\n');

% Pre-allocate a structure to hold the fit results
twoComponentFitToData{NSessionsMerged,NStimTypes}.paramsFit=[];
twoComponentFitToData{NSessionsMerged,NStimTypes}.fVal=[];
twoComponentFitToData{NSessionsMerged,NStimTypes}.fitResponse=[];

% Loop over subjects and stimulus types
% Skipping the attention task for now
for ss = 1:NSessionsMerged
    for mm = 1:NStimTypes-1
        % Update the user
        fprintf('* Subject, stimulus <strong>%g</strong> , <strong>%g</strong>', ss, mm);
        fprintf('\n');
        
        % Construct the model object
        temporalFit = tfeTPUP('verbosity','none');
        
        % Grab a single packet
        singlePacket=avgPackets{ss, mm};
        
        % Conduct the fit
        [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo, ...
            'paramLockMatrix',paramLockMatrix);
        
        % Store the fitResponse
        twoComponentFitToData{ss,mm}.paramsFit=paramsFit;
        twoComponentFitToData{ss,mm}.fVal=fVal;
        twoComponentFitToData{ss,mm}.modelResponseStruct=modelResponseStruct;
        
        % Clear the object
        delete(temporalFit);
        
    end % loop over subjects
end % loop over stimuli


end % function