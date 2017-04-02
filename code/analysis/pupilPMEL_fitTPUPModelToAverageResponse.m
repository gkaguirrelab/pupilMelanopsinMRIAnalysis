function [twoComponentFitToData] = pupilPMEL_fitTPUPModelToAverageResponse(mergedPacketCellArray, dropboxAnalysisDir)

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, 'lowFreqClean', true, 'aggregateMethod', 'median' );

%% SETUP

% Identify the number of (sessions x subjects) and the number of stimulus
% types from the dimensionality of the avgPackets cell array
NSessionsMerged=size(avgPackets,1);
NStimTypes=size(avgPackets,2);

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
    for mm = 1:NStimTypes-1 % Don't try to fit the attention task
        % Update the user
        fprintf('* Subject, stimulus <strong>%g</strong> , <strong>%g</strong>', ss, mm);
        fprintf('\n');
        
        % Construct the model object
        temporalFit = tfeTPUP('verbosity','none');
        
        % Grab a single packet
        singlePacket=avgPackets{ss, mm};
        
        % Scale the responses by 100 so that the fit params have roughly
        % equal weight in the minimization. This places the values in %
        % change units
        singlePacket.response.values=singlePacket.response.values*100;
        
        % Use the sem of the response to create a errorWeightVector
        % Increase the importance of fitting the first five seconds of data
        errorWeightVector=1./(100*singlePacket.response.sem);
        errorWeightVector(1:5000)=errorWeightVector(1:5000)*5;
        
        % Conduct the fit
        [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, ...
            'defaultParamsInfo', defaultParamsInfo, ...
            'errorWeightVector', errorWeightVector);
        
        % Store the time-series for each of the three model components
        tmpParams=paramsFit;
        tmpParams.paramMainMatrix(5:6)=0;
        tmp = temporalFit.computeResponse(tmpParams,singlePacket.stimulus,[]);
        modelResponseStruct.component1=tmp.values;
        tmpParams=paramsFit;
        tmpParams.paramMainMatrix([4 6])=0;
        tmp = temporalFit.computeResponse(tmpParams,singlePacket.stimulus,[]);
        modelResponseStruct.component2=tmp.values;
        tmpParams=paramsFit;
        tmpParams.paramMainMatrix(4:5)=0;
        tmp = temporalFit.computeResponse(tmpParams,singlePacket.stimulus,[]);
        modelResponseStruct.component3=tmp.values;
        
        paramsFit.paramMainMatrix
        
        % Store the fitResponse
        twoComponentFitToData{ss,mm}.paramsFit=paramsFit;
        twoComponentFitToData{ss,mm}.fVal=fVal;
        twoComponentFitToData{ss,mm}.modelResponseStruct=modelResponseStruct;
        
        % Clear the object
        delete(temporalFit);
        
    end % loop over subjects
end % loop over stimuli


end % function