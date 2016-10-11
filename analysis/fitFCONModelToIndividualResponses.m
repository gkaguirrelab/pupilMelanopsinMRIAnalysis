function [theResult] = fitFCONModelToIndividualResponses(mergedPacketCellArray, twoComponentFitToData)
%
% Description

% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=3; % zero center the initial period, % change units
splitParams.normalizationWindowMsecs=100; % define the size of the norm window

% instantiate the TPUP model
tpupModel=tfeTPUP('verbosity','none');

% instantiate the FCON model
fconModel=tfeFCON('verbosity','high');
defaultParamsInfo.nInstances=1;

% Loop over sessions
nSessions=size(mergedPacketCellArray,2);
for ss=1:nSessions
    
    % Build an fcon strcuture that is to be added to the stimulus struct
    % This contains the lookup table that maps effective contrast to the
    % expanded set of parameters for the TPUP model.
    
    contrastbase=[25,50,100,200,400];
    nContrasts=5;
    nParams=7; % number of params in the TPUP model
    observedParamMatrix=zeros(nParams,nContrasts);
    
    for cc=1:nContrasts
        observedParamMatrix(:,cc)=twoComponentFitToData{ss,cc}.paramsFit.paramMainMatrix;
    end
    
    % Spline interpolate the observedParamMatrix, following log spacing
    interpContrastBase=logspace(log10(25/2),log10(400*2),70);
    for pp=1:nParams
        interpParamMatrix(pp,:)=spline(log10(contrastbase),observedParamMatrix(pp,:),log10(interpContrastBase));
    end

    % Assemble the fcon structure for this session / subject
    fcon.contrastbase=interpContrastBase;
    fcon.observedParamMatrix=interpParamMatrix;
    fcon.modelObjHandle=tpupModel;
    fcon.logContrastFlag=false;
    
    
    % Loop over runs
    nRuns=size(mergedPacketCellArray{1,ss});
    for rr=1:nRuns
        
        % Get the packet for this run
        theRunPacket=mergedPacketCellArray{1,ss}{rr};
        
        % Loop over individual instances
        nInstances=size(theRunPacket.stimulus.values,1);
        for ii=1:5%nInstances
            
            % update the splitParams with the instance index
            splitParams.instanceIndex = ii;
            
            % grab the packet for this particular instance
            theInstancePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
            
            % make sure the kernel field is empty
            theInstancePacket.kernel=[];
            
            % add the fcon structure to the stimulus structure
            theInstancePacket.stimulus.fcon=fcon;
            
            
            % perform the fit
            [paramsFit,fVal,modelResponseStruct] = ...
                fconModel.fitResponse(theInstancePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',[]);
            
            % store the effective contrast from the paramsFit in an appropriate
            % variable
            
            myResultsVariable(ss,rr,ii)=paramsFit.paramMainMatrix; % this should be a scalar
            
            
        end % loop over instances
    end % loop over runs
end % loop over sessions


