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
defaultParamsInfo.nInstances=1;

% Set up the elements of the fcon stucture.
% This contains the lookup table that maps effective contrast to the
% expanded set of parameters for the TPUP model.
contrastbase=[0.25,0.50,1.0,2.0,4.0];
nContrasts=5;
nParams=7; % number of params in the TPUP model
observedParamMatrix=zeros(nParams,nContrasts);
interpContrastBase = ...
    logspace(log10(min(contrastbase)/2),log10(max(contrastbase)*2),(nContrasts+2)*100);
% set tpup lower bounds
vlb = [0.1, 0.1, eps, eps, eps, eps, 0.5];
vub = [0.4, 0.3, 1, 2, 1, 2, 6.0];

% instantiate the FCON model
fconModel=tfeFCON('verbosity','high');

% Loop over sessions
nSessions=size(mergedPacketCellArray,2);
for ss=1:nSessions
    
    % Build the paramLookUpMatrix for this session /subject
    for cc=1:nContrasts
        observedParamMatrix(:,cc)=twoComponentFitToData{ss,cc}.paramsFit.paramMainMatrix;
    end
    
    % Spline interpolate the observedParamMatrix 10x, following log spacing,
    % and extending one 2x spacing above and below the stimulus range
    for pp=1:nParams
        interpParams = ...
            spline(log10(contrastbase),observedParamMatrix(pp,:),log10(interpContrastBase));
        interpParams(find(interpParams < vlb(pp))) = vlb(pp);
        interpParams(find(interpParams > vub(pp))) = vub(pp);
        paramLookUpMatrix(pp,:)=interpParams;
    end
        
    % Assemble the fcon structure for this session / subject
    fcon.contrastbase=log10(interpContrastBase);
    fcon.paramLookUpMatrix=paramLookUpMatrix;
    fcon.modelObjHandle=tpupModel;
    fcon.defaultParams=tpupModel.defaultParams('defaultParamsInfo',defaultParamsInfo);

    
    % Loop over runs
    nRuns=size(mergedPacketCellArray{1,ss});
    for rr=1:nRuns
        
        % Get the packet for this run
        theRunPacket=mergedPacketCellArray{1,ss}{rr};
        
        % Loop over individual instances
        nInstances=size(theRunPacket.stimulus.values,1);
        for ii=1:5 %nInstances
            
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
                fconModel.fitResponse(theInstancePacket, 'defaultParamsInfo', defaultParamsInfo);
            
            % store the effective contrast from the paramsFit in an appropriate
            % variable
            
            myResultsVariable(ss,rr,ii)=paramsFit.paramMainMatrix; % this should be a scalar
            
            
        end % loop over instances
    end % loop over runs
end % loop over sessions


