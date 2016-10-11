function [theResult] = fitFCONModelToIndividualResponses(mergedPackets, twoComponentFitToData)
%
% Description

contrastbase=[25,50,100,200,400];

normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;

% instantiate the FCON model
fconModel=tfeFCON('verbosity','high');
defaultParamsInfo.nInstances=1;

% Loop over sessions
for ss=1:1 %nSessions
    
    
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
    
    % We could write some code here to expand the observed parameter matrix
    % along an interpolated contrast base
    
    % Instantiate the TPUP model object to be passed as the forward model in
    % the FCON implementation
    tpupModel=tfeTPUP('verbosity','none');
    
    % Assemble the fcon structure for this session / subject
    fcon.contrastbase=contrastbase;
    fcon.observedParamMatrix=observedParamMatrix;
    fcon.modelObjHandle=tpupModel;
    fcon.logContrastFlag=false;
    
    
    % Loop over individual trials (for a given subject / session)
    for ii=1:nInstances
        
        % grab the packet for this particular instance
        thePacket=mergedPackets(this session, this run, this instance)
        % keep the kernel field empty
        % stimulus structure needs the actual stimulus specification, not an
        % impulse
        
        % add the fcon structure to the stimulus structure
        thePacket.stimulus.fcon=fcon;
        
        % perform the fit
        [paramsFit,fVal,modelResponseStruct] = fconModel.fitResponse(thePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',[]);
        
        % store the effective contrast from the paramsFit in an appropriate
        % variable
        
        myResultsVariable(session,run,subject)=paramsFit.paramMainMatrix; % this should be a scalar
        
        
    end % loop over instances
    
    
end % loop over sessions

end % function