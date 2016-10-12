function [myResultsVariable] = fitFCONModelToIndividualResponses(mergedPacketCellArray, twoComponentFitToData)
%
% Description

% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=3; % zero center the initial period, change units
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
    logspace(log10(min(contrastbase)/2/2),log10(max(contrastbase)*2*2),(nContrasts+4)*100);

% set tpup lower bounds. Should derive this from a call to tpupModel
vlb = [0.1, 0.1, eps, eps, eps, eps, 0.5];
vub = [0.4, 0.3, 1, 2, 1, 2, 6.0];

% The DiffMinChange will be used to guide fmincon later
DiffMinChange=mean(diff(log10(interpContrastBase)));

% instantiate the FCON model
fconModel=tfeFCON('verbosity','none');

% Announce what we are about to do
fprintf('>> Fitting effective contrast (FCON) / TPUP model to individual trials\n');

% Loop over sessions
nSessions=size(mergedPacketCellArray,2);
for ss=1:nSessions
    
    % Build the paramLookUpMatrix for this session /subject
    for cc=1:nContrasts
        observedParamMatrix(:,cc)=twoComponentFitToData{ss,cc}.paramsFit.paramMainMatrix;
    end
    
    % Interpolate the observedParamMatrix to the interpContrastBase.
    % This implementation allows different fitting functions to be used for
    % the different tpup parameter types.
    % Plot the results in a display figure
    figure
    fitTypes={'nearestinterp','nearestinterp','exp1',...
        'nearestinterp','exp1','nearestinterp','nearestinterp'};

    for pp=1:nParams
        fitObject=fit(log10(contrastbase)',observedParamMatrix(pp,:)',fitTypes{pp});
        interpParams=fitObject(log10(interpContrastBase));
        interpParams(find(interpParams < vlb(pp))) = vlb(pp);
        interpParams(find(interpParams > vub(pp))) = vub(pp);
        paramLookUpMatrix(pp,:)=interpParams;
        subplot(ceil(nParams/2),2,pp);
        plot(log10(interpContrastBase),interpParams); hold on
        plot(log10(contrastbase),observedParamMatrix(pp,:),'bo'); hold off
    end
    
    % Assemble the fcon structure for this session / subject
    fcon.contrastbase=log10(interpContrastBase);
    fcon.paramLookUpMatrix=paramLookUpMatrix;
    fcon.modelObjHandle=tpupModel;
    fcon.defaultParams=tpupModel.defaultParams('defaultParamsInfo',defaultParamsInfo);
    
    % Loop over runs
    nRuns=size(mergedPacketCellArray{1,ss});
    for rr=1:nRuns
        
        fprintf('* Subject, run <strong>%g</strong> , <strong>%g</strong>\n', ss, rr);

        % Get the packet for this run
        theRunPacket=mergedPacketCellArray{1,ss}{rr};
        
        % Loop over individual instances
        nInstances=size(theRunPacket.stimulus.values,1);
        for ii=1:nInstances
            
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
                fconModel.fitResponse(theInstancePacket, ...
                'defaultParamsInfo', defaultParamsInfo, ...
                'DiffMinChange', DiffMinChange);
            
            % store the effective contrast from the paramsFit in an appropriate
            % variable
            
            myResultsVariable(ss,rr,ii)=paramsFit.paramMainMatrix; % this should be a scalar
            
            % if the param fit has hit the boundary of the avaialable
            % parameter space, report this
            if (myResultsVariable(ss,rr,ii) < min(log10(interpContrastBase)) || ...
                    myResultsVariable(ss,rr,ii) > max(log10(interpContrastBase)) )
                 warningText=['Hit effective contrast boundary'];
                 warning(warningText);
                 ss,rr,ii
%                fconModel.plot(modelResponseStruct)
%                fconModel.plot(theInstancePacket.response,'NewWindow',false)
            end
            
            % if this is the first run/instance, add a plot of the
            % available pupil responses across contrast
            if (rr==1 && ii==1)
                subplot(ceil(nParams/2),2,nParams+1); hold on;
                for cc=1:round(length(interpContrastBase)/10):length(interpContrastBase)
                    grayShade=0.1+([cc/length(interpContrastBase) cc/length(interpContrastBase) cc/length(interpContrastBase)].*0.9);
                    subclassParams=fcon.defaultParams;
                    subclassParams.paramMainMatrix = paramLookUpMatrix(:,cc)';
                    tpupResponseStruct = fcon.modelObjHandle.computeResponse(subclassParams,theInstancePacket.stimulus,[],'AddNoise',false);
                    plot(tpupResponseStruct.timebase,tpupResponseStruct.values,'Color',grayShade,'LineWidth',2);
                end
                hold off;
            end
        end % loop over instances
    end % loop over runs
end % loop over sessions

end % function
