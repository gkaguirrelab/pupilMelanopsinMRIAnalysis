function [myResultsVariable] = fitFCONModelToIndividualResponses(mergedPacketCellArray, twoComponentFitToData, dropboxAnalysisDir)
% function [myResultsVariable] = fitFCONModelToIndividualResponses(mergedPacketCellArray, twoComponentFitToData)
%
% This routine implements an effective contrast (FCON) model fit to the
% individual trials of the pupil response data.
%
% The routine takes as input the mergedPacketCellArray, which is a
% session x runs array of pupil response packets. Each run in the packet
% array corresponds to an fMRI scan during which pupil recording was
% performed. There are many events within each run.
%
% The second input is the twoComponentFitToData. This is the result of
% passing the average responses for each session / stimulus contrast to the
% TPUP model. The TPUP model returns seven parameters that define a fit to
% the mean evoked response, with a separate set of parameters for each of
% the stimulus contrast levels [25% --> 400%, log spaced].
%
% This routine loops through the sessions in the mergedPacketCellArray. For
% each session, an interpolated version of the TPUP parameters is created
% across contrast, extending the contrast range above and below that of the
% actual stimuli. The routine then loops over instances within each run and
% splits off a packet that contains just that instance. The instance packet
% is then fit with the FCON model, which searches over effective contrast
% levels, uses the interpolated parameters to obtain a set of TPUP model
% parameters corresponding to that contrast level, and then fits the data
% with the forward TPUP model given those parameters.
%
% The results are saved in myResultsVariable, which is a sessions x runs x
% instances matrix of calculated effective contrast levels. These measured
% levels are returned, and will be subjected to further analysis in other
% functions.
%
% A key hand-tuned property of this routine is to pick fitting functions
% that are used to interpolate and extend each of the TPUP parameter types
% across contrast levels.
%
% A property of the FCON model implementation is that the model expects an
% fcon structure to be included under thePacket.stimulus.metaData. This
% structure contains the information on how to expand effective contrast to
% a new set of parameters, and evaluate them with a different forward
% model.
%


%% SETUP

% Assign a name for this sub-analysis directory for saving plots and fits
subAnalysisDirectory='fitFCONModelToIndividualResponses';
outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Define the split params - These parameters define how single instances of
% response data are to be split off of the larger packet that contains the
% response data from a scanning run.
splitParams.instanceIndex=[]; % will hold the instance index, added later
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=3; % zero center the initial period, pct. change units
splitParams.normalizationWindowMsecs=100; % define the size of the norm window

% Instantiate the TPUP model - This will later be sent to FCON to use as a
% forward model. We also use it now to get information about the TPUP model params,
% which is used for setup and configuration.
tpupModel=tfeTPUP('verbosity','none');
defaultParamsInfo.nInstances=1;
[tpupModelDefaultParams, tpupVLB, tpupVUB] = ...
    tpupModel.defaultParams('defaultParamsInfo',defaultParamsInfo);
tpupVLB=tpupVLB.paramMainMatrix; % pull the upper and lower bound values
tpupVUB=tpupVUB.paramMainMatrix; %  out of the structure and into arrays

% Set up the elements of the fcon stucture.
% This contains the lookup table that maps effective contrast to the
% expanded set of parameters for the TPUP model.
contrastBase=[0.25,0.50,1.0,2.0,4.0]; % the contrast levels used in the experiment
nContrasts=length(contrastBase);
nParams=length(tpupModelDefaultParams.paramNameCell); % number of params in the TPUP model
observedParamMatrix=zeros(nParams,nContrasts); % pre-allocate this for speed later

% This is the expanded contrast that will be used for the interpolated
% parameter matrix. The range of contrasts is extended by two binary
% divisions above and below the original contrastBase.
interpContrastBase = ...
    logspace(log10(min(contrastBase)/2/2),log10(max(contrastBase)*2*2),(nContrasts+4)*100);

% fmincon search is performed by gradient descent. It will become stuck if
% small changes in the parameter (effective contrast) do not produce any
% change in the objective function. Because we are searching not over a
% continuous function but instead through a discontinuous look-up table of
% mappings between effective contrast and TPUP parameters, we need to know
% what is the size of the spacing between adjacent effective contrast
% values, and tell fmincon to have the minimum gradient step size in the
% search be at least this big. We calculate this now and pass it later to
% fmincon.
DiffMinChange=mean(diff(log10(interpContrastBase)));

%% LOOP OVER SESSIONS

% Announce what we are about to do
fprintf('>> Fitting effective contrast (FCON) / TPUP model to individual trials\n');

% instantiate the FCON model
fconModel=tfeFCON('verbosity','none');

% Loop over sessions
nSessions=size(mergedPacketCellArray,2);
for ss=1:nSessions
    
    %% Build the interpolated paramLookUpMatrix for this session /subject
    
    % First, we loop over contrast levels and grab the TPUP parameters that
    % are stored in the twoComponentFitToData, placing them in a matrix
    for cc=1:nContrasts
        observedParamMatrix(:,cc)=twoComponentFitToData{ss,cc}.paramsFit.paramMainMatrix;
    end
    
    % Interpolate the observedParamMatrix to the interpContrastBase.
    % This implementation allows different fitting functions to be used for
    % the different tpup parameter types. As we go, we plot the results in
    % a display figure so that we can evaluate the quality of the
    % interpolation.
    figParamHandle=figure('Name',['Model params, Session_' strtrim(num2str(ss))],'NumberTitle','off');
    fitTypes={'nearestinterp','nearestinterp','exp1',...
        'nearestinterp','exp1','nearestinterp','nearestinterp'};
    
    for pp=1:nParams        
        % For each param, create a fitObject which is the fit of a function
        % given in fitTypes) to the observed values of the TPUP parameter
        % across the studied contrast levels.
        fitObject=fit(log10(contrastBase)',observedParamMatrix(pp,:)',fitTypes{pp});
        
        % Use this fitObject to obtain the interpolated and extended
        % parameter values, and make sure they do not exceed the TPUP model
        % upper and lower bounds
        interpParams=fitObject(log10(interpContrastBase));
        interpParams(find(interpParams < tpupVLB(pp))) = tpupVLB(pp);
        interpParams(find(interpParams > tpupVUB(pp))) = tpupVUB(pp);
        paramLookUpMatrix(pp,:)=interpParams;
        
        % Make a plot of the params and interpolated params
        subplot(ceil(nParams/2),2,pp);
        plot(log10(interpContrastBase),interpParams); hold on
        title(tpupModelDefaultParams.paramNameCell{pp});
        plot(log10(contrastBase),observedParamMatrix(pp,:),'bo'); hold off
    end
    
    % Assemble the fcon structure for this session / subject.
    fcon.contrastbase=log10(interpContrastBase);
    fcon.paramLookUpMatrix=paramLookUpMatrix;
    fcon.modelObjHandle=tpupModel;
    fcon.defaultParams=tpupModel.defaultParams('defaultParamsInfo',defaultParamsInfo);
    
    % Loop over runs
    nRuns=size(mergedPacketCellArray{1,ss},2);
    figRunFitHandle=figure('Name',['Run Fits, Session_' strtrim(num2str(ss))],'NumberTitle','off');

    for rr=1:nRuns
        
        % Update the user on our progress
        fprintf('* Session/Subject, run <strong>%g</strong> , <strong>%g</strong>\n', ss, rr);
        
        % Get the packet for this run
        theRunPacket=mergedPacketCellArray{1,ss}{rr};
        
        % Adjust the response in the packet to remove the low-frequency
        % variation in pupil response
        lowFreqComponent=theRunPacket.response.metaData.lowFreqComponent;
        lowFreqComponent=lowFreqComponent-mean(lowFreqComponent);
        theRunPacket.response.values = ...
            theRunPacket.response.values - lowFreqComponent;
        
        % Pre-allocate a variable to hold the modeled response        
        modelRunResponse=theRunPacket.stimulus.timebase*0;
        
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
            
            % Add instance fit to the model response for the run.
            modelRunResponse(theInstancePacket.stimulus.metaData.idxToExtract) = ...
                modelRunResponse(theInstancePacket.stimulus.metaData.idxToExtract)+ ...
                modelResponseStruct.values;
            
            % store the effective contrast from the paramsFit
            myResultsVariable(ss,rr,ii)=paramsFit.paramMainMatrix;
            
            % if the param fit has hit the boundary of the avaialable
            % parameter space, report this. Shouldn't happen.
            if (myResultsVariable(ss,rr,ii) < min(log10(interpContrastBase)) || ...
                    myResultsVariable(ss,rr,ii) > max(log10(interpContrastBase)) )
                warningText=['Hit effective contrast boundary (ss-rr-ii): ' strtrim(num2str(ss)) '-' strtrim(num2str(rr)) '-' strtrim(num2str(ii)) ];
                warning(warningText);
            end
                        
            % If this is the first run/instance, add a plot of the
            % available pupil responses across contrast
            if (rr==1 && ii==1)
                figure(figParamHandle);
                subplot(ceil(nParams/2),2,nParams+1); hold on
                title('synthetic pupil responses');
                for cc=1:round(length(interpContrastBase)/10):length(interpContrastBase)
                    grayShade=0.1+([cc/length(interpContrastBase) cc/length(interpContrastBase) cc/length(interpContrastBase)].*0.9);
                    subclassParams=fcon.defaultParams;
                    subclassParams.paramMainMatrix = paramLookUpMatrix(:,cc)';
                    tpupResponseStruct = fcon.modelObjHandle.computeResponse(subclassParams,theInstancePacket.stimulus,[],'AddNoise',false);
                    plot(tpupResponseStruct.timebase,tpupResponseStruct.values,'Color',grayShade,'LineWidth',2);
                end
                hold off
            end
            
        end % loop over instances
        
        % Plot the data and fit for this run
        figure(figRunFitHandle);
        subplot(nRuns,1,rr);
        plot(theRunPacket.response.timebase,theRunPacket.response.values,'-b');
        hold on;
        meanPupilRunResponse=nanmean(theRunPacket.response.values);
        plot(theRunPacket.stimulus.timebase,modelRunResponse+meanPupilRunResponse,'-r');
        hold off;
        
    end % loop over runs
    
    % Save the parameter and simulated pupil response figure
    plotFileName=fullfile(dropboxAnalysisDir, subAnalysisDirectory, ['InterpParamsAndResponses_Sess' strtrim(num2str(ss)) '.pdf']);
    saveas(figParamHandle, plotFileName, 'pdf');
    close(figParamHandle);

    % Save the run time-series fits
    plotFileName=fullfile(dropboxAnalysisDir, subAnalysisDirectory, ['RunFits_Sess' strtrim(num2str(ss)) '.pdf']);
    saveas(figParamHandle, plotFileName, 'pdf');
    close(figParamHandle);

end % loop over sessions

% Save the fit values
dataFileName=fullfile(dropboxAnalysisDir, subAnalysisDirectory, ['effectiveContrast_SessxRunxEvent' '.mat']);
save(dataFileName,'myResultsVariable');

end % function
