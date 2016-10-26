function [beta, baselineSizeCombined, rsqCombined, slopeCombined, oneBefore] = fitIAMPModelToIndividualResponse(mergedPacketCellArray, dropboxAnalysisDir)


% Note about result structure: result{normFlag,filter}{subject,contrast}
% oneBefore is a little bit different: oneBefore is
% {result}{normFlag,filter}{subject}. each subject has as {stimuli,stimuli
% presented before} matrix, with each cell being result values of
% corresponding trials

subAnalysisDirectory='fitIAMPModelToAverageResponse';


subjectKey = [];
projectKey = [];
for ss = 1:size(mergedPacketCellArray,2);
    subjectKey{ss} = mergedPacketCellArray{ss}{1}.metaData.subjectName;
    projectKey{ss} = mergedPacketCellArray{ss}{1}.metaData.projectName;
end



filterStatus = {'unfiltered' 'filtered'};
normFlagStatus = {2,3};

% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normalizationWindowMsecs=100; % define the size of the norm window

% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];
IAMPFitToData = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;


% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

% Announce what we are about to do
fprintf('>> Fitting individual amplitude model to data (IAMP)\n');


% set up some variables
for x = 1:2 % for making oneBefore variable with 1 = beta, 2 = baseline size
    for nn = 1:length(normFlagStatus);
        for ff = 1:length(filterStatus);
            beta{nn,ff} = [];
            baselineSizeCombined{nn,ff} = [];
            for ss = 1:size(mergedPacketCellArray,2)
                for presentedStimuli = 1:6;
                    for stimuliOneBefore = 1:6;
                        oneBefore{x}{nn,ff}{ss}{presentedStimuli,stimuliOneBefore} = [];
                        oneBefore{x}{nn,3}{ss}{presentedStimuli,stimuliOneBefore} = [];
                    end
                end
            end
            
        end
    end
end


check = [];

for nn = 1:length(normFlagStatus);
    if normFlagStatus{nn} == 2;
        splitParams.normFlag = 2;
        normFlag = 2;
    end
    if normFlagStatus{nn} == 3;
        splitParams.normFlag = 3;
        normFlag = 3;
    end
    for ff = 1:length(filterStatus); % loop over all data filtered, then not filtered
        if strcmp(filterStatus{ff},'unfiltered');
            [ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, 'LowFreqClean', false, 'aggregateMethod', 'median');
        end
        if strcmp(filterStatus{ff},'filtered');
            [ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, 'LowFreqClean', true, 'aggregateMethod', 'median');
        end
        for ss = 1:size(mergedPacketCellArray,2); % Looping over subjects
            betaPerAmplitude{ss,1} = [];
            betaPerAmplitude{ss,2} = [];
            betaPerAmplitude{ss,3} = [];
            betaPerAmplitude{ss,4} = [];
            betaPerAmplitude{ss,5} = [];
            betaPerAmplitude{ss,6} = [];
            
            baselineSize{ss, 1} = [];
            baselineSize{ss, 2} = [];
            baselineSize{ss, 3} = [];
            baselineSize{ss, 4} = [];
            baselineSize{ss, 5} = [];
            baselineSize{ss, 6} = [];
            baselineSize_lowFreq{ss, 1} = [];
            baselineSize_lowFreq{ss, 2} = [];
            baselineSize_lowFreq{ss, 3} = [];
            baselineSize_lowFreq{ss, 4} = [];
            baselineSize_lowFreq{ss, 5} = [];
            baselineSize_lowFreq{ss, 6} = [];
            
            for rr = 1:(size(mergedPacketCellArray{ss},2)) % Looping over runs
                % Update the user
                fprintf('* Subject <strong>%g</strong>, Run <strong>%g</strong> of <strong>%g</strong>', ss, rr, size(mergedPacketCellArray{ss},2));
                fprintf('\n');
                
                % grab a packet that corresponds to a run for a given subject
                theRunPacket=mergedPacketCellArray{1,ss}{rr};
                if strcmp(filterStatus{ff},'filtered');
                    lowFreqComponent=theRunPacket.response.metaData.lowFreqComponent;
                    lowFreqComponent=lowFreqComponent-nanmean(lowFreqComponent);
                    theRunPacket.response.values = ...
                        theRunPacket.response.values - lowFreqComponent;
                end
                
                contrastList = [];
                % Loop over instances / events
                for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1)
                    
                    % Update user about instance
                    %                fprintf('***Instance <strong>%g</strong> of <strong>%g</strong>', ii, size(mergedPacketCellArray{ss}{rr}.stimulus.values,1));
                    %                fprintf('\n');
                    
                    % update the splitParams with the instance index
                    splitParams.instanceIndex = ii;
                    
                    % grab the packet for this particular instance
                    singlePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
                    
                    % identify contrast associated with that instance
                    contrast = theRunPacket.stimulus.metaData.stimTypes(ii);
                    
                    % store contrast value in contrastList for use with
                    % oneBefore variable
                    contrastList = [contrastList, contrast];
                    
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
                    
                    % determine baseline size
                   
                        sizeMean = nanmean(theRunPacket.response.metaData.lowFreqComponent);
                        baselineSize{ss,contrast} = [baselineSize{ss, contrast} (theRunPacket.response.values(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs))];
                        baselineSize_lowFreq{ss,contrast} = [baselineSize_lowFreq{ss, contrast} theRunPacket.response.metaData.lowFreqComponent(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs)];
                    
                        % store for oneBefore
                        
                        if ii == 1;
                        else
                            
                        oneBefore{2}{nn,ff}{ss}{contrast,contrastList(ii-1)} = [oneBefore{2}{nn,ff}{ss}{contrast,contrastList(ii-1)}, (theRunPacket.response.values(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs))];
                        oneBefore{2}{nn,3}{ss}{contrast,contrastList(ii-1)} = [oneBefore{2}{nn,3}{ss}{contrast,contrastList(ii-1)}, theRunPacket.response.metaData.lowFreqComponent(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs)];
                        end
                        
                        % Conduct the fit
                    [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix);
                    
                    % Store the fitResponse
                    % Save the paramsFit. You'll discover that this is a structure,
                    % with one field being "amplitude". Save the value associated with
                    % that.
                    IAMPFitToData{ss,contrast}{1,ii}.paramsFit=paramsFit;
                    IAMPFitToData{ss,contrast}{1,ii}.fVal=fVal;
                    IAMPFitToData{ss,contrast}{1,ii}.modelResponseStruct=modelResponseStruct;
                    betaPerAmplitude{ss,contrast} = [betaPerAmplitude{ss,contrast} paramsFit.paramMainMatrix];
                    
                    if ii == 1;
                        else
                            
                        oneBefore{1}{nn,ff}{ss}{contrast,contrastList(ii-1)} = [oneBefore{2}{nn,ff}{ss}{contrast,contrastList(ii-1)}, paramsFit.paramMainMatrix];
                    end
                    
                    check = [check theRunPacket.response.values(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs)-singlePacket.metaData.splitOffAnInstance.normValue];
                    
                end % loop over events
            end % loop over stimulus types
        end % loop over sessions
        
        beta{nn,ff} = betaPerAmplitude;
        
        baselineSizeCombined{nn,ff} = baselineSize;
        baselineSizeCombined{nn,3} = baselineSize_lowFreq;
        
    end % loop over filter status
end % loop over normFlag status

% Clear the object
delete(temporalFit);

% First determine the R2 values


subjectsLMS = [];
subjectsMel = [];
for x = 1:length(projectKey);
    if strcmp(projectKey(x), 'MelanopsinMRMaxLMSCRF');
        subjectsLMS = [subjectsLMS x];
    end
    if strcmp(projectKey(x), 'MelanopsinMRMaxMelCRF');
        subjectsMel = [subjectsMel x];
    end
end

% First determine R2 and slope
for nn = 1:length(normFlagStatus);
    for ff = 1:length(filterStatus);
       for s = 1:size(mergedPacketCellArray,2);
            for c = 1:6;
                rsqCombined{nn,ff}{s,c} = [];
                slopeCombined{nn,ff}{s,c} = [];
            end
        end
    end
end


for s = 1:size(mergedPacketCellArray,2);
    for n = 1:length(normFlagStatus);
        for f = 1:length(filterStatus);
            for c = 1:6
                
                x = baselineSizeCombined{n,f}{s,c};
                y = beta{n,f}{s,c};
                
                xnan = isnan(x);
                
                xnanlist = [];
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                x = [ones(1,length(x)); x];
                x = x';
                y = y';
                [B,BINT,R,RINT,STATS] = regress(y,x);
                slope = B(1)
                rsq = STATS(1);
                slopeCombined{n,f}{s,c} = [slopeCombined{n,f}{s,c}, slope];
                rsqCombined{n,f}{s,c} = [rsqCombined{n,f}{s,c}, rsq];
            end
        end
    end
end

%% Actually do the plotting. We want 4 sets of plots: filtered and not filtered, raw values and percent change

% to save simulated data according to simulation method

if isfield(mergedPacketCellArray{1}{1}.metaData, 'simulationStyle');
    if strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'arousal');
        simFileName = '_simulatedArousal';
    elseif strcmp(mergedPacketCellArray{1}{1}.metaData.simulationStyle, 'carryOver')
        simFileName = '_simulatedCarryOver';
    
    end
else
        simFileName = '';
end



% plot LMS subjects


for nn = 1:length(normFlagStatus);
    for f = 1:length(filterStatus);
        plotFig = figure;
        
        
        %title(['LMS Subjects with ' filterStatus{f} ' data, with normFlag ' num2str(normFlagStatus{nn})])
        
        for s = subjectsLMS;
            colorList = ['r', 'b', 'g', 'y', 'm'];
            contrastList = {'25%', '50%', '100%', '200%', '400%'};
            legendInfo = [[], [], [], [], []];
            
            hold on
            
            
            
            if length(subjectsLMS)>1
                subplot(2, length(subjectsLMS)/2, s)
            end
            for c = 1:5
                
                x = baselineSizeCombined{nn,f}{s,c};
                y = beta{nn,f}{s,c};
                
                xnan = isnan(x);
                
                xnanlist = [];
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                coeffs = polyfit(x, y, 1);
                fittedX = linspace(min(x), max(x), 200);
                fittedY = polyval(coeffs, fittedX);
                plot(x, y,'o', 'Color', colorList(c))
                hold on
                
                %pbaspect([1 1 1]);
                xlabel('Baseline Size');
                ylabel('Beta');
                
                title(['Subject: ', subjectKey{s}])
                
                
                
                legendInfo{c} = [contrastList{c} [' R2 = ' num2str(rsqCombined{nn,f}{s,c}(1))]];
                
                
                
                
            end
            %legend(legendInfo)
            for c = 1:5
                
                
                x = baselineSizeCombined{nn,f}{s,c};
                y = beta{nn,f}{s,c};
                
                xnan = isnan(x);
                
                xnanlist = [];
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                coeffs = polyfit(x, y, 1);
                fittedX = linspace(min(x), max(x), 200);
                fittedY = polyval(coeffs, fittedX);
                plot(fittedX, fittedY, 'LineWidth', 3, 'Color', colorList(c))
            end
        end
        outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        saveas(plotFig, fullfile(outDir, ['baselineSizeByAmplitude_LMS_normFlag_', num2str(normFlagStatus{nn}), '_', filterStatus{f}, simFileName, '.png']), 'png');
        close(plotFig);
        
        % now plotting the Melanopsin subjects
        % we're plotting only the data that has the low frequency component removed
        plotFig = figure;
        %title(['Mel Subjects with ' filterStatus{f} ' data, with normFlag ' num2str(normFlagStatus{nn})])

        for s = subjectsMel;
            colorList = ['r', 'b', 'g', 'y', 'm'];
            contrastList = {'25%', '50%', '100%', '200%', '400%'};
            legendInfo = [[], [], [], [], []];
            
            hold on
            
            
            
            if length(subjectsMel)>1
                subplot(2, length(subjectsMel)/2, find(subjectsMel==s))
            end
            for c = 1:5
                
                x = baselineSizeCombined{nn,f}{s,c};
                y = beta{nn,f}{s,c};
                
                xnan = isnan(x);
                
                xnanlist = [];
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                coeffs = polyfit(x, y, 1);
                fittedX = linspace(min(x), max(x), 200);
                fittedY = polyval(coeffs, fittedX);
                plot(x, y,'o', 'Color', colorList(c))
                hold on
                
                %pbaspect([1 1 1]);
                xlabel('Baseline Size');
                ylabel('Beta');
                
                title(['Subject: ', subjectKey{s}])
                
                
                
                legendInfo{c} = [contrastList{c} [' R2 = ' num2str(rsqCombined{nn,f}{s,c}(1))]];
                
                
                
                
            end
            %legend(legendInfo)
            for c = 1:5
                
                
                x = baselineSizeCombined{nn,f}{s,c};
                y = beta{nn,f}{s,c};
                
                xnan = isnan(x);
                
                xnanlist = [];
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                coeffs = polyfit(x, y, 1);
                fittedX = linspace(min(x), max(x), 200);
                fittedY = polyval(coeffs, fittedX);
                plot(fittedX, fittedY, 'LineWidth', 3, 'Color', colorList(c))
            end
        end
        
        outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        saveas(plotFig, fullfile(outDir, ['baselineSizeByAmplitude_Mel_normFlag_', num2str(normFlagStatus{nn}), '_', filterStatus{f}, simFileName, '.png']), 'png');
        close(plotFig);
    end
end


fileID = fopen([outDir '/summaryText.txt'],'w');
fprintf(fileID, 'The pupil response to a given light stimulus is related to the contrast of the light stimulus; more contrast \nevokes a larger constriction.\n\n');
fprintf(fileID, 'Baseline pupil size can also predict amplitude of the pupil response to a light stimulus. More precisely, \nsmaller pupils contract less to a light stimulus. If the baseline pupil size and amplitude of pupil \nresponse are related by some slope, the effect of increasing contrast appears to be increasing the \nmagnitude of the slope.\n\n');
fprintf(fileID, 'The relationship between baseline pupil size and amplitude of the pupil response to light is much stronger \non pupil data that has the low frequency component removed. One possible model for this observation \nis that two independent processes contribute to determining pupil size, a slow wave component and a \ncomponent that determines pupil size as a function of irradiance. When the brain decides how much to \nrespond to a brief flash of light, it scales the degree of response by how big the pupil should as determined \nby irradiance, not actual pupil size.\n');
fclose(fileID);

end % function