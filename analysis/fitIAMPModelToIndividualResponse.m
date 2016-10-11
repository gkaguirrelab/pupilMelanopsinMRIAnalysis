function [theResult] = fitIAMPModelToIndividualResponse(mergedPackets)

% Loop across stimulus types (contrast level)

% Build a single packet container
% the appropriate stim type
% the kernel (which is the average response across trials for that stim
% type, available in average packets)

% Loop over stimulus instances
% place the response into the packet
% perform the model fit
% record the amplitude
% record the intial value of the response

% Loop structure
%  sessions
%    stimTypes (contrast levels plus attention task)
%       [make a packet that has a given stimulus profile and kernel]
%        instances / events


normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;

subjectKey = {'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1';	'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1'};
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};




% wanted to have a modified accumStimTypesResp that had data in it's native
% units (diameter, mm) and that began each trial at 0. Make an
% accumStimTypesResp for both filtered and unfiltered data.

filterStatus = {'unfiltered' 'filtered'};
accumStimTypes{1} = [];
accumStimTypes{2} = [];
baselineSizeCombined{1} = [];
baselineSizeCombined{2} = [];
baselineSizeCombined{3} = [];

for ff = 1:length(filterStatus);
    accumStimTypesRespRaw = [];
    for ss = 1:length(mergedPackets);
        % update process
        fprintf('>> Processing session <strong>%g</strong>\n', ss);
        % create average run from all contrasts levels to be used as a model
        
        
        % create matrices to store results
        accumStimTypesRespRaw{ss, 1} = [];
        accumStimTypesRespRaw{ss, 2} = [];
        accumStimTypesRespRaw{ss, 3} = [];
        accumStimTypesRespRaw{ss, 4} = [];
        accumStimTypesRespRaw{ss, 5} = [];
        accumStimTypesRespRaw{ss, 6} = [];
        
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
        
        NRunsTotal = length(mergedPackets{ss});
        
        % make a model for each contrast level. avgPackets doesn't quite get
        % what we want for two reasons: it's shape will be based on percent
        % amplitude change, and a given run's starting point is not 0
        
        
        
        % loop over runs
        for ii = 1:NRunsTotal
            fprintf('\t* Run <strong>%g</strong> / <strong>%g</strong>\n', ii, NRunsTotal);
            % Find the stimulus onsets so that we can align the data to it. We
            % do that by finding a [0 1] edge from a difference operator.
            tmp = sum(mergedPackets{ss}{ii}.stimulus.values);
            tmp2 = diff(tmp);
            tmp2(tmp2 < 0) = 0;
            tmp2(tmp2 > 0) = 1;
            stimOnsets = strfind(tmp2, [0 1]);
            
            % Get the number of segments from the stimulus onsets
            NSegments = length(stimOnsets);
            
            t = (0:extractionDurInd)/1000; % tells how long each stimulus should be, based on information specified in MelanopsinMR_PupilAnalysis
            
            % loop over individual stimuli
            for jj = 1:length(stimOnsets)
                
                % figure out temporal boundaries of pupil response, adjusting
                %  if stimulus happens near the end of the run
                if (stimOnsets(jj)+extractionDurInd) <= length(mergedPackets{ss}{ii}.response.values);
                    idxToExtract = stimOnsets(jj):(stimOnsets(jj)+extractionDurInd);
                else
                    idxToExtract = stimOnsets(jj):length(mergedPackets{ss}{ii}.response.values);
                end
                thisPacket.packetType = 'pupil';
                thisPacket.stimulusFile = mergedPackets{ss}{ii}.metaData.stimulusFile;
                thisPacket.responseFile = mergedPackets{ss}{ii}.metaData.responseFile;
                
                 % for storing data, identify which type of contrast was
                % associated with the stimulus
                contrast = mergedPackets{ss}{ii}.stimulus.metaData.stimTypes(jj);
                
                
                % based on timing information, actually grab pupil
                % response. again, the actual data we're going to be
                % grabbing depends if we want the filtered or unfiltered
                % stuff
                if strcmp(filterStatus{ff},'unfiltered');
                    thisPacket.respValues = mergedPackets{ss}{ii}.response.values(idxToExtract);
                    accumStimTypesRespRaw{ss,contrast} = [accumStimTypesRespRaw{ss,contrast}; (thisPacket.respValues - mergedPackets{ss}{ii}.response.values(stimOnsets(jj)))];
                end
                if strcmp(filterStatus{ff},'filtered');
                    thisPacket.respValues = mergedPackets{ss}{ii}.response.values(idxToExtract) - mergedPackets{ss}{ii}.response.metaData.lowFreqComponent(idxToExtract);
                    accumStimTypesRespRaw{ss,contrast} = [accumStimTypesRespRaw{ss,contrast}; (thisPacket.respValues - thisPacket.respValues(1))];
                end
                
                % while we're at it, also grab the baseline size. again one
                % set of values for no filtering, one set of values with
                % filtering
                
                if strcmp(filterStatus{ff},'unfiltered');
                    sizeMean = nanmean(mergedPackets{ss}{ii}.response.values);
                    baselineSize{ss,contrast} = [baselineSize{ss, contrast} mergedPackets{1,ss}{1,ii}.response.values(1,stimOnsets(jj))];
                end
                if strcmp(filterStatus{ff},'filtered');
                    sizeMean = nanmean(mergedPackets{ss}{ii}.response.metaData.lowFreqComponent);
                    baselineSize{ss,contrast} = [baselineSize{ss, contrast} (mergedPackets{1,ss}{1,ii}.response.values(1,stimOnsets(jj))-mergedPackets{1,ss}{1,ii}.response.metaData.lowFreqComponent(1,stimOnsets(jj)))];
                    baselineSize_lowFreq{ss,contrast} = [baselineSize_lowFreq{ss, contrast} mergedPackets{1,ss}{1,ii}.response.metaData.lowFreqComponent(1,stimOnsets(jj))];
                end
         
                
                
               
                
                % make new dataBin for data that has following properties: data
                % from raw data in pupil diameter (mm) (not percent change),
                % and each response starts from 0
              
                
            end % end loop over stimuli per run
        end % end loop over runs
    end % end loop over subjects/sessions
    accumStimTypes{1,ff} = accumStimTypesRespRaw;
    baselineSizeCombined{1,ff} = baselineSize;
    baselineSizeCombined{1,3} = baselineSize_lowFreq;
end % end loop over filter status

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


beta{1} = [];
beta{2} = [];

baseline{1} = [];

for ff = 1:length(filterStatus); % loop over all data filtered, then not filtered
    for ss = 1:size(accumStimTypesRespRaw,1); % Looping over subjects
        betaPerAmplitude{ss,1} = [];
        betaPerAmplitude{ss,2} = [];
        betaPerAmplitude{ss,3} = [];
        betaPerAmplitude{ss,4} = [];
        betaPerAmplitude{ss,5} = [];
        betaPerAmplitude{ss,6} = [];
      
        for cc = 1:(size(accumStimTypesRespRaw,2)) % Looping over stimulus types
            % Update the user
            fprintf('* Subject <strong>%g</strong>, Contrast <strong>%g</strong> of <strong>%g</strong>', ss, cc, 6);
            fprintf('\n');
            
            % allocate room for results
            
            
            % Update the packet to have current kernel and stimulus profile
            % Note that the kernel should be scaled to have a unit excursion
            
            singlePacket.stimulus.values = zeros(1,13000);  % blip to be convolved with kernel; fixed per subject per contrast
            singlePacket.stimulus.values(1,1) = 1;
            singlePacket.stimulus.timebase = [0:12999];
            singlePacket.kernel.values = nanmean(accumStimTypesRespRaw{ss,cc});
            singlePacket.kernel.values = singlePacket.kernel.values/(abs(min(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
            if strcmp(cc, 6);
               singlePacket.kernel.values = singlePacket.kernel.values/(abs(max(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
            end
            singlePacket.kernel.timebase = [0:12999];
            singlePacket.response.values = []; % to be filled in for each individual trial
            singlePacket.response.timebase = [0:12999];
            singlePacket.metaData = mergedPackets{1,ss}{1,1}.metaData; % steals from one of the runs of the same subject
            
            
            % Loop over instances / events
            for ii = 1:size(accumStimTypes{ff}{ss,cc},1)
                
                % Update the packet to have the response values for this
                % event. Respones values correspond to either filtered or
                % unfiltered data
                if strcmp(filterStatus{ff}, 'unfiltered') == 1;
                    singlePacket.response.values = accumStimTypes{1}{ss,cc}(ii,:);
                end
                if strcmp(filterStatus{ff}, 'filtered') == 1;
                    singlePacket.response.values = accumStimTypes{2}{ss,cc}(ii,:);
                end
                % Conduct the fit
                [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix);
                
                % Store the fitResponse
                % Save the paramsFit. You'll discover that this is a structure,
                % with one field being "amplitude". Save the value associated with
                % that.
                IAMPFitToData{ss,cc}{1,ii}.paramsFit=paramsFit;
                IAMPFitToData{ss,cc}{1,ii}.fVal=fVal;
                IAMPFitToData{ss,cc}{1,ii}.modelResponseStruct=modelResponseStruct;
                betaPerAmplitude{ss,cc} = [betaPerAmplitude{ss,cc} paramsFit.paramMainMatrix];
                
                
            end % loop over events
        end % loop over stimulus types
    end % loop over sessions
    
    beta{ff} = betaPerAmplitude;
end % loop over filter status

% Clear the object
delete(temporalFit);


% Make some plots

% Comparing amplitudes obtained with and without filtering




% Plotting Amplitude vs. Baseline Size, comparing the results of filtered
% and unfiltered data

% First determine the R2 values

for s = 1:8;
    rsqCombined{1}{s,1} = [];
    rsqCombined{1}{s,2} = [];
    rsqCombined{1}{s,3} = [];
    rsqCombined{1}{s,4} = [];
    rsqCombined{1}{s,5} = [];
    rsqCombined{1}{s,6} = [];
    rsqCombined{2}{s,1} = [];
    rsqCombined{2}{s,2} = [];
    rsqCombined{2}{s,3} = [];
    rsqCombined{2}{s,4} = [];
    rsqCombined{2}{s,5} = [];
    rsqCombined{2}{s,6} = [];
    for f = 1:length(filterStatus);
        for c = 1:6
            
            x = baselineSizeCombined{1,f}{s,c};
            y = beta{1,f}{s,c};
            
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
            rsq = corrcoef(x,y);
            rsq = rsq(1,2)*rsq(1,2);
            rsqCombined{f}{s,c} = [rsqCombined{f}{s,c}, rsq];
        end
    end
end


for s = 1:8;
    colorList = ['r', 'b', 'g', 'y', 'm'];
    contrastList = {'25%', '50%', '100%', '200%', '400%'};
    legendInfo = [[], [], [], [], []];
    plotFig = figure
    hold on
    
    for f = 1:length(filterStatus);
        
        
        subplot(1, 2, f)
        for c = 1:5
            
            x = baselineSizeCombined{1,f}{s,c};
            y = beta{1,f}{s,c};
            
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
            
            title(['Subject: ', subjectKey{s}, ', Project: ', projectKey{s}, filterStatus{f}])
            
            
            
            legendInfo{c} = [contrastList{c} [' R2 = ' num2str(rsqCombined{f}{s,c}(1))]];
            
            
            
            
        end
        legend(legendInfo)
        for c = 1:5
            
            
            x = baselineSizeCombined{1,f}{s,c};
            y = beta{1,f}{s,c};
            
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
    %saveas(plotFig, ['~/Desktop/initialAnalyses/noFilter', subjectKey{s}, '-', projectKey{s}, '-sizeByBeta'], 'png')
end

plotFig = figure;
xlim([0.75 2.25]);
title('How Filtering Affects R2 For LMS/Mel Stimuli') 
set(gca,'Xtick',1:2,'XTickLabel',{'Unfiltered', 'Filtered'});

hold on
colorList = ['r', 'b', 'g', 'y', 'm'];
contrastList = {'25%', '50%', '100%', '200%', '400%'};
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};

projectKeyLine = [];
for x = 1:length(projectKey);
    if strcmp(projectKey{x},'LMS');
        projectKeyLine{x} = '-';
    else
        projectKeyLine{x} = ':';
    end
end


for s = 1:8;
    for c = 1:5;
        
        x1 = rsqCombined{1}{s,c}(1);
        x2 = rsqCombined{2}{s,c}(1);
        xs = [x1, x2];
        
        plot(xs, projectKeyLine{s}, 'Color', colorList(c))
        plot(xs, 'o', 'Color', colorList(c))
    end
end

% same plot, but with attention
plotFig = figure;
xlim([0.75 2.25]);
title('How Filtering Affects R2 in Attention Trials') 
set(gca,'Xtick',1:2,'XTickLabel',{'Unfiltered', 'Filtered'});

hold on
colorList = ['r', 'b', 'g', 'y', 'm', 'b'];
contrastList = {'25%', '50%', '100%', '200%', '400%'};
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};

projectKeyLine = [];
for x = 1:length(projectKey);
    if strcmp(projectKey{x},'LMS');
        projectKeyLine{x} = '-';
    else
        projectKeyLine{x} = ':';
    end
end


for s = 1:8;
    for c = 6;
        
        x1 = rsqCombined{1}{s,c}(1);
        x2 = rsqCombined{2}{s,c}(1);
        xs = [x1, x2];
        
        plot(xs, projectKeyLine{s}, 'Color', colorList(c))
        plot(xs, 'o', 'Color', colorList(c))
    end
end

plotFigure = figure
hold on
title('Baseline Size')
xlabel('Unfiltered')
ylabel('Filtered')

for s = 1:8;
    for c = 1:5;
        plot(baselineSizeCombined{1}{s,c}, baselineSizeCombined{2}{s,c}, 'o')
    end
end

plotFigure = figure
hold on
title('Beta')
xlabel('Unfiltered')
ylabel('Filtered')

for s = 1:8;
    for c = 1:5;
        plot(beta{1}{s,c}, beta{2}{s,c}, 'o')
    end
end




end % function