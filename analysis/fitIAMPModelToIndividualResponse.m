function [theResult] = fitIAMPModelToIndividualResponse(mergedPacketCellArray, avgPackets, dropboxAnalysisDir)

subAnalysisDirectory='fitIAMPModelToAverageResponse';

normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;

subjectKey = {'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1';	'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1'};
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};
filterStatus = {'unfiltered' 'filtered'};

% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=2; % zero center the initial period, change units
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


beta{1} = [];
beta{2} = [];

baselineSizeCombined{1} = [];
baselineSizeCombined{2} = [];

for ff = 1:length(filterStatus); % loop over all data filtered, then not filtered
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
            
            
            
            
            
            % Loop over instances / events
            for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1)
                
                % Update user about instance
                fprintf('***Instance <strong>%g</strong> of <strong>%g</strong>', ii, size(mergedPacketCellArray{ss}{rr}.stimulus.values,1));
                fprintf('\n');
                
                % update the splitParams with the instance index
                splitParams.instanceIndex = ii;
            
                % grab the packet for this particular instance
                singlePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
                
                % identify contrast associated with that instance
                contrast = theRunPacket.stimulus.metaData.stimTypes(ii);
              
                % create kernel, which is the average run for a contrast
                % level per subject
                singlePacket.kernel.values = avgPackets{ss,contrast}.response.values;
                singlePacket.kernel.values = singlePacket.kernel.values/(abs(min(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
                if strcmp(contrast, 6);
                    singlePacket.kernel.values = singlePacket.kernel.values/(abs(max(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
                end
                singlePacket.kernel.timebase = [0:12999];
                
                % determine baseline size
                 if strcmp(filterStatus{ff},'unfiltered');
                    sizeMean = nanmean(theRunPacket.response.values);
                    baselineSize{ss,contrast} = [baselineSize{ss, contrast} theRunPacket.response.values(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs)];
                end
                if strcmp(filterStatus{ff},'filtered');
                    sizeMean = nanmean(theRunPacket.response.metaData.lowFreqComponent);
                    baselineSize{ss,contrast} = [baselineSize{ss, contrast} (theRunPacket.response.values(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs)-theRunPacket.response.metaData.lowFreqComponent(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs))];
                    baselineSize_lowFreq{ss,contrast} = [baselineSize_lowFreq{ss, contrast} theRunPacket.response.metaData.lowFreqComponent(1,singlePacket.metaData.splitOffAnInstance.splitOnsetMsecs)];
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
                
                
            end % loop over events
        end % loop over stimulus types
    end % loop over sessions
    
    beta{ff} = betaPerAmplitude;
    
    baselineSizeCombined{1,ff} = baselineSize;
    baselineSizeCombined{1,3} = baselineSize_lowFreq;
    
end % loop over filter status

% Clear the object
delete(temporalFit);

% First determine the R2 values

subjectKey = {'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1';	'HERO_asb1';	'HERO_aso1';	'HERO_gka1';	'HERO_mxs1'};
projectKey = {'LMS'; 'LMS'; 'LMS'; 'LMS'; 'Melanopsin'; 'Melanopsin';	'Melanopsin';	'Melanopsin'};
subjectsLMS = [];
subjectsMel = [];
for x = 1:length(projectKey);
    if strcmp(projectKey(x), 'LMS');
        subjectsLMS = [subjectsLMS x];
    end
    if strcmp(projectKey(x), 'Melanopsin');
        subjectsMel = [subjectsMel x];
    end
end

% First determine R2
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


% plot LMS subjects
f = 2; % we're plotting only the data that has the low frequency component removed
plotFig = figure;
for s = subjectsLMS;
    colorList = ['r', 'b', 'g', 'y', 'm'];
    contrastList = {'25%', '50%', '100%', '200%', '400%'};
    legendInfo = [[], [], [], [], []];
    
    hold on
    
    
    
    
    subplot(1, 4, s)
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
        
        title(['Subject: ', subjectKey{s}, ', Project: ', projectKey{s}, ' ', filterStatus{f}])
        
        
        
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
outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, 'baselineSizeByAmplitude_filtered_LMS.png'), 'png');
close(plotFig);

% now plotting the Melanopsin subjects
f = 2; % we're plotting only the data that has the low frequency component removed
plotFig = figure;
for s = subjectsMel;
    colorList = ['r', 'b', 'g', 'y', 'm'];
    contrastList = {'25%', '50%', '100%', '200%', '400%'};
    legendInfo = [[], [], [], [], []];
    
    hold on
    
    
    
    
    subplot(1, 4, s-4)
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
        
        title(['Subject: ', subjectKey{s}, ', Project: ', projectKey{s}, ' ', filterStatus{f}])
        
        
        
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

outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
saveas(plotFig, fullfile(outDir, 'baselineSizeByAmplitude_filtered_Mel.png'), 'png');
close(plotFig);

fileID = fopen([outDir '/summaryText.txt'],'w');
fprintf(fileID, 'The pupil response to a given light stimulus is related to the contrast of the light stimulus; more contrast \nevokes a larger constriction.\n\n');
fprintf(fileID, 'Baseline pupil size can also predict amplitude of the pupil response to a light stimulus. More precisely, \nsmaller pupils contract less to a light stimulus. If the baseline pupil size and amplitude of pupil \nresponse are related by some slope, the effect of increasing contrast appears to be increasing the \nmagnitude of the slope.\n\n');
fprintf(fileID, 'The relationship between baseline pupil size and amplitude of the pupil response to light is much stronger \non pupil data that has the low frequency component removed. One possible model for this observation \nis that two independent processes contribute to determining pupil size, a slow wave component and a \ncomponent that determines pupil size as a function of irradiance. When the brain decides how much to \nrespond to a brief flash of light, it scales the degree of response by how big the pupil should as determined \nby irradiance, not actual pupil size.\n');
fclose(fileID);

end % function