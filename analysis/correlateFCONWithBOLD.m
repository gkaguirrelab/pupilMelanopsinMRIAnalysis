function [FCONCovariateValues, FCONCovariateTimeseries, rsqLag, meanRsqLag, meanRsqLagCombined, xcorrMean] = correlateFCONWithBOLD(mergedPacketCellArray, dropboxAnalysisDir)

subAnalysisDirectory = 'BOLD';

%% For a given run, create a vector that represents the the difference (in
% change units) between the FCON estimate of the response on that trial and
% the mean FCON estimate across all trials for that stimulus contrast

% load up the FCON results
effectiveContrast_sessxrunxevent = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/pupilMelanopsinMRIAnalysis/fitFCONModelToIndividualResponses/effectiveContrast_sessxrunxevent.mat');
FCON = effectiveContrast_sessxrunxevent.myResultsVariable;

% load up the MRI data
melBOLD = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/fmriMelanopsinMRIAnalysis/packetCache/MelanopsinMR_MelCRF_15a230960983db3a019db2778fe9f27d.mat');
melBOLD = melBOLD.packetCellArray;
lmsBOLD = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/fmriMelanopsinMRIAnalysis/packetCache/MelanopsinMR_LMSCRF_5bf6b9f0b02650f5734d3f0060e3d068.mat');
lmsBOLD = lmsBOLD.packetCellArray;

% pre-allocate room for data storage
for ss = 1:3; % loop over subjects. for now just looping over the first three subjects while there remains ambiguitiy over the fMRI data structure
    for project = 1:2; % loop over LMS (1) than melanopsin (2)
        if project == 1;
            adjustmentFactor = 0;
        elseif project == 2;
            adjustmentFactor = 4; % on mergedPacketCellArray, session 1 is LMS while session 5 is the same subject but melanopsin
        end
        
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            FCONCovariateValues{ss,project}{rr} = [];
            FCONCovariateTimeseries{ss,project}{rr} = zeros(1,length(mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.timebase));
            rsqLag{ss, rr} = [];
            rsqLag{ss+adjustmentFactor, rr} = [];
            xcorrCollapsed{ss+adjustmentFactor,rr} = [];
        end
    end
    
    
end




for ss = 1:3; % loop over subjects. for now just looping over the first three subjects while there remains ambiguitiy over the fMRI data structure
    for project = 1:2; % loop over LMS (1) than melanopsin (2)
        if project == 1;
            adjustmentFactor = 0;
        elseif project == 2;
            adjustmentFactor = 4; % on mergedPacketCellArray, session 1 is LMS while session 5 is the same subject but melanopsin
        end
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            
            for contrastTypes = 1:6;
                FCONPerRun{contrastTypes} = [];
            end
            % organize FCON results by presented stimulus, for purpose of
            % determining mean FCON value for a given stimulus type
            for ii = 1:size(mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.values,1);
                contrast = mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.metaData.stimTypes(ii);
                fconValue = FCON(ss+adjustmentFactor,rr,ii);
                
                FCONPerRun{contrast} = [FCONPerRun{contrast}, fconValue];
                
            end
            
            % determine mean FCON value for a given stimulus type
            for contrastTypes = 1:6;
                meanFCONPerRun(contrastTypes) = mean(FCONPerRun{contrastTypes});
            end
            % make covariate
            for ii = 1:size(mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.values,1);
                contrast = mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.metaData.stimTypes(ii);
                fconValue = FCON(ss+adjustmentFactor,rr,ii);
                
                FCONCovariateValues{ss,project}{rr} = [FCONCovariateValues{ss,project}{rr}, fconValue - meanFCONPerRun(contrast)];
                
                % convert these values to a timeseries covariate
                % first obtain time series of the stimulus profile for the
                % instance
                stimuliTimeseries = mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.values(ii,:);
                
                % adjust amplitude of event by multiplying by FCON
                % difference value
                adjustedStimuliTimeseries = FCONCovariateValues{ss,project}{rr}(ii) * stimuliTimeseries;
                
                % sum timeseries for
                if contrast == 6;
                    % do nothing
                else
                    FCONCovariateTimeseries{ss,project}{rr} = FCONCovariateTimeseries{ss,project}{rr} + adjustedStimuliTimeseries;
                end
            end
        end % end loop over runs
    end % end loop over projects
end % end loop over subjects

%% Now we have our covariate, time to move onto correlating it with BOLD data

for ss = 1:3; % loop over subjects. for now just looping over the first three subjects while there remains ambiguitiy over the fMRI data structure
    for project = 1:2; % loop over LMS (1) than melanopsin (2)
        if project == 1;
            adjustmentFactor = 0;
            boldPacket = lmsBOLD;
        elseif project == 2;
            adjustmentFactor = 4; % on mergedPacketCellArray, session 1 is LMS while session 5 is the same subject but melanopsin
            boldPacket = melBOLD;
        end
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            % perform convolution
            % first setup single run packet
            responseStruct.timebase = mergedPacketCellArray{1,ss+adjustmentFactor}{1,rr}.response.timebase;
            responseStruct.metaData = mergedPacketCellArray{1,ss+adjustmentFactor}{1,rr}.response.metaData;
            responseStruct.values = FCONCovariateTimeseries{ss,project}{rr};
            % perform convolution
            temporalFit = tfeIAMP('verbosity','none');
            convolutionStruct = temporalFit.applyKernel(responseStruct, boldPacket{ss,rr}.kernel);
            % downsample convolution to match TR of MRI data
            for x = 1:length(boldPacket{ss,rr}.response.timebase);
                convolutionStruct.timebaseResampled = boldPacket{ss,rr}.response.timebase;
                index = find(convolutionStruct.timebase==round(convolutionStruct.timebaseResampled(x)));
                
                % for convolution, use the next two lines:
                convolutionStruct.valuesResampled(x) = convolutionStruct.values(index);
                
            end
            
            % perform correlation with BOLD data
            [acorr, lag] = xcorr(convolutionStruct.valuesResampled,boldPacket{ss,rr}.response.values);
            xcorrCollapsed{ss+adjustmentFactor,rr} = acorr;
            
            % perform lag correlation analyses by hand
            lagRange = -50:50;
            for lagIndex = lagRange;
                laggedFCONcovariate = [];
                for x = 1:length(convolutionStruct.valuesResampled);
                    if x-lagIndex > 0 && (x-lagIndex) < (length(convolutionStruct.valuesResampled)+1);
                        laggedFCONcovariate(x) = convolutionStruct.valuesResampled(x-lagIndex);
                        % so a lagIndex of +1 means the FCON covariate
                        % shifted forward in time (new value at time t = 2
                        % corresponds to original t = 1)
                        if lagIndex > 0;
                            laggedFCONcovariate(1:lagIndex) = NaN;
                        end
                        if lagIndex < 0;
                            laggedFCONcovariate(length(convolutionStruct.valuesResampled)+1+lagIndex:length(convolutionStruct.valuesResampled)) = NaN;
                        end
                    end
                end
                % calculate R2 for a given lagIndex
                x = laggedFCONcovariate;
                y = boldPacket{ss,rr}.response.values;
                xnan = isnan(x);
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                
                
                rsq = corr2(x,y);
                rsqLag{(ss+adjustmentFactor), rr}(1,(find(lagRange==lagIndex))) = lagIndex;
                rsqLag{(ss+adjustmentFactor), rr}(2,(find(lagRange==lagIndex))) = rsq;
            end
            
            
        end
    end
end


%% do some plotting
% collapse across runs
for ss = 1:3;
    for project = 1:2;
        if project == 1;
            adjustmentFactor = 0;
        elseif project == 2;
            adjustmentFactor = 4;
        end
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            
            for x = lagRange;
                rsqLagCollapsed{ss+adjustmentFactor}{find(lagRange==x)} = [];
                
            end
            for x = 1:1119;
                xcorrCollapsedAcrossRuns{ss+adjustmentFactor}{x} = [];
            end
        end
    end
end
for ss = 1:3;
    for project = 1:2;
        if project == 1;
            adjustmentFactor = 0;
        elseif project == 2;
            adjustmentFactor = 4;
        end
        
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            
            for x = lagRange;
                rsqLagCollapsed{ss+adjustmentFactor}{find(lagRange==x)} = [rsqLagCollapsed{ss+adjustmentFactor}{(find(lagRange==x))}, rsqLag{ss+adjustmentFactor,rr}(2,(find(lagRange==x)))];
                
            end
            
        end
    end
end
% find mean average pearson correlation coefficient at each lag time
meanRsqLag = zeros(8,length(lagRange));
for ss = 1:3;
    for x = lagRange;
        meanRsqLag(ss,find(lagRange==x))= mean(rsqLagCollapsed{ss}{find(lagRange==x)});
        meanRsqLag(ss+4,find(lagRange==x))= mean(rsqLagCollapsed{ss+4}{find(lagRange==x)});
        
    end
end
for x = lagRange;
    meanRsqLagCombined(1,find(lagRange==x)) = mean(meanRsqLag(1:3,find(lagRange==x)));
    meanRsqLagCombined(2,find(lagRange==x)) = mean(meanRsqLag(4:7,find(lagRange==x)));
end

% confirm with xcorr method
for ss = 1:3;
    for project = 1:2;
        if project == 1;
            adjustmentFactor = 0;
        elseif project == 2;
            adjustmentFactor = 4;
        end
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            for x = 1:length(xcorrCollapsed{ss+adjustmentFactor,rr});
                xcorrCollapsedAcrossRuns{ss+adjustmentFactor}{x} = [xcorrCollapsedAcrossRuns{ss+adjustmentFactor}{x}, xcorrCollapsed{ss+adjustmentFactor,rr}(x)];
            end
        end
    end
end

for ss = 1:3;
    for project = 1:2;
        if project == 1;
            adjustmentFactor = 0;
        elseif project == 2;
            adjustmentFactor = 4;
        end
        for x = 1:1119;
            xcorrMeanSubject(ss+adjustmentFactor,x) = mean(xcorrCollapsedAcrossRuns{ss+adjustmentFactor}{x});
        end
    end
end
for ss = 1:3;
    for x = 1:1119;
        xcorrMean(1,x) = mean(xcorrMeanSubject(1:3,x));
        xcorrMean(2,x) = mean(xcorrMeanSubject(4:7,x));
    end
end

% actually doing the plotting
% first average overall
plotFig = figure;
plot(lagRange*0.8,meanRsqLagCombined(1,:), 'Color', 'r')
hold on
plot(lagRange*0.8,meanRsqLagCombined(2,:), 'Color', 'b')
xlabel('Lag (s)');
ylabel('Pearson Correlation Coeffcient');
legend('LMS', 'Mel');
outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

saveas(plotFig, fullfile(outDir, ['correlateFCONWithBOLD' '.png']), 'png');
close(plotFig);

% second displaying per subject
plotFig = figure;
hold on
for ss = 1:3;
    for project = 1:2;
        if project == 1;
            adjustmentFactor = 0;
            color = 'b';
        elseif project == 2;
            adjustmentFactor = 4;
            color = 'r';
        end
        plot(lagRange*0.8,meanRsqLag(ss+adjustmentFactor,:), 'Color', color');
    end
end
xlabel('Lag (s)');
ylabel('Pearson Correlation Coeffcient');
legend('LMS', 'Mel');
saveas(plotFig, fullfile(outDir, ['correlateFCONWithBOLD_perSubject' '.png']), 'png');
close(plotFig);

end % end function






