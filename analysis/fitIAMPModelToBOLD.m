function [percentSignalChange, percentSignalChange_byStimulus] = fitIAMPModelToBOLD(mergedPacketCellArray, dropboxAnalysisDir)



close all;
projectList = {'lms' 'mel'};
subAnalysisDirectory = 'BOLD/expanded';
% load in the BOLD data
melBOLD = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/fmriMelanopsinMRIAnalysis/packetCache/MelanopsinMR_MelCRF_15a230960983db3a019db2778fe9f27d.mat');
melBOLD = melBOLD.packetCellArray;
lmsBOLD = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/fmriMelanopsinMRIAnalysis/packetCache/MelanopsinMR_LMSCRF_5bf6b9f0b02650f5734d3f0060e3d068.mat');
lmsBOLD = lmsBOLD.packetCellArray;

% pre-allocate space for results saving
for ss = 1:3;
    for project = 1:length(projectList);
        if project == 1; % corresponding to LMS
            adjustmentFactor = 0;
            
        elseif project == 2; % corresponding to mel
            adjustmentFactor = 4;
            
        end
        for contrast = 1:6;
            percentSignalChange_byStimulus{ss+adjustmentFactor,contrast} = [];
        end
    end
end

for ss = 1:3;
    for project = 1:length(projectList);
        if project == 1; % corresponding to LMS
            adjustmentFactor = 0;
            boldPacket = lmsBOLD;
        elseif project == 2; % corresponding to mel
            adjustmentFactor = 4;
            boldPacket = melBOLD;
        end
        
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            
            thePacket=boldPacket{ss,rr};
            
            % Construct the model object
            temporalFit = tfeIAMP('verbosity','none');
            
            % grab the average hrf and prepare it as a kernel
            check = diff(thePacket.response.timebase);
            responseDeltaT = check(1);
            nSamples = ceil((thePacket.kernel.timebase(end)-thePacket.kernel.timebase(1))/responseDeltaT);
            newKernelTimebase = thePacket.kernel.timebase(1):responseDeltaT:(thePacket.kernel.timebase(1)+nSamples*responseDeltaT);
            hrfKernelStruct = temporalFit.resampleTimebase(thePacket.kernel,newKernelTimebase);
            thePacket.kernel=prepareHRFKernel(hrfKernelStruct);
            
            % Put the TR into the packet metaData
            thePacket.response.metaData.TRmsecs=800;
            
            % downsample the stimulus values to 100 ms deltaT to speed things up
            totalResponseDuration=thePacket.response.metaData.TRmsecs * ...
                length(thePacket.response.values);
            newStimulusTimebase=linspace(0,totalResponseDuration-100,totalResponseDuration/100);
            thePacket.stimulus=temporalFit.resampleTimebase(thePacket.stimulus,newStimulusTimebase);
            
            % Convert the response values to mean-centered % signal change
            response=thePacket.response.values;
            response=(response-mean(response))/mean(response);
            thePacket.response.values=response;
            
            % How many trial instances are in this packet?
            defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);
            
            % Define an empty locking matrix
            paramLockMatrix=[];
            
            % Perform the fit
            [paramsFit,FVal,modelResponseStruct] = ...
                temporalFit.fitResponse(thePacket,...
                'defaultParamsInfo', defaultParamsInfo, ...
                'paramLockMatrix',paramLockMatrix, ...
                'searchMethod','linearRegression', ...
                'errorType','1-r2');
            
            rSquared=1-FVal;
            
            
            % plot the data and the fit
            
            % store the results in two variables. one variable organizes
            % results according to (ss,rr,ii). the other variable organizes
            % result according to stimulus identity (ss,contrast)
            for ii = 1:length(paramsFit.paramMainMatrix) % loop over instances
                percentSignalChange(ss+adjustmentFactor,rr,ii) = paramsFit.paramMainMatrix(ii);
                contrastLevel = mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.metaData.stimTypes(ii);
                percentSignalChange_byStimulus{ss+adjustmentFactor,contrastLevel} = [percentSignalChange_byStimulus{ss+adjustmentFactor,contrastLevel}, paramsFit.paramMainMatrix(ii)];
            end
            %figure
            %plot(modelResponseStruct.timebase,modelResponseStruct.values)
            %hold on
            %plot(thePacket.response.timebase,thePacket.response.values)
            
            % dump out the estimated percent signal change for each trial
            paramsFit.paramMainMatrix;
        end % end loop over runs
    end % end loop over project
end % end loop over subjects

%% Now to do some plotting
% First, plot mean evoked response by contrast for each subject
% make mean percent signal change for a given contrast level


for ss = 1:3;
    plotFig = figure;
    
    
    for project = 1:length(projectList);
        if project == 1; % lms
            adjustmentFactor = 0;
            color = 'b';
        elseif project == 2; % mel
            adjustmentFactor = 4;
            color = 'r';
        end
        for contrastLevel = 1:5;
            meanPercentSignalChange(ss+adjustmentFactor,contrastLevel) = mean(percentSignalChange_byStimulus{ss+adjustmentFactor,contrastLevel});
            SEMPercentSignalChange(ss+adjustmentFactor,contrastLevel) = std(percentSignalChange_byStimulus{ss+adjustmentFactor,contrastLevel})/sqrt(length(percentSignalChange_byStimulus{ss+adjustmentFactor,contrastLevel}));

        end
        subplot(1, 2, project)
        plot([25 50 100 200 400], meanPercentSignalChange(ss+adjustmentFactor,:)*100, 'Color', color);
        h = errorbar([25 50 100 200 400],meanPercentSignalChange(ss+adjustmentFactor,1:5),SEMPercentSignalChange(ss+adjustmentFactor,1:5), color); set(h,'linestyle','-');
        xlabel('Presented Contrast (%)');
        ylabel('BOLD Percent Signal Change (%)');
        title(['Subject: ',  mergedPacketCellArray{ss+adjustmentFactor}{1}.metaData.subjectName, 'Project: ', projectList(project)])
        
    end
    outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    saveas(plotFig, fullfile(outDir, ['BOLDPercentSignalChangeByContrast_', mergedPacketCellArray{ss+adjustmentFactor}{1}.metaData.subjectName, '.png']), 'png');
    close(plotFig);
end

% make combined mean

for contrastLevel = 1:5;
    meanPercentSignalChangeCollapsed(1,contrastLevel) = sum(meanPercentSignalChange(1:3,contrastLevel))/3;
    meanPercentSignalChangeCollapsed(2,contrastLevel) = sum(meanPercentSignalChange(4:7,contrastLevel))/3;

end

plotFig = figure;
plot([25 50 100 200 400], meanPercentSignalChangeCollapsed(1,:)*100);
xlabel('Presented Contrast (%)');
ylabel('BOLD Percent Signal Change (%)');
saveas(plotFig, fullfile(outDir, ['BOLDPercentSignalChangeByContrast_LMS_mean.png']), 'png');
close(plotFig);

plotFig = figure;
plot([25 50 100 200 400], meanPercentSignalChangeCollapsed(2,:)*100);
xlabel('Presented Contrast (%)');
ylabel('BOLD Percent Signal Change (%)');
saveas(plotFig, fullfile(outDir, ['BOLDPercentSignalChangeByContrast_mel_mean.png']), 'png');
close(plotFig);


% Now plot by FCON value by BOLD percent signal change
% load up FCON results
effectiveContrast_sessxrunxevent = load('~/Dropbox-Aguirre-Brainard-Lab/MELA_analysis/pupilMelanopsinMRIAnalysis/fitFCONModelToIndividualResponses/expanded/effectiveContrast_sessxrunxevent.mat');
FCON = effectiveContrast_sessxrunxevent.myResultsVariable;

colorList = ['r', 'b', 'g', 'y', 'm'];
contrastList = {'25%', '50%', '100%', '200%', '400%'};

for ss = 1:3;
    
    
    for project = 1:length(projectList);
        if project == 1; % corresponding to LMS
            adjustmentFactor = 0;
            saveAs = 'lms';
            
        elseif project == 2; % corresponding to mel
            adjustmentFactor = 4;
            saveAs = 'mel';
            
        end
        for c = 1:6;
            x{c}=[];
            y{c}=[];
        end
        plotFig = figure;
        for rr = 1:(size(mergedPacketCellArray{ss+adjustmentFactor},2));
            for ii = 1:size(mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.values,1)
                contrastLevel = mergedPacketCellArray{ss+adjustmentFactor}{rr}.stimulus.metaData.stimTypes(ii);
                
                FCONValue = FCON(ss+adjustmentFactor,rr,ii);
                x{contrastLevel} = [x{contrastLevel}, FCONValue];
                
                BOLDValue = percentSignalChange(ss+adjustmentFactor, rr, ii)*100;
                y{contrastLevel} = [y{contrastLevel}, BOLDValue];
                if contrastLevel == 6;
                    % do nothing for attention trials
                else
                    subplot(2,3,contrastLevel);
                    title([contrastList(contrastLevel)])
                    hold on
                    plot(FCONValue,BOLDValue, 'o', 'Color', colorList(contrastLevel));
                    subplot(2,3,6);
                    title('Combined')
                    hold on
                    plot(FCONValue,BOLDValue, 'o', 'Color', colorList(contrastLevel))
                    
                end
            end
        end
        % determine max percent signal change to set axes consistent across
        % presented stimulus types
        boldBucket = [];
        for contrastLevel = 1:5;
            boldBucket = [boldBucket, y{contrastLevel}];
            maxBOLD = max(boldBucket);
            minBOLD = min(boldBucket);
        end
        
        for c = 1:5;
            coeffs = polyfit(x{c}, y{c}, 1);
            fittedX = linspace(min(x{c}), max(x{c}), 200);
            fittedY = polyval(coeffs, fittedX);
            subplot(2,3,c);
            plot(fittedX, fittedY, 'LineWidth', 3, 'Color', colorList(c))
            xlabel('FCON (log(%))');
            ylabel('BOLD Percent Signal Change (%)');
            hold on
            ylim([minBOLD maxBOLD]);
            
            
        end
        outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        saveas(plotFig, fullfile(outDir, ['BOLDPercentSignalChangeByFCON_', mergedPacketCellArray{ss+adjustmentFactor}{1}.metaData.subjectName, '_', saveAs, '.png']), 'png');
        close(plotFig);
    end
end

end

