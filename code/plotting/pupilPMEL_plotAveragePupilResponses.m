function [ ] = pupilPMEL_plotAveragePupilResponses( mergedPacketCellArray, twoComponentFitToData, lockedTimeTPUPFitToData, analysisDir)
% Create plots of average evoked pupil responses within and across subjects

%% GET THE AVERAGE PACKETS
normFlag=3; % zero center the initial period, % change units
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, 'lowFreqClean', true, 'aggregateMethod', 'mean');

%% Plot the within-subject data
nSessionsMerged=size(mergedPacketCellArray,2);
nStimTypes=5;

% Set a plot offset for the model components
plotOffset=30;

% First, calculate the subjectScaler
subjectScaler=zeros(1,nSessionsMerged/2);
for ss = 1:nSessionsMerged/2 % loop over subjects
    for ee = 1:2 % loop over experiments, LMS and Mel
        thisPacketIndex=(ee-1)*(nSessionsMerged/2)+ss;
        for mm=1:nStimTypes
            subjectScaler(ss)=subjectScaler(ss)+sum(abs(avgPackets{thisPacketIndex, mm}.response.values));
        end
    end
end
subjectScaler=subjectScaler ./ mean(subjectScaler);

%% Plot the average response across subjects, along with the average fits
for ee = 1:2
    figHandle = figure;
    set(gcf, 'PaperSize', [8.5 11]);
    experimentName=mergedPacketCellArray{(ee-1)*(nSessionsMerged/2)+1}{1}.metaData.projectName;
    switch experimentName
        case 'MaxLMSCRF'
            lineColorBase=[.25 .25 .25];
        case 'MaxMelCRF'
            lineColorBase=[0 0 1];
    end
    for mm=1:nStimTypes
        subplotHandle=subplot(1,nStimTypes,mm);
        for ss = 1:nSessionsMerged/2 % loop over subjects
            thisPacketIndex=(ee-1)*(nSessionsMerged/2)+ss;
            acrossSubResponseMatrix(ss,:)=avgPackets{thisPacketIndex, mm}.response.values ./ subjectScaler(ss);
            acrossSubFitMatrix(ss,:)=twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.values./ subjectScaler(ss);
            acrossSubComponent1Matrx(ss,:)=twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.component1./ subjectScaler(ss);
            acrossSubComponent2Matrx(ss,:)=twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.component2./ subjectScaler(ss);
            acrossSubComponent3Matrx(ss,:)=twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.component3./ subjectScaler(ss);
        end
        acrossSubResponse.values=nanmean(acrossSubResponseMatrix);
        acrossSubResponse.sem=nanstd(acrossSubResponseMatrix)/sqrt(nSessionsMerged/2);
        acrossSubResponse.fit=nanmean(acrossSubFitMatrix);
        acrossSubResponse.component1=nanmean(acrossSubComponent1Matrx);
        acrossSubResponse.component2=nanmean(acrossSubComponent2Matrx);
        acrossSubResponse.component3=nanmean(acrossSubComponent3Matrx);
        acrossSubResponse.timebase=avgPackets{thisPacketIndex, mm}.response.timebase;
        if mm==1
            dataOnly=false;
        else
            dataOnly=true;
        end
        % plot the response and SEM across subjects
        pupilPMEL_PlotEvokedResponse( subplotHandle, ...
            acrossSubResponse.timebase, ...
            100*acrossSubResponse.values, ...
            100*acrossSubResponse.sem, ...
            'dataOnly', dataOnly, ...
            'ylim', [-50 30], ...
            'xAxisAspect', 1, 'yAxisAspect', 2, ...
            'lineColor', lineColorBase, 'plotTitle', experimentName);
        hold on
        % Plot the model fit
        pupilPMEL_PlotEvokedResponse( subplotHandle, ...
            acrossSubResponse.timebase, ...
            acrossSubResponse.fit, ...
            [], ...
            'dataOnly', dataOnly, ...
            'ylim', [-50 30], ...
            'xAxisAspect', 1, 'yAxisAspect', 2, ...
            'lineColor', [1 0 0],'plotTitle', experimentName);
        % Add the three model components
        plotOffset=30;
        pupilPMEL_PlotEvokedResponse( subplotHandle, ...
            acrossSubResponse.timebase, ...
            acrossSubResponse.component1+plotOffset, ...
            [], ...
            'dataOnly', dataOnly, ...
            'ylim', [-50 30], ...
            'lineWidth',1, 'lineStyle','-',...
            'xAxisAspect', 1, 'yAxisAspect', 2, ...
            'lineColor', [1 0.75 0.75],'plotTitle', experimentName);
        pupilPMEL_PlotEvokedResponse( subplotHandle, ...
            acrossSubResponse.timebase, ...
            acrossSubResponse.component2+plotOffset, ...
            [], ...
            'dataOnly', dataOnly, ...
            'ylim', [-50 30], ...
            'lineWidth',1, 'lineStyle','-',...
            'xAxisAspect', 1, 'yAxisAspect', 2, ...
            'lineColor', [1 0.5 0.5],'plotTitle', experimentName);
        pupilPMEL_PlotEvokedResponse( subplotHandle, ...
            acrossSubResponse.timebase, ...
            acrossSubResponse.component3+plotOffset, ...
            [], ...
            'dataOnly', dataOnly, ...
            'ylim', [-50 30], ...
            'lineWidth',0.5, 'lineStyle','-',...
            'xAxisAspect', 1, 'yAxisAspect', 2, ...
            'lineColor', [1 0.25 0.25],'plotTitle', experimentName);
        set(subplotHandle,'FontSize',8);
        hold off
    end
    % Save the plots
    if strcmp([experimentName '_PupilAcrossSubjectResponses'],'MaxMelCRF_PupilAcrossSubjectResponses')
        plotFileName=fullfile(analysisDir, 'figures', 'Figure4b.pdf');
    elseif strcmp([experimentName '_PupilAcrossSubjectResponses'],'MaxLMSCRF_PupilAcrossSubjectResponses')
        plotFileName=fullfile(analysisDir, 'figures', 'Figure4a.pdf');
    else
        plotFileName=fullfile(analysisDir, 'figures', [experimentName '_PupilAcrossSubjectResponses.pdf']);
    end
    pupilPMEL_suptitle(figHandle,['Mean subject pupil responses +-SEM (w subject scaler)']);
    set(figHandle,'Renderer','painters');
    print(figHandle, plotFileName, '-dpdf', '-fillpage');
    close(figHandle);
end

%% Plot CRFs of the pupil response model fit components
contrastLabels={'25','50','100','200','400'};
paramLabels=twoComponentFitToData{1,1}.paramsFit.paramNameCell;
paramYLimLower=[0 0 0 -25 -50 -100];
paramYLimUpper=[300 300 30 0 0 0];
yDirSettings={'normal','normal','normal','reverse','reverse','reverse'};
yLabels={'msecs','msecs','seconds','area [%change * secs]','area [%change * secs]','area [%change * secs]'};

figHandle = figure;
set(gcf, 'PaperSize', [8.5 11]);
for pp=1:length(paramLabels)
    subPlotHandles{pp}=subplot(2,3,pp);
end
for ee = 1:2
    experimentName=mergedPacketCellArray{(ee-1)*(nSessionsMerged/2)+1}{1}.metaData.projectName;
    switch experimentName
        case 'MaxLMSCRF'
            lineColorBase=[.25 .25 .25];
        case 'MaxMelCRF'
            lineColorBase=[0 0 1];
    end
    for mm=1:nStimTypes
        for ss = 1:nSessionsMerged/2 % loop over subjects
            thisPacketIndex=(ee-1)*(nSessionsMerged/2)+ss;
            allParamVals(ss,mm,:)=twoComponentFitToData{thisPacketIndex,mm}.paramsFit.paramMainMatrix;
            % Adjust the gain parameters by the SubjectScaler
            allParamVals(ss,mm,4:6)=allParamVals(ss,mm,4:6) ./ subjectScaler(ss);
        end
    end
    
    meanParamVals=squeeze(mean(allParamVals,1));
    semParamVals=squeeze(std(allParamVals,1))/sqrt(nSessionsMerged/2);
    for pp=1:length(paramLabels)
        plotTitle=[paramLabels{pp} ' subject mean �SEM'];
        % plot the response and SEM across subjects
        pupilPMEL_PlotCRF( subPlotHandles{pp}, 1:5, meanParamVals(:,pp), semParamVals(:,pp), ...
            'xTickLabels',contrastLabels,...
            'xlim',[0 6],...
            'ylim',[paramYLimLower(pp) paramYLimUpper(pp)],...
            'lineColor',lineColorBase,...
            'markerColor',lineColorBase,...
            'errorColor',lineColorBase,...
            'plotTitle',plotTitle,...
            'xLabel','contrast [%]',...
            'yLabel',yLabels{pp},...
            'yDir',yDirSettings{pp});
    end % loop over params
end % loop over experiments
% Save the plots
plotFileName=fullfile(analysisDir, 'figures', 'Figure4c_FigureS10.pdf'); %PupilFitParamCRFs
pupilPMEL_suptitle(figHandle,['Mean pupil fit params +-SEM (w subject scaler for amplitudes)']);
set(figHandle,'Renderer','painters');
print(figHandle, plotFileName, '-dpdf', '-fillpage');
close(figHandle);


%% Plot CRFs of the pupil response model fit components using the locked temporal components
contrastLabels={'25','50','100','200','400'};
paramLabels=lockedTimeTPUPFitToData{1,1}.paramsFit.paramNameCell;
paramYLimLower=[0 0 0 -25 -50 -100];
paramYLimUpper=[300 300 30 0 0 0];
yDirSettings={'normal','normal','normal','reverse','reverse','reverse'};
yLabels={'msecs','msecs','seconds','area [%change * secs]','area [%change * secs]','area [%change * secs]'};

figHandle = figure;
set(gcf, 'PaperSize', [8.5 11]);
for pp=1:length(paramLabels)
    subPlotHandles{pp}=subplot(2,3,pp);
end
for ee = 1:2
    experimentName=mergedPacketCellArray{(ee-1)*(nSessionsMerged/2)+1}{1}.metaData.projectName;
    switch experimentName
        case 'MaxLMSCRF'
            lineColorBase=[.25 .25 .25];
        case 'MaxMelCRF'
            lineColorBase=[0 0 1];
    end
    for mm=1:nStimTypes
        for ss = 1:nSessionsMerged/2 % loop over subjects
            thisPacketIndex=(ee-1)*(nSessionsMerged/2)+ss;
            allParamVals(ss,mm,:)=lockedTimeTPUPFitToData{thisPacketIndex,mm}.paramsFit.paramMainMatrix;
            % Adjust the gain parameters by the SubjectScaler
            allParamVals(ss,mm,4:6)=allParamVals(ss,mm,4:6) ./ subjectScaler(ss);
        end
    end
    
    meanParamVals=squeeze(mean(allParamVals,1));
    semParamVals=squeeze(std(allParamVals,1))/sqrt(nSessionsMerged/2);
    for pp=1:length(paramLabels)
        plotTitle=[paramLabels{pp} ' subject mean �SEM'];
        % plot the response and SEM across subjects
        pupilPMEL_PlotCRF( subPlotHandles{pp}, 1:5, meanParamVals(:,pp), semParamVals(:,pp), ...
            'xTickLabels',contrastLabels,...
            'xlim',[0 6],...
            'ylim',[paramYLimLower(pp) paramYLimUpper(pp)],...
            'lineColor',lineColorBase,...
            'markerColor',lineColorBase,...
            'errorColor',lineColorBase,...
            'plotTitle',plotTitle,...
            'xLabel','contrast [%]',...
            'yLabel',yLabels{pp},...
            'yDir',yDirSettings{pp});
    end % loop over params
end % loop over experiments
% Save the plots
plotFileName=fullfile(analysisDir, 'figures', 'FigureReviewer2_LockedPupilParams.pdf'); %PupilFitParamCRFs
pupilPMEL_suptitle(figHandle,['Mean pupil fit params with locked temporal params +-SEM (w subject scaler for amplitudes)']);
set(figHandle,'Renderer','painters');
print(figHandle, plotFileName, '-dpdf', '-fillpage');
close(figHandle);


%% Plot evoked responses for each subject / stimulus
for ee = 1:2 % loop over experiments, LMS and Mel
    figHandle = figure;
    set(gcf, 'PaperSize', [8.5 11]);
    for ss = 1:nSessionsMerged/2 % loop over subjects
        
        thisPacketIndex=(ee-1)*(nSessionsMerged/2)+ss;
        % Identify the stimulus direction (LMS or Mel)
        experimentName=mergedPacketCellArray{thisPacketIndex}{1}.metaData.projectName;
        lineColorBase=[.25 .25 .25];
        switch experimentName
            case 'MaxLMSCRF'
                lineColorBase=[.25 .25 .25];
            case 'MaxMelCRF'
                lineColorBase=[0 0 1];
        end
                
        for mm = 1:nStimTypes
            if mm==1
                dataOnly=false;
            else
                dataOnly=true;
            end
            % Plot the response and SEM error
            subplotHandle=subplot(nSessionsMerged/2,nStimTypes,(ss-1)*nStimTypes+mm);
            plotTitle=mergedPacketCellArray{thisPacketIndex}{1}.metaData.subjectName;
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                100*avgPackets{thisPacketIndex, mm}.response.values, ...
                100*avgPackets{thisPacketIndex, mm}.response.sem, ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 30], ...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', lineColorBase, 'plotTitle', plotTitle);
            hold on
            % Plot the model fit
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.values, ...
                [], ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 30], ...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', [1 0 0],'plotTitle', plotTitle);
            % Add the three model components
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.component1+plotOffset, ...
                [], ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 30], ...
                'lineWidth',1, 'lineStyle','-',...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', [1 0.75 0.75],'plotTitle', plotTitle);
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.component2+plotOffset, ...
                [], ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 30], ...
                'lineWidth',1, 'lineStyle','-',...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', [1 0.5 0.5],'plotTitle', plotTitle);
            pupilPMEL_PlotEvokedResponse( subplotHandle, ...
                avgPackets{thisPacketIndex, mm}.response.timebase, ...
                twoComponentFitToData{thisPacketIndex,mm}.modelResponseStruct.component3+plotOffset, ...
                [], ...
                'dataOnly', dataOnly, ...
                'ylim', [-50 30], ...
                'lineWidth',0.5, 'lineStyle','-',...
                'xAxisAspect', 1, 'yAxisAspect', 2, ...
                'lineColor', [1 0.25 0.25],'plotTitle', plotTitle);
            set(subplotHandle,'FontSize',8);
        end % loop over stim types
        
    end % loop over subjects
    
    % Save the plots
    if strcmp ([experimentName '_PupilTrialMeanResponses'],'MaxLMSCRF_PupilTrialMeanResponses')
        plotFileName=fullfile(analysisDir, 'figures', 'FigureS8a.pdf');
    elseif strcmp ([experimentName '_PupilTrialMeanResponses'],'MaxMelCRF_PupilTrialMeanResponses')
        plotFileName=fullfile(analysisDir, 'figures', 'FigureS8b.pdf');
    else
        plotFileName=fullfile(analysisDir, 'figures', [experimentName '_PupilTrialMeanResponses.pdf']);
    end
    pupilPMEL_suptitle(figHandle,[experimentName ' - Pupil responses']);
    set(figHandle,'Renderer','painters');
    print(figHandle, plotFileName, '-dpdf', '-fillpage');
    close(figHandle);
    
end % loop over session sets (Mel and LMS)

end % function