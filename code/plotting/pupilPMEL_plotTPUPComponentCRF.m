function [ ] = pupilPMEL_plotTPUPComponentCRF( twoComponentFitToData, dropboxAnalysisDir )
% Create plots of average evoked pupil responses within and across subjects


nStimTypes=5;
nSubjects=4;

    figHandle = figure;
    set(gcf, 'PaperSize', [8.5 11]);
    sumComponents=nan(2,nSubjects,nStimTypes);
    ratioComponents=nan(2,nSubjects,nStimTypes);

    
for ee = 1:2 % loop over experiments, LMS and Mel
    for ss = 1:nSubjects % loop over subjects        
        for mm = 1:nStimTypes
        thisPacketIndex=(ee-1)*nSubjects+ss;
        paramMainMatrix=twoComponentFitToData{thisPacketIndex,mm}.paramsFit.paramMainMatrix;
        sumComponents(ee,ss,mm)=sum(paramMainMatrix(6));
        ratioComponents(ee,ss,mm)=paramMainMatrix(5)/paramMainMatrix(6);
        end % loop over stim types
         end % loop over subjects
end % loop over experiments

sumCom=squeeze(nanmedian(sumComponents,2));
ratioCom=squeeze(nanmedian(ratioComponents,2));
plot(sumCom(1,:),'-ok');
hold on
plot(sumCom(2,:),'-ob');
figure
plot(ratioCom(1,:),'-ok');
hold on
plot(ratioCom(2,:),'-ob');


%     % Save the plots
%     plotFileName=fullfile(dropboxAnalysisDir, 'Figures', [experimentName '_PupilTrialMeanResponses.pdf']);
%     pupilPMEL_suptitle(figHandle,[experimentName ' - Pupil responses']);
%     set(gca,'FontSize',6);
%     set(figHandle,'Renderer','painters');
%     print(figHandle, plotFileName, '-dpdf', '-fillpage');
%     close(figHandle);
    

end % function