%% matching runs from MRI data to runs from pupil data
% mergedPacketCellArray{1,1}{1,1}.metaData.responseFile corresponds to packetCellArray{1,1}.metaData.stimulusFile

% LMS
% subject 1: HERO_asb1/060816
% subject 2: HERO_aso1/060116
% subject 3: HERO_gka1/060616
% subject 4: HERO_mxs1/062816

% Mel
% subject 1: HERO_asb1/060716
% subject 2: HERO_aso1/053116
% subject 3: HERO_gka1/060216
% subject 4: HERO_mxs1/060916

[~, userName] = system('whoami');
userName = strtrim(userName);
dropboxAnalysisDir = ...
    fullfile('/Users', userName, ...
    '/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/pupilMelanopsinMRIAnalysis');
subAnalysisDirectory = 'BOLD';


for ss = 1:8;
    for analysis = 1:2;
        rsqLag{ss,analysis} = [];
    end
end

for analysis = 1:2; % low frequency component, total response
    if analysis == 1;
        pupilType = 'Low Frequency Component';
        saveAs = 'lowFreq';
    elseif analysis == 2;
        pupilType = 'Unfiltered Response';
        saveAs = 'unfiltered';
    end
    for ss = 1:3;
        plotFig = figure;
        for pp = 1:2;
            if pp == 1; % LMS
                boldPacket = LMSpacketCellArray;
                adjustment = 0;
            elseif pp == 2; % mel
                boldPacket = MELpacketCellArray;
                adjustment = 4; % mel subjects are 4 later (subject 1 in LMS corresponds to subject 5 in mel)
            end
            
            bold = [];
            pupil = [];
            for rr = 1:length(mergedPacketCellArray{ss+adjustment}); % loop over runs
                
                
                %% resample pupil data to frequency of MR data
                % MR data -> TR = 0.8 s
                mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.timebaseResampled = boldPacket{ss,rr}.response.timebase;
                mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.lowFreqComponentResampled = [];
                
                % perform convolution
                if analysis == 1;
                    responseStruct.values = mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.metaData.lowFreqComponent;
                elseif analysis == 2;
                    responseStruct.values = mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.values;
                end
                responseStruct.timebase = mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.timebase;
                responseStruct.metaData = mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.metaData;
                % convolution
                temporalFit = tfeIAMP('verbosity','none');
                convolutionStruct = temporalFit.applyKernel(responseStruct, boldPacket{ss,rr}.kernel);
                
                
                
                for x = 1:length(boldPacket{ss,rr}.response.timebase);
                    index = find(mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.timebase==round(mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.timebaseResampled(x)));
                    % for no convolution, use the next two lines:
                    %mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.valuesResampled(x) = mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.values(index);
                    %mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.metaData.lowFreqComponentResampled(x) = mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.metaData.lowFreqComponent(index);
                    
                    % for convolution, use the next two lines:
                    mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.valuesResampled(x) = convolutionStruct.values(index);
                    mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.metaData.lowFreqComponentResampled(x) = convolutionStruct.values(index);
                    
                    
                end
                
                bold = [bold, boldPacket{ss,rr}.response.values];
                
                
                if analysis == 1;
                    pupil = [pupil, mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.metaData.lowFreqComponentResampled];
                elseif analysis == 2;
                    pupil = [pupil, mergedPacketCellArray{1,ss+adjustment}{1,rr}.response.valuesResampled];
                end
                
            end
            
            % test auto-correlations
            %pupil = bold;
            % deleting the first 5 TRs, so we should have a 4 second lag
            %pupil = pupil(1,6:length(pupil));
            %pupil(1,length(pupil)+1:length(pupil)+5) = NaN;
            
            % shift pupil data in time according to lagIndex
            
            lagRange = -20:50;
            for lagIndex = lagRange;
                laggedPupilTimeSeries = [];
                for x = 1:length(pupil);
                    if x-lagIndex > 0 && (x-lagIndex) < (length(pupil)+1);
                        laggedPupilTimeSeries(x) = pupil(x-lagIndex);
                        %laggedPupilTimeSeries(x) = mergedPacketCellArray{1,ss}{1,rr}.response.valuesResampled(x-lagIndex);
                        if lagIndex > 0;
                            laggedPupilTimeSeries(1:lagIndex) = NaN;
                        end
                        if lagIndex < 0;
                            laggedPupilTimeSeries(length(pupil)+1+lagIndex:length(pupil)) = NaN;
                        end
                    end
                end
                % calculate R2 for a given lagIndex
                x = laggedPupilTimeSeries;
                y = bold;
                xnan = isnan(x);
                hits = 0;
                for xx = 1:length(x);
                    if xnan(xx) == 1;
                        hits = hits+1;
                        x(xx-(hits-1)) = [];
                        y(xx-(hits-1)) = [];
                        
                    end
                end
                
                %x = [ones(1,length(x)); x];
                %x = x';
                %y = y';
                %[B,BINT,R,RINT,STATS] = regress(y,x);
                %slope = B(1);
                %rsq = STATS(1);
                rsq = corr2(x,y);
                rsqLag{(ss+adjustment), analysis}(1,(find(lagRange==lagIndex))) = lagIndex;
                rsqLag{(ss+adjustment), analysis}(2,(find(lagRange==lagIndex))) = rsq;
            end
            subplot(2,1,pp);
            if pp == 1;
                plot(rsqLag{(ss+adjustment),analysis}(1,:)*0.8,rsqLag{(ss+adjustment),analysis}(2,:), 'Color', 'r')
                title(['Subject: ',  mergedPacketCellArray{ss}{1}.metaData.subjectName, ' LMS ', pupilType]);
                
            elseif pp == 2;
                plot(rsqLag{(ss+adjustment),analysis}(1,:)*0.8,rsqLag{(ss+adjustment),analysis}(2,:), 'Color', 'b')
                title(['Subject: ',  mergedPacketCellArray{ss}{1}.metaData.subjectName, ' Melanopsin ', pupilType]);
                
            end
            
        end % end loop over project
        outDir = fullfile(dropboxAnalysisDir,subAnalysisDirectory);
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        saveas(plotFig, fullfile(outDir, [mergedPacketCellArray{ss}{1}.metaData.subjectName, '_', saveAs, '.png']), 'png');
        close(plotFig);
    end % end loop over subjects
end % loop over analysis types


% create average plots
for analysis = 1:2;
    for x = lagRange;
        rsqLagCollapsed{1,(find(lagRange==x))} = x;
        rsqLagCollapsed{2,(find(lagRange==x))} = [];
        rsqLagCollapsed{3,(find(lagRange==x))} = [];
        
    end
    
    for ss = 1:3;
        
        for x = lagRange;
            rsqLagCollapsed{analysis+1,(find(lagRange==x))} = [rsqLagCollapsed{2,(find(lagRange==x))}, rsqLag{ss,analysis}(2,(find(lagRange==x)))];
            rsqLagCollapsed{analysis+1,(find(lagRange==x))} = [rsqLagCollapsed{2,(find(lagRange==x))}, rsqLag{ss+4,analysis}(2,(find(lagRange==x)))];
            
        end
    end
    for x = lagRange;
        meanRsqLag(1,find(lagRange==x)) = x;
        meanRsqLag(analysis+1,find(lagRange==x))= mean(rsqLagCollapsed{analysis+1,(find(lagRange==x))});
    end
end


plotFig = figure;
plot(meanRsqLag(1,:),meanRsqLag(2,:));
title('Lag Plot on Low Frequency Component')
saveas(plotFig, fullfile(outDir, ['averageLowFreq.png']), 'png');
close(plotFig);

plotFig = figure;
plot(meanRsqLag(1,:),meanRsqLag(3,:));
title('Lag Plot on Unfiltered Pupil Data')
saveas(plotFig, fullfile(outDir, ['averageUnfiltered.png']), 'png');
close(plotFig);
