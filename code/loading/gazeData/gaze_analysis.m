% Load the packet
%load('/Users/spitschan/Dropbox (Aguirre-Brainard Lab)/MELA_analysis/pupilMelanopsinMRIAnalysis/packetCache/maxMelLSM_CRF_Pupil_f27f85474cfa53c837113ec15f56c17c.mat')
for obs = 1:8
    gapMs = 4000;
    NRuns = length(mergedPacketCellArray{obs});
    NStimuli = 5;
    if obs < 5
        theStimulusLabels = {'MelanopsinMRMaxLMS_PulseMaxLMS_3s_CRF16sSegment_25%' ...
            'MelanopsinMRMaxLMS_PulseMaxLMS_3s_CRF16sSegment_50%' ...
            'MelanopsinMRMaxLMS_PulseMaxLMS_3s_CRF16sSegment_100%' ...
            'MelanopsinMRMaxLMS_PulseMaxLMS_3s_CRF16sSegment_200%' ...
            'MelanopsinMRMaxLMS_PulseMaxLMS_3s_CRF16sSegment_400%'}
    else
        theStimulusLabels = {'MelanopsinMRMaxMel_PulseMaxMel_3s_CRF16sSegment_25%' ...
            'MelanopsinMRMaxMel_PulseMaxMel_3s_CRF16sSegment_50%' ...
            'MelanopsinMRMaxMel_PulseMaxMel_3s_CRF16sSegment_100%' ...
            'MelanopsinMRMaxMel_PulseMaxMel_3s_CRF16sSegment_200%' ...
            'MelanopsinMRMaxMel_PulseMaxMel_3s_CRF16sSegment_400%'}
    end
    
    % Populate an empty vector
    for kk = 1:NStimuli
        gazeX_aggregateVarOn{kk} = [];
        gazeX_aggregateVarOff{kk} = [];
        gazeY_aggregateVarOn{kk} = [];
        gazeY_aggregateVarOff{kk} = [];
    end
    
    for ii = 1:NRuns
        % Determine which trial types we care about
        % validTrialTypes = find(~cellfun('isempty',strfind(mergedPacketCellArray{obs}{ii}.stimulus.metaData.stimLabels,'CRF')));
        
        % Construct a master cell array of all stimulus labels
        allStimLabels = {mergedPacketCellArray{obs}{ii}.stimulus.metaData.stimLabels{mergedPacketCellArray{obs}{ii}.stimulus.metaData.stimTypes}};
        
        
        % Iterate over 16 s segments
        for mm = 1:size(mergedPacketCellArray{obs}{ii}.stimulus.values, 1)
            idx_stimOn = find(mergedPacketCellArray{obs}{ii}.stimulus.values(mm, :) == 1);
            idx_stimOff = idx_stimOn(end)+4000:idx_stimOn(end)+6000;
            gazeX_theVarOn(ii, mm) = var(mergedPacketCellArray{obs}{ii}.response.gazeX(idx_stimOn));
            gazeX_theVarOff(ii, mm) = var(mergedPacketCellArray{obs}{ii}.response.gazeX(idx_stimOff));
            gazeY_theVarOn(ii, mm) = var(mergedPacketCellArray{obs}{ii}.response.gazeY(idx_stimOn));
            gazeY_theVarOff(ii, mm) = var(mergedPacketCellArray{obs}{ii}.response.gazeY(idx_stimOff));
        end
        
        % Iterate over the stimulus types
        for kk = 1:NStimuli
            theStimulusLabels{kk}
            theIdx = find(ismember(allStimLabels, theStimulusLabels{kk}))
            gazeX_aggregateVarOn{kk} = [gazeX_aggregateVarOn{kk} gazeX_theVarOn(ii, theIdx)];
            gazeX_aggregateVarOff{kk} = [gazeX_aggregateVarOff{kk} gazeX_theVarOff(ii, theIdx)];
            gazeY_aggregateVarOn{kk} = [gazeY_aggregateVarOn{kk} gazeY_theVarOn(ii, theIdx)];
            gazeY_aggregateVarOff{kk} = [gazeY_aggregateVarOff{kk} gazeY_theVarOff(ii, theIdx)];
        end
    end
    
    % Plot the data
    figure;
    subplot(1, 2, 1);
    for kk = 1:NStimuli
        errorbar(kk-0.1, mean(gazeX_aggregateVarOff{kk}), std(gazeX_aggregateVarOff{kk}), 'sk', 'MarkerFaceColor', 'k'); hold on;
        errorbar(kk+0.1, mean(gazeX_aggregateVarOn{kk}), std(gazeX_aggregateVarOn{kk}), 'sr', 'MarkerFaceColor', 'r');
    end
    ylabel('X variance'); xlabel('Stimulus contrast');
    pbaspect([1 1 1]);
    
    subplot(1, 2, 2);
    for kk = 1:NStimuli
        errorbar(kk-0.1, mean(gazeY_aggregateVarOff{kk}), std(gazeY_aggregateVarOff{kk}), 'sk', 'MarkerFaceColor', 'k'); hold on;
        errorbar(kk+0.1, mean(gazeY_aggregateVarOn{kk}), std(gazeY_aggregateVarOn{kk}), 'sr', 'MarkerFaceColor', 'r');
    end
    ylabel('Y variance'); xlabel('Stimulus contrast');
    pbaspect([1 1 1]);
end