[~, userName] = system('whoami');
userName = strtrim(userName);

dropboxDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data');
sessDirs = {'MelanopsinMRMaxLMSCRF/HERO_asb1/060816' ...
    'MelanopsinMRMaxLMSCRF/HERO_aso1/060116' ...
    'MelanopsinMRMaxLMSCRF/HERO_gka1/060616' ...
    'MelanopsinMRMaxLMSCRF/HERO_mxs1/062416' ...
    'MelanopsinMRMaxLMSCRF/HERO_mxs1/062816' ...
    'MelanopsinMRMaxMelCRF/HERO_asb1/060716' ...
    'MelanopsinMRMaxMelCRF/HERO_aso1/053116' ...
    'MelanopsinMRMaxMelCRF/HERO_gka1/060216' ...
    'MelanopsinMRMaxMelCRF/HERO_mxs1/060916' ...
    'MelanopsinMRMaxMelCRF/HERO_mxs1/061016'};

for ss = 1:length(sessDirs)
    Data_Per_Segment = [];
    uniqueCombos = [];
    packets = makePackets(fullfile(dropboxDir, sessDirs{ss}), 'pupil', 'pupil', [], false);
    %%
    NRuns = length(packets);
    allData = [];
    allIndices = [];
    %
    for ii = 1:NRuns;
        
        tmp = sum(packets{ii}.stimulus.values);
        tmp2 = diff(tmp);
        tmp2(tmp2 < 0) = 0;
        tmp2(tmp2 > 0) = 1;
        
        stimOnsets = strfind(tmp2, [0 1]);
        
        %
        NSegments = length(stimOnsets);
        durIdx = 12999;
        t = (0:durIdx)/1000;
        for jj = 1:length(stimOnsets)
            %figure
            %packets{1}.stimulus.metaData.params.theContrastRelMaxIndices(jj)
            if (stimOnsets(jj)+durIdx) <= length(packets{ii}.response.values)
                Data_Per_Segment(:, jj) = (packets{ii}.response.values(stimOnsets(jj):(stimOnsets(jj)+durIdx)));
            else
                Data_Per_Segment(:, jj) = (packets{ii}.response.values(stimOnsets(jj):end));
            end
        end
        
        normalizationTimeSecs = 0.1;
        normalizationDurId = normalizationTimeSecs*1000-1;
        
        % Normalize by the dead time
        for jj  = 1:NSegments
            Data_Per_Segment(:, jj) = (Data_Per_Segment(:, jj) - nanmean(Data_Per_Segment(1:normalizationDurId, jj)))./nanmean(Data_Per_Segment(1:normalizationDurId, jj));
        end
        allData = [allData Data_Per_Segment];
        
        %
        [uniqueCombos{ii}, ~, idx] = unique([packets{ii}.stimulus.metaData.params.theDirections' packets{ii}.stimulus.metaData.params.theContrastRelMaxIndices'], 'rows');
        allIndices = [allIndices idx];
    end
    %%
    
    % Let's just make sure that in every run we have the same conditions
    for ii = 1:NRuns
        logTest(ii) = all(all(uniqueCombos{ii} == uniqueCombos{1}));
    end
    if all(logTest)
        uniqueCombos1 = uniqueCombos{1};
    else
        error('This session is heterogenous');
    end
    
    
    %%
    NCombos = size(uniqueCombos1, 1);
    plotFig = figure;
    for jj = 1:NCombos
        subplot(1, NCombos, jj);
        data = allData(:, allIndices == jj);
        data_mean = nanmean(data, 2);
        data_sem = nanstd(data, [], 2)/sqrt(size(data, 2)); hold on;
        plot([t(1) t(end)], [0 0], '-k');
        shadedErrorBar(t, data_mean, data_sem);
        xlim([t(1) t(end)]);
        ylim([-0.5 0.5]);
        pbaspect([1 1 1]);
        xlabel('Time [s]');
        ylabel('Amplitude [%]');
        if uniqueCombos1(jj, 1) == 1
            title([num2str(100*packets{ii}.stimulus.metaData.params.theContrastMax*packets{ii}.stimulus.metaData.params.theContrastsPct(jj), '%g') '%'])
        else
            title('Attention task');
        end
    end
    adjustPlot(plotFig)
    fileName = strrep(sessDirs{ss}, '/', '_');
    
    set(plotFig, 'PaperPosition', [0 0 13 3]);
    set(plotFig, 'PaperSize', [13 3]);
    saveas(plotFig, ['~/Desktop/' fileName '.png'], 'png');
    %pause;
end