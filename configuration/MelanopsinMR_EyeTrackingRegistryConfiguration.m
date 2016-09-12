%% OLFlickerSensitivityRegistryConfig
%
% Declare the toolboxes we need for the IBIOColorDetect project and
% write them into a JSON file.  This will let us use the ToolboxToolbox to
% deliver unto us the perfect runtime environment for this project.
%
% 2016 benjamin.heasly@gmail.com

% Clear
clear;

%% Declare some toolboxes we want.
config = [ ...  
    tbToolboxRecord( ...
    'name', 'BrainardLabToolbox', ...
    'type', 'git', ...
    'url', 'https://github.com/brainardlab/BrainardLabToolbox.git'), ...
    tbToolboxRecord( ...
    'name', 'MelanopsinMR_EyeTracking', ...
    'type', 'git', ...
    'url', 'https://github.com/gkaguirrelab/MelanopsinMR_EyeTracking.git'), ...
    tbToolboxRecord( ...
    'name', 'mriTemporalFitting', ...
    'type', 'git', ...
    'url', 'https://github.com/gkaguirrelab/mriTemporalFitting.git')];
    
%% Write the config to a JSON file.
configPath = 'MelanopsinMR_EyeTrackingConfig.json';
tbWriteConfig(config, 'configPath', configPath);
