function [ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, filterStatus )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% comment out or change here to change filter status
filterStatus = 'filtered';

% Set some parameters we need
normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;

% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=3; % zero center the initial period, change units
splitParams.normalizationWindowMsecs=100; % define the size of the norm window

NSessionsMerged=size(mergedPacketCellArray,2);

%% Assemble the different stimTypes
for ss = 1:NSessionsMerged
    % Make an 'accumulator'
    NStimTypes = 6;
    for mm = 1:NStimTypes
        accumStimTypesResp{ss, mm} = [];
        accumStimTypesStim{ss, mm} = [];
    end
    
    NRunsTotal = length(mergedPacketCellArray{ss});
    for rr = 1:NRunsTotal
        
        % grab a packet that corresponds to a run for a given subject
        theRunPacket=mergedPacketCellArray{1,ss}{rr};
        
        % if indicated, filter out the low frequency component
        if strcmp(filterStatus, 'filtered');
            
            lowFreqComponent=theRunPacket.response.metaData.lowFreqComponent;
            lowFreqComponent=lowFreqComponent-mean(lowFreqComponent);
            theRunPacket.response.values = ...
                theRunPacket.response.values - lowFreqComponent;
        end
        
        
        for ii = 1:size(mergedPacketCellArray{ss}{rr}.stimulus.values,1)
            splitParams.instanceIndex = ii;
            
            % grab the packet for this particular instance
            singlePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
            
            % identify contrast associated with that instance
            contrast = theRunPacket.stimulus.metaData.stimTypes(ii);
            accumStimTypesResp{ss, contrast} = [accumStimTypesResp{ss, contrast} ; singlePacket.response.values];
            accumStimTypesStim{ss, contrast} = [accumStimTypesStim{ss, contrast} ; singlePacket.stimulus.values];
        end
    end
end

%% Make average packets per subject
for ss = 1:NSessionsMerged
    for mm = 1:NStimTypes
        thisPacket = [];
        thisPacket.packetType = 'pupil';
        thisPacket.sessionDir = '';
        thisPacket.stimulusFile = mergedPacketCellArray{ss}{1}.metaData.stimulusFile;
        thisPacket.responseFile = mergedPacketCellArray{ss}{1}.metaData.responseFile;
        thisPacket.respValues =  nanmean(accumStimTypesResp{ss, mm});
        thisPacket.respTimeBase = [0:12999];
        thisPacket.stimValues = max(accumStimTypesStim{ss, mm});
        thisPacket.stimTimeBase = [0:12999];
        thisPacket.stimMetaData.stimTypes = mergedPacketCellArray{ss}{1}.stimulus.metaData.stimTypes;
        thisPacket.stimMetaData.stimLabels = mergedPacketCellArray{ss}{1}.stimulus.metaData.stimLabels;
        avgPackets{ss, mm} = makePacket(thisPacket);
    end
end

end