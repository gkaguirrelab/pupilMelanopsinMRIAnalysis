function [ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Set some parameters we need
normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;

NSessionsMerged=size(mergedPacketCellArray,1);

%% Assemble the different stimTypes
for ss = 1:NSessionsMerged
    % Make an 'accumulator'
    NStimTypes = 6;
    for mm = 1:NStimTypes
        accumStimTypesResp{ss, mm} = [];
        accumStimTypesStim{ss, mm} = [];
    end
    
    NRunsTotal = length(mergedPacketCellArray{ss});
    for ii = 1:NRunsTotal
        % Find the stimulus onsets so that we can align the data to it. We
        % do that by finding a [0 1] edge from a difference operator.
        tmp = sum(mergedPacketCellArray{ss}{ii}.stimulus.values);
        tmp2 = diff(tmp);
        tmp2(tmp2 < 0) = 0;
        tmp2(tmp2 > 0) = 1;
        stimOnsets = strfind(tmp2, [0 1]);
        
        % Get the number of segments from the stimulus onsets
        NSegments = length(stimOnsets);
        t = (0:extractionDurInd)/1000;
        for jj = 1:length(stimOnsets)
            if (stimOnsets(jj)+extractionDurInd) <= length(mergedPacketCellArray{ss}{ii}.response.values)
                idxToExtract = stimOnsets(jj):(stimOnsets(jj)+extractionDurInd);
            else
                idxToExtract = stimOnsets(jj):length(mergedPacketCellArray{ss}{ii}.response.values);
            end
            thisPacket.packetType = 'pupil';
            thisPacket.stimulusFile = mergedPacketCellArray{ss}{ii}.metaData.stimulusFile;
            thisPacket.responseFile = mergedPacketCellArray{ss}{ii}.metaData.responseFile;
            thisPacket.respValues = mergedPacketCellArray{ss}{ii}.response.values(idxToExtract);
            % Normalize the pupil data
            thisPacket.respValues = (thisPacket.respValues - nanmean(thisPacket.respValues(1:normalizationDurInd)))./nanmean(thisPacket.respValues(1:normalizationDurInd));
            thisPacket.respTimeBase = mergedPacketCellArray{ss}{ii}.response.timebase(idxToExtract);
            thisPacket.respTimeBase = thisPacket.respTimeBase-thisPacket.respTimeBase(1);
            thisPacket.stimValues = mergedPacketCellArray{ss}{ii}.stimulus.values(jj, idxToExtract);
            thisPacket.stimTimeBase = mergedPacketCellArray{ss}{ii}.stimulus.timebase(idxToExtract);
            thisPacket.stimTimeBase =  thisPacket.stimTimeBase - thisPacket.stimTimeBase(1);
            thisPacket.stimMetaData.stimTypes = mergedPacketCellArray{ss}{ii}.stimulus.metaData.stimTypes(jj);
            thisPacket.stimMetaData.stimLabels = mergedPacketCellArray{ss}{ii}.stimulus.metaData.stimLabels;
            
            % Could make packets here for each event, but not doing it...
            % Just making the 'accumStimTypes' variable for now.
            accumStimTypesResp{ss, thisPacket.stimMetaData.stimTypes} = [accumStimTypesResp{ss, thisPacket.stimMetaData.stimTypes} ; thisPacket.respValues];
            accumStimTypesStim{ss, thisPacket.stimMetaData.stimTypes} = [accumStimTypesStim{ss, thisPacket.stimMetaData.stimTypes} ; thisPacket.stimValues];
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
        thisPacket.respTimeBase = mergedPacketCellArray{ss}{1}.response.timebase(idxToExtract);
        thisPacket.respTimeBase = thisPacket.respTimeBase-thisPacket.respTimeBase(1);
        thisPacket.stimValues = max(accumStimTypesStim{ss, mm});
        thisPacket.stimTimeBase = mergedPacketCellArray{ss}{1}.stimulus.timebase(idxToExtract);
        thisPacket.stimTimeBase = thisPacket.stimTimeBase-thisPacket.stimTimeBase(1);
        thisPacket.stimMetaData.stimTypes = mergedPacketCellArray{ss}{ii}.stimulus.metaData.stimTypes(jj);
        thisPacket.stimMetaData.stimLabels = mergedPacketCellArray{ss}{ii}.stimulus.metaData.stimLabels;
        avgPackets{ss, mm} = makePacket(thisPacket);
    end
end

end

