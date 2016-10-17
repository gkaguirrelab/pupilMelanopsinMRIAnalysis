function [ avgPackets ] = pupilPMEL_makeAverageResponsePackets( mergedPacketCellArray, normFlag, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Parse vargin for options passed here
%
% Setting 'KeepUmatched' to true means that we can pass the varargin{:})
% along from a calling routine without an error here, if the key/value
% pairs recognized by the calling routine are not needed here.
p = inputParser; p.KeepUnmatched = true;
p.addRequired('mergedPacketCellArray',@iscell);
p.addRequired('normFlag',@isnumeric);
p.addParameter('lowFreqClean','false',@islogical);
p.addParameter('aggregateMethod','mean',@ischar);
p.parse(mergedPacketCellArray, normFlag,varargin{:});


normalizationTimeSecs = 0.1;
normalizationDurInd = normalizationTimeSecs*1000-1;
extractionTimeSecs = 13;
extractionDurInd = extractionTimeSecs*1000-1;



% define the split params
splitParams.instanceIndex=[]; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normFlag=normFlag; % zero center the initial period, % change units
splitParams.normalizationWindowMsecs=100; % define the size of the norm window

NStimTypes=6;

%% Loop through the sessions
NSessions=size(mergedPacketCellArray,2);

% pre-allocate the avgPackets cell array
avgPackets=cell(NSessions,NStimTypes);

for ss = 1:NSessions
    
    NRuns=size(mergedPacketCellArray{ss},2);
    
    for rr = 1:NRuns
        
        % grab a packet that corresponds to a run for a given subject
        theRunPacket=mergedPacketCellArray{ss}{rr};
        
        % Adjust the response in the packet to remove the low-frequency
        % variation in pupil response
        if p.Results.lowFreqClean
            lowFreqComponent=theRunPacket.response.metaData.lowFreqComponent;
            lowFreqComponent=lowFreqComponent-mean(lowFreqComponent);
            theRunPacket.response.values = ...
                theRunPacket.response.values - lowFreqComponent;
        end
        
        for ii = 1:size(theRunPacket.stimulus.values,1)
            
            % Put the instance index into the splitParams
            splitParams.instanceIndex = ii;
            
            % grab the packet for this particular instance
            singlePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
            
            % identify the stim type for this instance
            thisStimType = theRunPacket.stimulus.metaData.stimTypes(ii);
            
            % Add to (or create) this stimType in the avgPackets
            if isempty(avgPackets{ss,thisStimType})
                avgPackets{ss,thisStimType}=singlePacket;
                stimTypeCounter(thisStimType)=1;
            else
                stimTypeCounter(thisStimType)=stimTypeCounter(thisStimType)+1;
                avgPackets{ss,thisStimType}.response.values(stimTypeCounter(thisStimType),:)= ...
                    singlePacket.response.values;
                avgPackets{ss,thisStimType}.stimulus.values(stimTypeCounter(thisStimType),:)= ...
                    singlePacket.stimulus.values;
            end % check for the existence of this stimType entry
        end % loop over instances
    end % loop over runs
    
    % Take the nanmean of the matrix of values for each simulus/response
    for cc=1:length(stimTypeCounter)
        if ~stimTypeCounter(cc)==0
            switch (p.Results.aggregateMethod)
                case 'mean'
                    avgPackets{ss,cc}.response.values = nanmean(avgPackets{ss,cc}.response.values);
                    avgPackets{ss,cc}.stimulus.values = nanmean(avgPackets{ss,cc}.stimulus.values);
                case 'median'
                    avgPackets{ss,cc}.response.values = nanmedian(avgPackets{ss,cc}.response.values);
                    avgPackets{ss,cc}.stimulus.values = nanmedian(avgPackets{ss,cc}.stimulus.values);
                otherwise
                    error('Unknown aggregation method requested');
            end % switch aggregateMethod
        end % check that an instances was found
    end % loop over stimTypes
end % loop over sessions

end % function