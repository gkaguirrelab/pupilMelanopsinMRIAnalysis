[simulatedMergedPacketCellArray] = pupilPMEL_makeSimulatedData()

%subject 1, run 4, instance 3 is one of the instances affected.
%specifically applies to data that has been simulated with carryOver.
ss = 1;
rr = 4;
ii = 3;

% apply the relevant parameters, directly from fitIAMP
splitParams.instanceIndex=ii; % will hold the instance index
splitParams.splitDurationMsecs=13000; % Grab 13 second windows
splitParams.normalizationWindowMsecs=100; % define the size of the norm window
splitParams.normFlag=3; % percent change
% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];
IAMPFitToData = [];
% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;
temporalFit = tfeIAMP('verbosity','none');
normFlag = 3;




%grab the avgPacket which serves as the model
[ avgPackets ] = pupilPMEL_makeAverageResponsePackets( simulatedMergedPacketCellArray, normFlag, 'LowFreqClean', false, 'aggregateMethod', 'median');


% assemble the instance packet
theRunPacket=simulatedMergedPacketCellArray{1,ss}{rr};
singlePacket=splitOffAnInstancePacket(theRunPacket,splitParams);
contrast = theRunPacket.stimulus.metaData.stimTypes(ii);
% create stimulus profile -> has to be a blip with this
% configuration of IAMP (it convolves the stimulus profile
% with the kernel)
singlePacket.stimulus.values = zeros(1,13000);  % blip to be convolved with kernel; fixed per subject per contrast
singlePacket.stimulus.values(1,1) = 1;
% create kernel, which is the average run for a contrast
% level per subject
singlePacket.kernel.values = avgPackets{ss,contrast}.response.values;
singlePacket.kernel.values = singlePacket.kernel.values/(abs(min(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
if strcmp(contrast, 6);
    singlePacket.kernel.values = singlePacket.kernel.values/(abs(max(singlePacket.kernel.values))); % the average of runs of a given contrast level; fixed per subject per contrast
end
singlePacket.kernel.timebase = [0:12999];

% do the fit
[paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix);


% plot to show what's wrong
plot(singlePacket.response.values)
hold on
plot(modelResponseStruct.values)
title('Without LinearRegression Flag')
hold off

% now do the fit with linearRegression flag on
[paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix, 'searchMethod','linearRegression');
figure
plot(singlePacket.response.values)
hold on
plot(modelResponseStruct.values)
title('With LinearRegression Flag')
hold off
