% Demo of IAMP analysis on fMRI data

close all;

for pp = 1:9

thePacket=packetCellArray{1,pp};

% Construct the model object
temporalFit = tfeIAMP('verbosity','none');

% grab the average hrf and prepare it as a kernel
check = diff(thePacket.response.timebase);
responseDeltaT = check(1);
nSamples = ceil((thePacket.kernel.timebase(end)-thePacket.kernel.timebase(1))/responseDeltaT);
newKernelTimebase = thePacket.kernel.timebase(1):responseDeltaT:(thePacket.kernel.timebase(1)+nSamples*responseDeltaT);
hrfKernelStruct = temporalFit.resampleTimebase(thePacket.kernel,newKernelTimebase);
thePacket.kernel=prepareHRFKernel(hrfKernelStruct);

% Put the TR into the packet metaData
thePacket.response.metaData.TRmsecs=800;

% downsample the stimulus values to 100 ms deltaT to speed things up
totalResponseDuration=thePacket.response.metaData.TRmsecs * ...
    length(thePacket.response.values);
newStimulusTimebase=linspace(0,totalResponseDuration-100,totalResponseDuration/100);
thePacket.stimulus=temporalFit.resampleTimebase(thePacket.stimulus,newStimulusTimebase);

% Convert the response values to mean-centered % signal change
response=thePacket.response.values;
response=(response-mean(response))/mean(response);
thePacket.response.values=response;

% How many trial instances are in this packet?
defaultParamsInfo.nInstances = size(thePacket.stimulus.values,1);

% Define an empty locking matrix
paramLockMatrix=[];

% Perform the fit
[paramsFit,FVal,modelResponseStruct] = ...
    temporalFit.fitResponse(thePacket,...
    'defaultParamsInfo', defaultParamsInfo, ...
    'paramLockMatrix',paramLockMatrix, ...
    'searchMethod','linearRegression', ...
    'errorType','1-r2');

rSquared=1-FVal;

% plot the data and the fit

figure
plot(modelResponseStruct.timebase,modelResponseStruct.values)
hold on
plot(thePacket.response.timebase,thePacket.response.values)

% dump out the estimated percent signal change for each trial
theResult(pp,:)=paramsFit.paramMainMatrix;
theTrials(pp,:)=thePacket.stimulus.metaData.stimTypes;
end % loop over packets

% Aggregate data across trial types
for i=1:6
idx=find(theTrials==i);
meanResponseVal(i)=mean(mean(theResult(idx)));
semResponseVal(i)=std(theResult(idx))/sqrt(length(idx));
end

figure
plot(meanResponseVal(1:5)*100,'-k');
hold on
plot(meanResponseVal(1:5)*100+semResponseVal(1:5)*100,'Color',[0.5 0.5 0.5]);
plot(meanResponseVal(1:5)*100-semResponseVal(1:5)*100,'Color',[0.5 0.5 0.5]);
plot(meanResponseVal(1:5)*100,'ro');




