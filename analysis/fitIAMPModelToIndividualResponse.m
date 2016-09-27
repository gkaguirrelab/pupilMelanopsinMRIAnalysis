function [theResult] = fitIAMPModelToIndividualResponse(accumStimTypesResponse,accumStimTypesStimulus, avgPackets)

% Loop across stimulus types (contrast level)

% Build a single packet container
  % the appropriate stim type
  % the kernel (which is the average response across trials for that stim
  % type, available in average packets)
  
    % Loop over stimulus instances
      % place the response into the packet
      % perform the model fit
      % record the amplitude
      % record the intial value of the response

% Loop structure
%  sessions
%    stimTypes (contrast levels plus attention task)
%       [make a packet that has a given stimulus profile and kernel]
%        instances / events
      




% Define a parameter lock matrix, which in this case is empty
paramLockMatrix = [];

% We will fit each average response as a single stimulus in a packet, so
% each packet therefore contains a single stimulus instance.
defaultParamsInfo.nInstances = 1;

        % Construct the model object
        temporalFit = tfeIAMP('verbosity','high');

% Announce what we are about to do
fprintf('>> Fitting individual amplitude model to data (IAMP)\n');

% Loop here
for ss = 1:NSessionsMerged
    for mm = 1:NStimTypes-1
        % Update the user
        fprintf('* Subject, stimulus <strong>%g</strong> , <strong>%g</strong>', ss, mm);
        fprintf('\n');
        
        % Update the packet to have current kernel and stimulus profile
        % Note that the kernel should be scaled to have a unit excursion
        singlePacket=avgPackets{ss, mm};

        % Loop over instances / events
        for ii = 1:NEvents

            % Update the packet to have the response values for this event
            singlePacket.response.values=WHATEVER;
            
         % Conduct the fit
           [paramsFit,fVal,modelResponseStruct] = temporalFit.fitResponse(singlePacket, 'defaultParamsInfo', defaultParamsInfo,'paramLockMatrix',paramLockMatrix);
        
        % Store the fitResponse
        % Save the paramsFit. You'll discover that this is a structure,
        % with one field being "amplitude". Save the value associated with
        % that.
        twoComponentFitToData{ss,mm}.paramsFit=paramsFit;
        twoComponentFitToData{ss,mm}.fVal=fVal;
        twoComponentFitToData{ss,mm}.modelResponseStruct=modelResponseStruct;
        
        end % loop over events
        
    end % loop over stimulus types
end % loop over sessions

        % Clear the object
        delete(temporalFit);

end % function