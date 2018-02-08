function [ rawData ] = getLabJackRawData( channelList, dataRate, duration, resolutionIndex )


% scanChannel: cell array of channel names
% dataRate: data frequency in Hz
% duration: data capture duration in seconds

% resolutionIndex = 3;


% Check: duration
% make sure duration (sec) input is an positive integer
if round(dataRate) < 1
    dataRate = 100;
    warning('dataRate must be positive. Assuming 100Hz')
else
    dataRate = round(dataRate);
end


% Check: duration
% make sure duration (sec) input is an positive integer
if round(duration) < 1
    duration = 1;
    warning('duration must be positive.')
else
    duration = round(duration);
end


numScans = duration;



% LABJACK MAGIC

%GETLABJACKDATA Summary of this function goes here
%   Detailed explanation goes here
ljmAsm = NET.addAssembly('LabJack.LJM'); %Make the LJM .NET assembly visible in MATLAB

t = ljmAsm.AssemblyHandle.GetType('LabJack.LJM+CONSTANTS');
LJM_CONSTANTS = System.Activator.CreateInstance(t); %creating an object to nested class LabJack.LJM.CONSTANTS

dispErr = true;
handle = 0;

try
    %Open first found LabJack
    [ljmError, handle] = LabJack.LJM.OpenS('ANY', 'ANY', 'ANY', handle);
%     
%     [ljmError, deviceType, connType, serialNumber, ipAddress, port, maxBytesPerMB] = LabJack.LJM.GetHandleInfo(handle, 0, 0, 0, 0, 0, 0);
%     ipAddrStr = '';
%     [ljmError, ipAddrStr] = LabJack.LJM.NumberToIP(ipAddress, ipAddrStr);
%     disp(['Opened a LabJack with Device type: ' num2str(deviceType) ', Connection type: ' num2str(connType) ','])
%     %     disp(['Serial number: ' num2str(serialNumber) ', IP address: ' char(ipAddrStr)  ', Port: ' num2str(port) ','])
%     disp(['Max bytes per MB: ' num2str(maxBytesPerMB)])
    

LabJack.LJM.eWriteName(handle, 'STREAM_RESOLUTION_INDEX', resolutionIndex);


[ljmError, value] = LabJack.LJM.eReadName(handle, 'STREAM_RESOLUTION_INDEX', 0);
fprintf('STREAM_RESOLUTION_INDEX: %i\n',value);







% Get channel addresses
    numAddresses = length(channelList);
    

    aScanListNames = NET.createArray('System.String', numAddresses); %Scan list names to stream.
    
    for n = 1:numAddresses
        aScanListNames(n) = channelList{n};
    end

    
    aScanList = NET.createArray('System.Int32', numAddresses); %Scan list addresses to stream.
    aTypes = NET.createArray('System.Int32', numAddresses); %Dummy array for aTypes parameter
    LabJack.LJM.NamesToAddresses(numAddresses, aScanListNames, aScanList, aTypes);
    scanRate = double(dataRate); %Scans per second 117647 1ch
    scansPerRead = int32(scanRate/1); % scans per read, read duraion = 1/n secs
    
    aData = NET.createArray('System.Double', numAddresses*scansPerRead);
    
    
    % pre-allocate output (values)
    rawData = nan(1,numAddresses*scanRate*duration);
    

    try
        %Configure and start stream
        [ljmError, scanRate] = LabJack.LJM.eStreamStart(handle, scansPerRead, numAddresses, aScanList, scanRate);
%         disp(['Stream started with a scan rate of ' num2str(scanRate) ' Hz.'])
    catch e
        disp(e.message)
    end
    
    
    
    for countScans = 1:numScans
        
        try
            % read (scansPerRead) from LabJack
            [ljmError, deviceScanBacklog, ljmScanBacklog] = LabJack.LJM.eStreamRead(handle, aData, 0, 0);
            
            % convert aData from .NET object to double array
            % write double array into rawData
            idxS = (countScans-1)*numAddresses*dataRate + 1;
            idxE = idxS + numAddresses*dataRate - 1;

            rawData(1,idxS:idxE) = double(aData);

            
            disp(sprintf('Scan %d/%d', countScans, numScans))
    
        catch e
            disp(e.message)
        end
    end
    
    
    % stop stream
    LabJack.LJM.eStreamStop(handle);
    
    
    
    
    
catch e
    if dispErr
        disp(e.message)
    end
end

try
    % Close handle
    LabJack.LJM.Close(handle);
catch e
    disp(e.message)
end


rawData = reshape(rawData,numAddresses,[])';



end

