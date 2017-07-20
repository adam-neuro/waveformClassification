function [waveforms,times,units]=readWaveforms(packetLocations,samples,filename)
% [waveforms,times,units] = readWaveforms(packetLocations,samples,filename)
%
% Code to read in waveforms from .nev file.
% 
global FileInfo
packetSize = FileInfo(1).PacketSize; % this is also in FileInfo, but I'm choosing not to pass it in %Changed from the magic number '112' to FileInfo.PacketSize. This is because we encountered a file where the assumption that the packet size was 112 bytes was incorrect, resulting in an error. -ACS 09Sep2014
% packet format:
% 4b(time) + 2b(junk) + 1b(unit) + 1b(junk) + 104b(waveform samples * 2) = 112

spikeCount = numel(packetLocations);

fid = fopen(filename,'r');

if (fid==-1)
    error(['File ',filename,' not found.']);
end

% preallocate for speed
tempData = zeros(packetSize,spikeCount,'uint8');

fseek(fid,packetLocations(1),'bof');
pl2 = diff(packetLocations);

readStr = [num2str(packetSize),'*uint8=>uint8'];
for I=1:spikeCount-1
    tempData(:,I) = fread(fid,packetSize,readStr,pl2(I)-packetSize);
end
tempData(:,spikeCount)=fread(fid,packetSize,'uint8=>uint8');

times = typecast(reshape(tempData(1:4,:),numel(tempData(1:4,:)),1),'uint32');
units = tempData(7,:)';
waveforms = typecast(reshape(tempData(9:packetSize,:),numel(tempData(9:packetSize,:)),1),'int16');
waveforms = reshape(waveforms,samples,spikeCount);

fclose(fid);
