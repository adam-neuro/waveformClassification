function waveforms = nevWaveforms(filenames,channels)
%function waveforms = nevWaveforms(filenames,channels)
%
% filenames can be a single string or a cell array of strings
% channels is optional, uses all channels by default
%
clear global FileInfo
global loc;
global FileInfo;
global WaveformInfo;

FileInfo=struct('filename',[],'format','nev','HeaderSize',0,...
    'PacketSize',0,'ActiveChannels',[],'PacketOrder',uint8([]),...
    'SpikesNumber',[],'BytesPerSample',0); %initialize FileInfo structure
        
if iscell(filenames)
    for I=1:length(filenames)
        FileInfo(I).filename = filenames{I};
    end
else
    FileInfo(1).filename = filenames;
end

ActiveChannelList=[];

for i = 1:length(FileInfo)
    nevscan(i);

    ActiveChannelList = union(ActiveChannelList,FileInfo(i).ActiveChannels);
end

if nargin < 2 % changed from 1:96 - 2013.08.23 MAS
    channels = ActiveChannelList;
end

if length(unique([FileInfo(:).PacketSize]))~=1
    error('Variable Packet Sizes');
end

WaveformSize=(FileInfo(1).PacketSize-8)/FileInfo(1).BytesPerSample;
ByteLength=['int' num2str(FileInfo(1).BytesPerSample*8)];

numSamples = (FileInfo(1).PacketSize-8)/2;

for channelIndex = 1:length(channels)
    WaveformInfo(channelIndex).Waveforms = [];
    WaveformInfo(channelIndex).Unit = [];
    WaveformInfo(channelIndex).Times = [];
    for fileIndex = 1:length(FileInfo)
      
      PacketNumbers=find(FileInfo(fileIndex).PacketOrder==channels(channelIndex));
        loc = FileInfo(fileIndex).HeaderSize + FileInfo(fileIndex).Locations(PacketNumbers);
        [wav, times, units] = readWaveforms(loc,numSamples, FileInfo(fileIndex).filename);
        WaveformInfo(channelIndex).Channel = channels(channelIndex);
        WaveformInfo(channelIndex).Waveforms = [WaveformInfo(channelIndex).Waveforms; wav'];
        WaveformInfo(channelIndex).Unit = [WaveformInfo(channelIndex).Unit; units];
        WaveformInfo(channelIndex).Times = [WaveformInfo(channelIndex).Times; double(times)/FileInfo(fileIndex).TimeResolutionTimeStamps*1000];
    end
end

waveforms = WaveformInfo;
