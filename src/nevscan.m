function nevscan(j)

% Reads header and packet location information from .nev file.
global FileInfo;
fid=fopen([FileInfo(j).filename],'r','l');

%read general header
FileType=fread(fid,8,'char');
Version=fread(fid,2,'uchar');
FileFormatAdditional=fread(fid,2,'char');
HeaderSize=fread(fid,1,'uint32');
PacketSize=fread(fid,1,'uint32');
FileInfo(j).TimeResolutionTimeStamps=fread(fid,1,'uint32');
FileInfo(j).TimeResolutionSamples=fread(fid,1,'uint32');
%unsure about actualtype of TimeOrigin
TimeOrigin=fread(fid,8,'uint16');
Application=fread(fid,32,'char');
Comment=fread(fid,256,'uchar');
ExtendedHeaderNumber=fread(fid,1,'ulong');

BytesPerSample = 2;

% read extended headers
for i=1:ExtendedHeaderNumber
    Identifier=char(fread(fid,8,'char'))';
    switch Identifier
        case 'NEUEVWAV'
            ElecID=fread(fid,1,'uint16');
            PhysConnect=fread(fid,1,'uchar');
            PhysConnectPin=fread(fid,1,'uchar');
            FileInfo(j).nVperBit(ElecID)=fread(fid,1,'uint16');
            EnergyThresh=fread(fid,1,'uint16');
            FileInfo(j).HighThresh(ElecID)=fread(fid,1,'int16');
            FileInfo(j).LowThresh(ElecID)=fread(fid,1,'int16');
            SortedUnits=fread(fid,1,'uchar');
            BytesPerSample=((fread(fid,1,'uchar'))>1)+1;
            temp=fread(fid,10,'uchar');
        otherwise, % added26/7/05 after identifying bug in reading extended headers
            temp=fread(fid,24,'uchar');
    end
end

% Calculate number of packets
fseek(fid,0,1);
FileSize=ftell(fid);
PacketNumber=(FileSize-HeaderSize)/PacketSize;

% initialization
fseek(fid,HeaderSize,-1);
fread(fid,1,'uint32');


% read the data
FileInfo(j).PacketOrder=fread(fid,PacketNumber,'uint16=>uint16',PacketSize-2);%read the channel identifiers %changed from nested inside a uint16() call to just doing the casting in the call to fread -ACS 27Oct2014
fseek(fid,HeaderSize,-1);
% Times=fread(fid,PacketNumber,'uint32',PacketSize-4);%read the packet timestamps %'Times' is not used anywhere. I commented this out because it was a waste of time. -ACS 27Oct2014
fclose(fid);
FileInfo(j).Locations=[0:PacketSize:PacketSize*(PacketNumber-1)]';%The location of packets.

%Determine active channels and number of spikes on each
maxActiveChannel = max(FileInfo(j).PacketOrder(:)); %Changed from magic numbers to the active channel with the greatest index 31AUG2012 -ACS
FileInfo(j).SpikesNumber=zeros(1,maxActiveChannel);

spikesRead = 0;
for i1=1:maxActiveChannel  
    FileInfo(j).SpikesNumber(i1)=length(find(FileInfo(j).PacketOrder==i1));
    spikesRead = spikesRead + FileInfo(j).SpikesNumber(i1);
end

%Set file information global
FileInfo(j).SamplingRate=FileInfo(j).TimeResolutionTimeStamps/1000;
FileInfo(j).ActiveChannels=find(FileInfo(j).SpikesNumber);
FileInfo(j).HeaderSize=HeaderSize;
FileInfo(j).PacketSize=PacketSize;
FileInfo(j).BytesPerSample=BytesPerSample; %This can have electrode dependent values, but here only set once
FileInfo(j).NumSamples = (PacketSize - 8)/BytesPerSample;

