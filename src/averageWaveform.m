function [avgwaves] = averageWaveform(datafiles, datapath)
%
% [avgwaves] = averageWaveform(datafiles, datapath);
%
% Computes average waveform for each unit in spike-sorted set of .nev files.
%
% N - temporal samples per spike waveform
%
% INPUTS:
%
%  datafiles - cell array of .nev files for all days and experiments
% 
%  datapath  - string specifying directory of .nev files
%
% OUTPUTS:
%
%  avgwaves - [N x # neurons] matrix where columns contain the average waveforms
%                              of each neuron
%
% @ 2016   Adam Snyder   adam@adamcsnyder.com
%          Sean Bittner  sbittner@andrew.cmu.edu
%

global FileInfo

chansToRead = 1:96;
validSortCodes = 1:254;

if ~iscell(datafiles),datafiles={datafiles};end;
datafiles = sort(datafiles);
curdir = pwd;
cd(datapath);

FileInfo = struct;
wav = cell(max(chansToRead),numel(datafiles));
uni = cell(max(chansToRead),numel(datafiles));
times = cell(max(chansToRead),numel(datafiles));
for fx = 1:numel(datafiles),
    FileInfo(fx).filename = datafiles{fx};
    D = dir(FileInfo(fx).filename);
    FileInfo(fx).bytes = D.bytes;
    nevscan(fx);
    for ch = chansToRead,
        nums = find(FileInfo(fx).PacketOrder == ch);
        if ~isempty(nums),
            loc = FileInfo(fx).HeaderSize + FileInfo(fx).Locations(nums);
            [wav{ch,fx},times{ch,fx},uni{ch,fx}] = readWaveforms(loc,FileInfo(fx).NumSamples,FileInfo(fx).filename);
            uni{ch,fx} = uni{ch,fx}';
        else
            wav{ch,fx} = int16([]); uni{ch,fx} = uint8([]);
        end;
    end;
end

%%
wav = mat2cell(wav,ones(size(wav,1),1),size(wav,2));
uni = mat2cell(uni,ones(size(uni,1),1),size(uni,2));
wav = cellfun(@cell2mat,wav,'uni',0);
uni = cellfun(@cell2mat,uni,'uni',0);
numWaves = nan(max(chansToRead),max(validSortCodes));
avgWaveform = cell(max(chansToRead),max(validSortCodes));
for ch = chansToRead,
    unitList = unique(uni{ch}(ismember(uni{ch},validSortCodes)));
    for ux = 1:numel(unitList);
        avgWaveform{ch,ux} = mean(wav{ch}(:,uni{ch}==unitList(ux)),2);               
        numWaves(ch,ux) = sum(uni{ch}==unitList(ux));
    end;
end;   

d = numel(avgWaveform);
avgwaves = [];
for j=1:d
  wave = avgWaveform{j};
  if (~isnan(wave))
    avgwaves = [avgwaves, wave];
  end
end
    
cd(curdir);

