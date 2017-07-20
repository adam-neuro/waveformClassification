function M = justSNR(filename)
% [M] = justSNR(filename)
%
% justSNR(filename) takes an NEV file as input and returns a list of
% channel numbers, sort codes, SNR values, and spike counts from that file
%
%
  clear global WaveformInfo
  global WaveformInfo;
  
  nevWaveforms(filename);

  k = 1;
  
  waveforms = struct;
  waveInd = 1;
  for i = 1:length(WaveformInfo) 
    units = unique(WaveformInfo(i).Unit);

    for j = 1:length(units)
      wav = WaveformInfo(i).Waveforms(find(WaveformInfo(i).Unit == ...
                                           units(j)),:);
      M(k,:) = [WaveformInfo(i).Channel double(units(j)) getSNR(double(wav)) size(wav,1)];
      k = k + 1;
    end
  end
end

