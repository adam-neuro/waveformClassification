function [coefVals] = fitWaveformParameters(avgwaves, start_opt_raw, lb_opt_raw, ub_opt_raw, ...
                                            start_opt_der, lb_opt_der, ub_opt_der, ...
                                            sampsPerMsec, normAvgWaveforms)
%
% [coefVals] = fitWaveformParameters(avgwaves, start_opt_raw, lb_opt_raw, ub_opt_raw, ...
%                                    start_opt_der, lb_opt_der, ub_opt_der, ...
%                                    sampsPerMsec, normAvgWaveforms)
%
% Fit waveform parameters via two gaussians to each average waveform, and 
% derivative average waveform.  
%
% INPUTS:
%
% avgwaves - [N x # neurons] matrix where columns contain the average waveforms
%                            of each neuron 
% start_opt_raw - vector of average waveform parameter fitting initialization values
%             [amplitude, mean, standard dev of leading gaussian, ...
%              amplitude, mean, standard dev of trailing gaussian].
% lb_raw - vector of lower bounds for raw average waveform parameters
% ub_raw - vector of upper bounds for raw average waveform parameters
%
% start_opt_der - vector of derivative average waveform parameter fitting initialization values
%             [amplitude, mean, standard dev of leading gaussian, ...
%              amplitude, mean, standard dev of trailing gaussian].
% lb_der - vector of lower bounds for derivative average waveform parameters
% ub_der - vector of upper bounds for derivative average waveform parameters
% 
% sampsPerMsec - integer - samples per milisecond
% normAvgWaveform - boolean indicating whether average waveforms should be
%                        set to have max. deviation of +/- 1
%
% OUTPUTS
%
% coefVals - [# neurons x # parameters (12)] matrix where rows are the fitted
%            waveform parameters of each neuron.
%
% @ 2016 Adam Snyder    adam@adamcsnyder.com
%        Sean Bittner   sbittner@andrew.cmu.edu
%


  tl = (1:size(avgwaves,1))'./sampsPerMsec; %timeline
  if (normAvgWaveforms) 
    %set maximum deviation of waveform to +/- 1.
    avgwaves = bsxfun(@rdivide,avgwaves,max(abs(avgwaves)));
  end;
  diffwaves = diff(avgwaves);
  dvTl = tl(:)<=1;
  nParams = 3*2*2; %three terms times two gaussians times two fits...
  coefVals = zeros(size(avgwaves,2),nParams); %preallocate
  
  % Fit the raw waveform
  opt_raw = fitoptions('method','nonlinearleastsquares','lower',lb_opt_raw,'upper', ub_opt_raw,'startpoint',start_opt_raw);
  % Fit the derivative waveform
  opt_der = fitoptions('method','nonlinearleastsquares','lower',lb_opt_der,'upper', ub_opt_der,'startpoint',start_opt_der); fprintf('fitting wave     ');
  for w = 1:size(avgwaves,2),
      fprintf('\b\b\b\b%4d',w);
      fObj_raw = fit(tl,avgwaves(:,w),'gauss2',opt_raw);
      fObj_der = fit(tl(dvTl),diffwaves(dvTl,w),'gauss2',opt_der);
      coefVals(w,:) = [coeffvalues(fObj_raw),coeffvalues(fObj_der)];
  end
  fprintf('\n');
end
