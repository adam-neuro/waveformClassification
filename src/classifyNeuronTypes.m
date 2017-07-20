function [Pall,MoG] = classifyNeuronTypes(coefVals, okInds, MoG_prior)
%
% [Pall,MoG] = classifyNeuronTypes(coefVals, okInds)
%
% Fits a mixture of two gaussians to the average waveform parameters of neurons
% with satisfactory signal-to-noise ratios and shape correlations.  The average
% waveform parameters used are time of AHP (coefVals(:,5)), 
% width of AHP (coefVals(:,6)), and rise + fall rate (coefVals(:,7)+coefVals(:,10).
%
% INPUTS:
% 
% coefVals - [# neurons x # parameters (12)] matrix where rows are the fitted
%            waveform parameters of each neuron.
% okInds - [# neurons x 1] boolean vector indicating neurons with satisfactor 
%          signal-to-noise ratio and shape correlation.
%
% OUTPUTS:
%
% Pall - [# neurons, 2] posterior probability of each neuron belonging to a
%        particular class (1st column inhibitory (fast spiking, 2nd column 
%        excitatory (regual spiking)).
% MoG - contains properties of mixture of gaussians fit to satisfactory neurons
%
% @ 2016 Adam Snyder    adam@adamcsnyder.com
%        Sean Bittner   sbittner@andrew.cmu.edu
%
  X = [coefVals(:,5),coefVals(:,6),coefVals(:,10)+coefVals(:,7)];
  if (isempty(MoG_prior))
    gmdistopt = statset('MaxIter', 100);
    MoG = gmdistribution.fit(X((okInds==1),:),2,'start','randsample','replicates',10); 
  else
    MoG = MoG_prior;
  end
  Pall = posterior(MoG,X); 
  mu = MoG.mu;
   % shorter time of AHP indicates the fast class 
  [~,fastComp] = min(mu(:,1));
  [~,slowComp] = max(mu(:,1));
  %use fastComp above to sort columns of P for output. -acs... e.g.:
  Pall = Pall(:,[fastComp,slowComp]);
end  

