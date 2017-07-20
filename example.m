function [] = example(datapath, MoG)
% example(datapath, MoG)
%
% Example of how to classify neurons using mixture of gaussians (MoG)
% model, as described in the following publication:
%
% Snyder AC, Morais MJ, Smith MA (2016). Dynamics of excitatory and
% inhibitory networks are differentially altered by selective attention.
% Journal of Neurophysiology, V(I): pp-pp. doi:11.1152/jn.00343.2016.
%
% http://jn.physiology.org/content/early/2016/07/22/jn.00343.2016.abstract
%
% INPUTS:
%
% datapath - path to .nev files. Set to current working dir when not
%            provided. Each NEV file should contain all the data from a
%            single recording session. If multiple files were collected
%            within a single session, these should be appended somehow into 
%            a single file before running this routine.
%
% MoG - mixture of gaussians distribution for the two neuron types. Mixture
%       of gaussians is fit to the files in datapath if a model is not
%       provided.
%
% -----------------------------------------
% This routine saves a MAT file with the results of intermediate steps of
% the classification procedure. The steps are as follows:
% 
% 1) Calculate the signal-to-noise ratio of the waveforms for each unit in
%    each NEV file. This step saves the channel index, sort code, and
%    calculated SNR for each unit.
%
% 2) Compute the average waveform for each unit from all the individual
%    waveform events. These are saved as 'avgwaves' in the MAT file.
%
% 3) Fit curves to the average waveforms. This mitigates discretization
%    errors in the estimates of waveform shape parameters. Saved as
%    'coefVals' in the MAT file.
%
% 4) Fit a multidimensional mixture of Gaussians to the waveform shape
%    parameters of the best-isolated (i.e., highest SNR) units across all 
%    files. The parameters of the Gaussian mixture model are saved in a
%    file named 'model.mat' as an object of the gmdistribution class named
%    'MoG'.
%
% 5) Use the fitted Gaussian mixture model to classify all the neurons in
%    sample (not just the well-isolated ones). This is expressed as a
%    posterior probability of belonging to one or the other component of
%    the Gaussian mixture. The posteriors are saved in the single-session 
%    MAT file as a variable named 'P'; the first column of P is the 
%    posterior probability that the given neuron came from the
%    fast-/narrow-spiking class, and the second column of P is the
%    posterior probability that the given neuron came from the
%    regular-/broad-spiking class.
%
% 
% This code was prepared by Sean Bittner (sbittner@andrew.cmu.edu) and Adam
% Snyder, based on the methodology developed for Snyder et al (2016).
%

fprintf(['\nExample code using mixture of multidimensional Gaussians to classify\n',...
    'neurons based on waveform shapes as described in:\n\n',...
    'Snyder AC, Morais MJ, Smith MA (2016). Dynamics of excitatory and\n',...
    'inhibitory networks are differentially altered by selective attention.\n'...
    'Journal of Neurophysiology, V(I): pp-pp. doi:11.1152/jn.00343.2016.\n\n']);

%add the directory of source code to the MATLAB search path (this assumes
%that 'src' is a subfolder of the location of this function):
addpath(fullfile(fileparts(which(mfilename)),'src'));

%Defauly behaviors for missing input arguments:
if nargin<1,
    basepath = pwd;
    datapath = [basepath, filesep, 'data', filesep];
    fprintf('Using default data path \n%s\n\n', datapath);
end
if nargin<2,
    MoG = [];
    fprintf('Fitting mixture of gaussians for classification to neurons in \n%s\n\n', datapath);
end
new_model_flg = isempty(MoG);

% Set waveform classification parameters (these are choices to be made by
% the user):
snrRatio = .25; %proportion of units to use for training Gaussian mixture model (e.g., 0.25 means use top 25% of best-isolated units)
shapeCorrelationThresh = .5; %exclude units with average waveforms that are not on average at least this correlated with the other units (catches non-neuronal waveforms)
sampsPerMsec = 30;  % 30 samps per ms (get this from the nev file)
normalizeAvgWaveforms = true; %scale waveforms to be bounded by +/- 1 (merely for visualization, doesn't affect the classification result)
validSortCodes = 1:254; %sort codes in the NEV file that correspond to putative neurons
validChannelCodes = 1:96; %channels in the NEV file to use for the analysis

% Initialization and bounds of parameters in average raw waveform fit.
% These are fairly permissive and just describe basic spike shapes (e.g.,
% negative phase precedes positive phase, etc.):
start_opt_raw = [-1 0.6 0.1, 1 1 0.3];
lb_opt_raw = [-1000 .5 0    0  .6 0];
ub_opt_raw = [0     .7 1 1000 1.4 1];

% Initialization and bounds of parameters in derivative average raw
% waveform fit.
start_opt_der = [-45 0.5 0.1, 45 0.5 0.1];
lb_opt_der = [-500 .4 0   0 .4 0];
ub_opt_der = [   0 .7 1 500 .9 1];

% Setup file I/O
assert(ischar(datapath),'First parameter "datapath" must be string.');
if exist(datapath,'dir') == 7
    if ~strcmp(datapath(end), filesep)
        datapath = [datapath, filesep];
    end
    fnames = dir([datapath, '*.nev']);
else
    error('Directory %s does not exist.\n', datapath);
end
nfiles = length(fnames);

%preallocate:
allCoefVals = [];
allOkInds = [];
numneurons = zeros(nfiles, 1);
classfnames = cell(nfiles, 1);

% Fit waveform parameters for each recording session
for i=1:nfiles
    
    %This function saves out the results of each step in a MAT file. First,
    %check if the MAT file already exists and load it if so:
    fname = fnames(i).name(1:(end-4));
    nevfname = [datapath, fnames(i).name];
    classfname = [datapath, 'class_', fname, '.mat'];
    if exist(classfname,'file') == 2
        M = load(classfname);
    else
        M = [];
    end
    
    % Compute signal-to-noise ratio for each neuron if it hasn't been done
    % already:
    if (isfield(M,'channels') ~= 1 || isfield(M,'units') ~= 1 || isfield(M,'SNR') ~= 1)
        fprintf('Computing signal-to-noise ratio for %s.\n', fname);
        X = justSNR(nevfname);
        X = X(ismember(X(:,1),validChannelCodes), :);
        X = X(ismember(X(:,2),validChannelCodes), :);
        channels = X(:,1); %#ok<*NASGU>
        units = X(:, 2);
        SNR = X(:,3);
    else
        channels = M.channels;
        units = M.units;
        SNR = M.SNR;
    end
    nneurons = length(SNR);
    fprintf('%d neurons in file %s.\n',nneurons,fname);
    numneurons(i) = nneurons;
    
    % Compute the average waveform for each neuron if it hasn't been done
    % already:
    if ~isfield(M, 'avgwaves')
        fprintf('Computing average waveforms for %s.\n', fname);
        avgwaves = averageWaveform(fnames(i).name, datapath);
    else
        avgwaves = M.avgwaves;
    end
    
    % Fit average waveform parameters. Fitting curves to the waveforms
    % mitigates quantization/discretization errors.
    if (isfield(M,'coefVals') ~= 1)
        fprintf('Computing average waveform fit parameters for %s.\n', fname);
        coefVals = fitWaveformParameters(avgwaves, start_opt_raw, lb_opt_raw, ub_opt_raw, ...
            start_opt_der, lb_opt_der, ub_opt_der, sampsPerMsec, normalizeAvgWaveforms);
    else
        coefVals = M.coefVals;
    end
    allCoefVals = [allCoefVals; coefVals]; %#ok<AGROW>
    
    % Only use neurons with satisfactory SNR and waveform shape to train
    % waveform classifier.
    sorted_SNR = sort(SNR);
    snrThresh = sorted_SNR(floor(snrRatio*nneurons));
    okSnr = SNR > snrThresh;
    r = mean(corr(avgwaves))';
    okShape = r > shapeCorrelationThresh;
    okInds = okSnr&okShape;
    allOkInds = [allOkInds; okInds]; %#ok<AGROW>
    save(classfname, 'channels', 'units', 'SNR', 'avgwaves', 'coefVals');
    classfnames{i} = classfname;
end

% Fit mixture of gaussians to neurons from all recording sessions and
% classify neurons Pall contains probabilities of belonging to fast spiking
% class in the first column, and probabilities of belonging to
% broad-spiking class in the second column.
[Pall,MoG] = classifyNeuronTypes(allCoefVals, allOkInds, MoG);
if (new_model_flg)
    modelfname = [datapath, 'model.mat'];
    save(modelfname, 'MoG');
end

% Save classified neurons by recording session
dataind = 1;
for i=1:nfiles
    nneurons = numneurons(i);
    P = Pall(dataind:(dataind+nneurons-1), :);
    fprintf('Saving %s\n', classfnames{i});
    save(classfnames{i}, 'P', 'fnames', '-append');
    dataind = dataind + nneurons;
end

% Visualize the results of the classification:
fprintf('Plotting averaged waveforms colored by classification\n');
for i=1:nfiles
    X = load(classfnames{i});
    figure;
    nsamps = size(X.avgwaves,1);
    xvals = (0:(nsamps-1))/sampsPerMsec;
    h1 = plot(xvals, bsxfun(@rdivide,X.avgwaves,max(abs(X.avgwaves))));
    box off;
    set(gca,'tickdir','out');
    cmap = colormap(jet(256));
    set(h1,{'color'},mat2cell(cmap(round(255*X.P(:,1))+1,:),ones(numel(h1),1),3));
    title(fnames(i).name, 'Interpreter', 'none');
    ylabel('normalized average waveform amplitude');
    xlabel('time (ms)');
    set(gca, 'FontSize', 12);
    h2 = colorbar;
    ylabel(h2, 'Probability of being the fast-spiking class');
    set(h2,'ytick',linspace(1,256,5),'yticklabel',linspace(0,1,5));
end

