function [stats] = fitModPoissModel(unitName, varargin)
%%
% fitModPoissModel fits the modulated Poisson model presented in Goris,
%   Movshon & Simoncelli (2014) to spike count data.
%   
%   Input
%   unitName is a string specifying the label of a single unit
%   fitModPoissModel(unitName, 'loadPath', path) specifies the directory
%   from which to load data files. Default is the current directory.
% 
%   Output 
%   Stats is a structure with subfields exp and fit
%   exp.stimValue --- all unique stimulus values 
%   exp.stimResp  --- mean response rates, one for every condition 
%   exp.countMean --- spike count mean
%   exp.countvar  --- spike count variance
%   exp.countCorr --- spike count correlation
%
%   fit.sigG2     --- variance of the gain
%   fit.NLL       --- negative log likelihood of the spike counts
%   fit.cpeNLL    --- cumulative probability estimate of log likelihood
%   fit.SS.stim   --- sum of squares due to stimulus
%   fit.SS.PP     --- sum of squares due to point process
%   fit.SS.gain   --- sum of squares due to gain fluctuations
%   fit.corr.PP   --- estimated point process correlation
%   fit.corr.gain --- estimated gain correlation
%   fit.corr.pred --- predicted correlation for every stimulus condition
%%

%% Relation to the paper
% In the paper, the negative binomial distribution is first introduced in equation (4):
% pN1 = (gamma(N + 1/sigG2)./(gamma(N + 1).*(gamma(1/sigG2)))) .* (((sigG2*fS*dT).^N)./((sigG2*fS*dT + 1).^(N + 1/sigG2)))
% where sigG2 is the variance of the gain, fS is firing rate expressed in
% spikes per s, dT is the stimulus duration in s, and N is spike count

% In the methods, the gamma distribution is given as:
% pG1 = ((G.^(r-1)).*exp(-G/s))./((s.^r).*gamma(r))
% where r = 1/sigG2, and s = sigG2,

% Subsequently, the negative binomial distribution is derived in the methods as:
% pN2 = (gamma(N + r)./(gamma(N + 1).*gamma(r))) .* ((1/(1+s)).^r) .* (s./(1+s)).^N
% where r = 1/sigG2, and s = sigG2*fS*dT

% For both the gamma distribution and the negative binomial distribution,
% more convenient Matlab implementations can be used:
% pN3 = nbinpdf(N, r, p)
% where mu = fS*dT, and p = r./(r + mu);
% and
% pG2 = gampdf(G, 1/sigG2, sigG2);
% where G is gain
%%



%% Begin analysis
% Get input values from varargin or assign default values
loadPath    = GetNamedInput(varargin, 'loadPath', pwd);
currentPath = pwd;

% Load data
cd(loadPath)
load(unitName)
cd(currentPath)

% Set constants
nUnits = size(S.spikeCounts, 1);

% Get stimulus values and corresponding response mean, variance, and correlation
condIDs = unique(S.condVec(isfinite(S.condVec)));
for iC = 1:numel(condIDs)
    condRepeats(iC)       = sum(S.condVec == condIDs(iC));
    condMeanRates(:,iC)   = mean(S.spikeCounts(:, S.condVec == condIDs(iC))./repmat(S.trialDur(S.condVec == condIDs(iC))'/1000, [nUnits 1]), 2);
    condMeanCounts(:,iC)  = mean(S.spikeCounts(:, S.condVec == condIDs(iC)), 2);
    condVarCounts(:,iC)   = var(S.spikeCounts(:, S.condVec == condIDs(iC)), [], 2);
    condCorCounts(:,:,iC) = corr(S.spikeCounts(:, S.condVec == condIDs(iC))');
end
%%


%% I. Model-based analysis of response variance
% Set fit options
options = optimset('Display', 'off', 'Maxiter', 10^5, 'MaxFuneval', 10^5);

% Select appropriate trials -- leave out NaN's
passCondID    = S.condVec(isfinite(S.condVec));
passDurations = S.trialDur(isfinite(S.condVec));

for iU = 1:nUnits
    fprintf('Fitting responses of unit %d of %d... \n', iU, nUnits);
    
    % Fit the modulated Poisson model
    passSpikeCounts      = S.spikeCounts(iU, isfinite(S.condVec))';
    obFun                = @(param) giveNLL(param, condIDs, condMeanRates(iU,:), passCondID, passDurations, passSpikeCounts);
    [sigG2(iU), NLL(iU)] = fminsearch(obFun, rand, options);
    
    % Analysis of variance under modulated Poisson model -- summarize predicted variance due to stimulus, point process, and gain fluctuations
    SSstim(iU) = sum(condRepeats.*(condMeanCounts(iU,:) - mean(S.spikeCounts(iU, isfinite(S.condVec)))).^2);
    SSpp(iU)   = sum(condRepeats.*condMeanCounts(iU,:));
    SSgain(iU) = sum(condRepeats.*(sigG2(iU)*(condMeanCounts(iU,:).^2)));
    
    % Evaluate absolute goodness of fit through a 1000-fold parametric bootstrap
    [~, loc] = ismember(passCondID, condIDs);
    mu       = max(.001, .001*passDurations.*condMeanRates(iU, loc)');
    r        = 1/sigG2(iU);
    p        = r./(r + mu);
    
    LLBootstrap = NaN + zeros(1000, 1);
    for iR = 1:1000
        bootCounts      = nbinrnd(r, p);
        llhBootCounts   = nbinpdf(bootCounts, r, p);
        LLBootstrap(iR) = sum(log(llhBootCounts));
    end
    
    % Express goodness of fit with a cumulative probability estimate (0 = all simulated data sets are more probable than the observed data; 1 = all simulated data sets are less probable than the observed data)
    cpeLL(iU) = mean(LLBootstrap <= -NLL(iU));
end
%%



%% II. Model-based analysis of response correlation
% Loop throught all pairs, infer gain- and point-process correlation
corrPP        = nan + zeros(nUnits, nUnits);
corrGain      = nan + zeros(nUnits, nUnits);
corrCountPred = nan + zeros(nUnits, nUnits, numel(condIDs));
for iU = 1:nUnits
    for jU = 1:nUnits
        
        if jU < iU
            
            fprintf('Fitting response correlations of unit %d and %d of %d... \n', iU, jU, nUnits);
            
            % Fit the modulated Poisson model
            mu_i       = condMeanCounts(iU,:);
            mu_j       = condMeanCounts(jU,:);
            sigG2_i    = sigG2(iU);
            sigG2_j    = sigG2(jU);
            corrObs_ij = squeeze(condCorCounts(iU,jU,:))';
            obFun      = @(params) giveSSE(params, mu_i, mu_j, sigG2_i, sigG2_j, corrObs_ij);
            paramsEst  = fminsearch(obFun, rand(1,2), options);
            
            [~, corrPred_ij]       = giveSSE(paramsEst, mu_i, mu_j, sigG2_i, sigG2_j, corrObs_ij);
            corrPP(iU,jU)          = (exp(2*paramsEst(1))-1)./((exp(2*paramsEst(1))+1));
            corrGain(iU,jU)        = (exp(2*paramsEst(2))-1)./((exp(2*paramsEst(2))+1));
            corrCountPred(iU,jU,:) = corrPred_ij;            
        end
    end
end
condCorCounts(~isfinite(corrCountPred)) = NaN;
%%



%% Summarize useful statistics
stats.exp.stimValue = condIDs;
stats.exp.stimResp  = condMeanRates;
stats.exp.countMean = condMeanCounts;
stats.exp.countVar  = condVarCounts;
stats.exp.countCorr = condCorCounts;
stats.fit.sigG2     = sigG2;
stats.fit.NLL       = NLL;
stats.fit.cpeLL     = cpeLL;
stats.fit.SS.stim   = SSstim;
stats.fit.SS.pp     = SSpp;
stats.fit.SS.gain   = SSgain;
stats.fit.corr.PP   = corrPP;
stats.fit.corr.gain = corrGain;
stats.fit.corr.pred = corrCountPred;

%%
end






%%
function [NLL] = giveNLL(param, condIDs, condMeanRates, passCondID, passDurations, passSpikeCounts)
% the negative binomial distribution can be parameterized with its mean and
% variance. The maximum likelihood estimator for the mean is the observed
% sample mean. To maximize the likelihood of the observed spike counts
% under the model, we only need to find the optimal variance of the gain.

% Set the variance of the gain
sigG2 = param;

% Get the predicted spike count distribution for every trial
[~, loc] = ismember(passCondID, condIDs);
mu       = max(.001, .001*passDurations.*condMeanRates(loc)');
r        = 1/sigG2;
p        = r./(r + mu);

% Evaluate the likelihood of each spike count
llh = nbinpdf(passSpikeCounts, r, p);

% Get the negative log-likelihood of the whole data-set
NLL = sum(-log(llh));
end
%%



%%
function [SSE, corrPred_ij] = giveSSE(params, mu_i, mu_j, sigG2_i, sigG2_j, corrObs_ij)
% there is no multivariate extension of the negative binomial distribution,
% but we can decompose the spike count covariance into a point process and
% a gain term. To estimate both components, we evaluate how spike count
% correlation depends on the pairwise response strength.

% Set the point process and gain correlation
r_Pij = (exp(2*params(1))-1)./((exp(2*params(1))+1));
r_Gij = (exp(2*params(2))-1)./((exp(2*params(2))+1));

% Get the predicted spike count variance for every stimulus condition
varPred_i  = mu_i + sigG2_i*mu_i.^2;
varPred_j  = mu_j + sigG2_j*mu_j.^2;

% Get the predicted spike count covariance for every stimulus condition
covPred_ij = r_Pij*sqrt(mu_i.*mu_j) + r_Gij*sqrt(sigG2_i)*sqrt(sigG2_j)*(mu_i.*mu_j);

% Get the predicted spike count correlation for every stimulus condition
corrPred_ij = covPred_ij./sqrt(varPred_i.*varPred_j);

% Fisher transform the correlations
fishObs_ij  = .5*log((1 + corrObs_ij)./(1 - corrObs_ij));
fishPred_ij = .5*log((1 + corrPred_ij)./(1 - corrPred_ij));

% Get the SSE
SSE = sum((fishObs_ij - fishPred_ij).^2);
end
%%



%%
function y = GetNamedInput(C, varName, varDefault)
% looks for the string varName in varargin, and returns the following entry
% in varargin. If varName is named more than once, a cell array is
% returned. If it is not found, varDefault is returned.

y = varDefault;

k=0;
for i = 1:(length(C)-1)
    if strcmpi(C{i},varName)
        k = k+1; % increment k every time the varName is found in varargin
        if k>1
            y{k} = C{i+1};
        else
            y = C{i+1};
        end
    end
end
end
%%

