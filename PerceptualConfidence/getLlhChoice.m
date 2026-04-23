function [choiceLlh] = getLlhChoice(stimValue, modelParams, calcPrecision, asymFlag)

% Decode function arguments
stimVal     = stimValue;                   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
noiseSens   = 1;                           % If the sensory noise is set to 1, then distributions of decision variable and confidence variable can be compared directly
guessRate   = modelParams(1);              % The fraction of guesses
stimSens    = modelParams(2);              % Stimulus sensitvity parameter, higher values produce a steeper psychometric function, strictly positive
stimCrit    = modelParams(3);              % The sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
noiseMeta   = modelParams(4);              % Meta-noise: the second stage noise parameter, only affects confidence judgments, strictly positive
confCrit    = cumsum(modelParams(5:end));  % The confidence criteria, unitless

% Set calculation precision
sampleRate = calcPrecision(1);             % Higher values produce slower, more precise estimates. Precision saturates after ~25

%% Compute model prediction
% Step 0 - rescale sensory representation by sensitivity parameter
sensMean = stimVal*stimSens;
sensCrit = stimCrit*stimSens;

for iC = 1:numel(stimVal)
    
          
    %% Compute llh of each response alternative
    % Step 1 - sample decision variable denominator in steps of constant cumulative density
    muLogN    = log((noiseSens.^2)./sqrt(noiseMeta.^2 + noiseSens.^2));
    sigmaLogN = sqrt(log((noiseMeta.^2)./(noiseSens.^2) + 1));
    dv_Den_x  = logninv(linspace(.5/sampleRate, 1-(.5/sampleRate), sampleRate), muLogN, sigmaLogN);
    
    % Step 2 - compute choice distribution under each scaled sensory distribution
    % Crucial property: linear transformation of normal variable is itself normal variable
    % Trick: we take inverse of denominator to work with products instead of ratios
    mu    = (1./dv_Den_x').*(sensMean(iC) - sensCrit);                                                     
    sigma = (1./dv_Den_x').*noiseSens;
    x     = sort([-confCrit, 0, confCrit]);
    try if asymFlag == 1; confCrit = modelParams(5:end); x = sort([-cumsum(confCrit(1:numel(confCrit)/2)), 0, cumsum(confCrit(numel(confCrit)/2 + 1:end))]); end; catch; end
    P     = normcdf(repmat(x, [sampleRate 1]), repmat(mu, [1 numel(x)]), repmat(sigma, [1 numel(x)]));
    
    % Step 3 - average across all scaled sensory distributions to get likelihood functions
    ratio_dist_p  = mean(P);   
    
    for iX = 1:numel(x)+1
        if iX == 1
            llhC{iX}(iC) = (guessRate/(numel(x)+1)) + (1 - guessRate)*ratio_dist_p(1);
        elseif (iX > 1 && iX <= numel(x))
            llhC{iX}(iC) = (guessRate/(numel(x)+1)) + (1 - guessRate)*(ratio_dist_p(iX) - ratio_dist_p(iX-1));
        elseif iX == (numel(x)+1)
            llhC{iX}(iC) = (guessRate/(numel(x)+1)) + (1 - guessRate)*(1 - ratio_dist_p(numel(x)));
        end
    end
end
choiceLlh = cell2mat(llhC');

