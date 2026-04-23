% L1: Use getLlhChoice to compute the choice pdfs for a 4-AFC confidence task
% L2: Compute the associated psychometric function and confidence function
% L3: Compute the psychometric function conditioned on the confidence report
% L4: Vary the level of meta-uncertainty, evaluate the effect on these plots
% L5: Use the model as generative model to simulate 500 trials in the task
% L6: Fit the model to the synthetic data, evaluate the parameter recovery

% L1 
% stimulus conditions
stimValue = linspace(-3, 3, 101); % moves from strong evidence for A to strong evi for B
calcPrecision = [200]; 

% to use symmetric confidence criteria which means confidence criteria are symmetric
asymFlag = 0; 

% model params
guessRate = 0.02;   
stimSens  = 1.5;    
stimCrit  = 0;      
noise = 0.25;    
conf  = 1.0;    
modelParams = [guessRate, stimSens, stimCrit, noise, conf];


choiceLlh = getLlhChoice(stimValue, modelParams, calcPrecision, asymFlag); % the response likelihoods

figure;
plot(stimValue, choiceLlh(1,:), 'LineWidth', 2); hold on;
plot(stimValue, choiceLlh(2,:), 'LineWidth', 2);
plot(stimValue, choiceLlh(3,:), 'LineWidth', 2);
plot(stimValue, choiceLlh(4,:), 'LineWidth', 2);

xlabel('stimulus value');
ylabel('probability');
title('choice PDFs');

grid on;


% L2
psychometric = choiceLlh(3,:) + choiceLlh(4,:);
confFunction = choiceLlh(1,:) + choiceLlh(4,:);


% L3

% total probability of low-confidence responses
pLow = choiceLlh(2,:) + choiceLlh(3,:);

% total probability of high-confidence responses
pHigh = choiceLlh(1,:) + choiceLlh(4,:);

% conditional psychometric functions
psychLow  = choiceLlh(3,:) ./ max(pLow, eps);
psychHigh = choiceLlh(4,:) ./ max(pHigh, eps);

% Plot
figure;
plot(stimValue, psychLow, 'LineWidth', 2); hold on;
plot(stimValue, psychHigh, 'LineWidth', 2);
xlabel('Stimulus value');
ylabel('P(choose B | confidence)');
title('psychometric Function conditioned on Confidence');
legend('Low confidence', 'High confidence', 'Location', 'best');
grid on;


% L5: 500 trials 

nTrials = 500;

% sample stimulus values for each trial
stimIdx = randi(numel(stimValue), [nTrials, 1]); % (randomly pick 500 stimulus conditions)
stimTrial = stimValue(stimIdx)';

% storing responses
respTrial = zeros(nTrials, 1);

for t = 1:nTrials
    p = choiceLlh(:, stimIdx(t));   % probabilities for this stimulus
    
    % cumulative distribution
    cdf = cumsum(p);
    
    u = rand;
    
    % sample response (1 to 4)
    respTrial(t) = find(u <= cdf, 1, 'first');
end



% L6

function nll = negLogLikelihood(params, stimTrial, respTrial, calcPrecision, asymFlag)
    
    % getting predicted probabilities
    llh = getLlhChoice(stimTrial, params, calcPrecision, asymFlag);
    
    % then extract probability of actual responses
    probs = zeros(length(respTrial),1);
    
    for t = 1:length(respTrial)
        probs(t) = llh(respTrial(t), t);
    end
    
    probs = max(probs, eps);
    
    nll = -sum(log(probs));
end
    


% initial guess (can be slightly off from true params)
initParams = [0.05, 1.0, 0.1, 0.3, 0.8];

% bounds
lb = [0, 0, -inf, 0, 0];
ub = [1, inf, inf, inf, inf];

% optimization
options = optimset('Display','iter');

fitParams = fmincon(@(p) negLogLikelihood(p, stimTrial, respTrial, calcPrecision, asymFlag), ...
                   initParams, [], [], [], [], lb, ub, [], options);


disp('True params:');
disp(modelParams);

disp('Recovered params:');
disp(fitParams);