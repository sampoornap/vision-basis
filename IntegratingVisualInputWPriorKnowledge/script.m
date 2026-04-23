% L1: Get the posterior for a Gaussian prior and likelihood with same spread but different means
% L2: Vary the spread of the likelihood, compute the mode of the posterior
% L3: Vary the spread of the prior, compute the mode of the posterior
% L4: For all these cases, compute the spread of the posterior, compare with the likelihood spread
% L5: Simulate a maximum-likelihood cue-integration strategy for a 2AFC discrimination task
% L6: Compute the resulting within-modality and cross-modal psychometric functions


mu_prior = 2;
mu_like  = 6;
sigma = 2;

% x-axis for plotting
x = linspace(-2, 10, 1000);

% Prior and likelihood
prior = normpdf(x, mu_prior, sigma);
like  = normpdf(x, mu_like, sigma);

% posterior parameters
mu_post = (mu_prior + mu_like) / 2;
sigma_post = sigma / sqrt(2); % weighing by reciprocal of variance 

% posterior curve
post = normpdf(x, mu_post, sigma_post);


figure;
plot(x, prior, 'LineWidth', 2); hold on;
plot(x, like, 'LineWidth', 2);
plot(x, post, 'LineWidth', 2);

legend('Prior', 'Likelihood', 'Posterior');
xlabel('x');
ylabel('Density');
title('Posterior for Gaussian Prior and Likelihood w Same Spread');
grid on;


% L2 

% new prior
mu_prior = 2;
sigma_prior = 2;

% new mean
mu_like = 8;

% different likelihood spreads
sigma_like_vals = [1 2 3 5 8];

% for storing posterior modes
post_modes = zeros(size(sigma_like_vals));
posterior_spreads_L2 = zeros(size(sigma_like_vals));

for i = 1:length(sigma_like_vals)
    sigma_like = sigma_like_vals(i);

    % posterior variance
    sigma_post_sq = 1 / (1/sigma_prior^2 + 1/sigma_like^2);
    posterior_spreads_L2(i) = sigma_post_sq;

    % posterior mean 
    mu_post = sigma_post_sq * (mu_prior/sigma_prior^2 + mu_like/sigma_like^2);

    post_modes(i) = mu_post;
end

% res
disp(table(sigma_like_vals', post_modes', 'VariableNames', {'likelihoodSpread', 'posteriorMode'}));


figure;
plot(sigma_like_vals, post_modes, '-o', 'LineWidth', 2);
xlabel('likelihood spread');
ylabel('posterior mode');
title('Effects of likelihood spread on posterior mode');
grid on;

% res L4

disp(table(sigma_like_vals', posterior_spreads_L2', 'VariableNames', {'likelihoodSpread', 'posteriorSpread'}));


figure;
plot(sigma_like_vals, posterior_spreads_L2, '-o', 'LineWidth', 2);
xlabel('likelihood spread');
ylabel('posterior spread');
title('Effects of likelihood spread on posterior spread L2');
grid on;



% L3 

% new likelihood
mu_like = 8;
sigma_like = 2;

% new prior
mu_prior = 2;
sigma_prior_vals = [1 2 3 5 8];

% for storing posterior modes
post_modes = zeros(size(sigma_prior_vals));
posterior_spreads_L3 = zeros(size(sigma_prior_vals));
likelihood_spreads = zeros(size(sigma_prior_vals));

for i = 1:length(sigma_prior_vals)
    sigma_prior = sigma_prior_vals(i);

    % posterior variance
    sigma_post_sq = 1 / (1/sigma_prior^2 + 1/sigma_like^2);
    posterior_spreads_L3(i) = sigma_post_sq;
    likelihood_spreads(i) = sigma_like;

    % posterior mean 
    mu_post = sigma_post_sq * (mu_prior/sigma_prior^2 + mu_like/sigma_like^2);

    post_modes(i) = mu_post;
end

% res
disp(table(sigma_prior_vals', post_modes', 'VariableNames', {'priorSpread', 'posteriorMode'}));


figure;
plot(sigma_like_vals, post_modes, '-o', 'LineWidth', 2);
xlabel('prior spread');
ylabel('posterior mode');
title('Effects of prior spread on posterior mode');
grid on;

% res L4

disp(table(likelihood_spreads', posterior_spreads_L3', 'VariableNames', {'likelihoodSpread', 'posteriorSpread'}));


figure;
plot(likelihood_spreads, posterior_spreads_L3, '-o', 'LineWidth', 2);
xlabel('likelihood spread');
ylabel('posterior spread');
title('Effects of likelihood spread on posterior spread L3');
grid on;


% L4 - For all these cases, compute the spread of the posterior, compare with the likelihood spread

% done 

% L5: Simulate a maximum-likelihood cue-integration strategy for a 2AFC discrimination task

standard = 0;
comparisons = -5:1:5;

% noise
sigma1 = 2;   
sigma2 = 1;   % less noisy (like vision)

% MLE weights
w1 = (1/sigma1^2) / (1/sigma1^2 + 1/sigma2^2);
w2 = (1/sigma2^2) / (1/sigma1^2 + 1/sigma2^2);

nTrials = 1000;
p_choose_comparison = zeros(size(comparisons)); % to store frac of times comparison is taller
for i = 1:length(comparisons) % 11
    comp = comparisons(i);
    chooseComp = 0;

    for t = 1:nTrials
        % noisy cue measurements for standard
        x1_s = standard + sigma1*randn;
        x2_s = standard + sigma2*randn;

        % noisy cue measurements for comparison
        x1_c = comp + sigma1*randn;
        x2_c = comp + sigma2*randn;

        % MLE combined estimates
        s_hat = w1*x1_s + w2*x2_s;
        c_hat = w1*x1_c + w2*x2_c;

        % 2AFC final decision
        if c_hat > s_hat
            chooseComp = chooseComp + 1;
        end
    end

    p_choose_comparison(i) = chooseComp / nTrials;
end

figure;
plot(comparisons, p_choose_comparison, '-o', 'LineWidth', 2);
xlabel('comparison stimulus height');
ylabel('Pr(comparison is taller)');
title('2AFC with MLE cues integration');
grid on;

% L6: Compute the resulting within-modality and cross-modal psychometric functions

standard = 0;
comparisons = -5:1:5;

% noise
sigma1 = 2;   
sigma2 = 1;   % less noisy (like vision)

% MLE weights
w1 = (1/sigma1^2) / (1/sigma1^2 + 1/sigma2^2);
w2 = (1/sigma2^2) / (1/sigma1^2 + 1/sigma2^2);

nTrials = 1000;
p_choose_comparison = zeros(size(comparisons)); % to store frac of times comparison is taller
for i = 1:length(comparisons) % 11
    comp = comparisons(i);
    chooseComp = 0;

    for t = 1:nTrials
        % noisy cue measurements for standard
        x1_s = standard + sigma1*randn;
        % x2_s = standard + sigma2*randn;

        % noisy cue measurements for comparison
        x1_c = comp + sigma1*randn;
        % x2_c = comp + sigma2*randn;

        % MLE combined estimates
        s_hat = x1_s;
        c_hat = x1_c;

        % 2AFC final decision
        if c_hat > s_hat
            chooseComp = chooseComp + 1;
        end
    end

    p_choose_comparison(i) = chooseComp / nTrials;
end

figure;
plot(comparisons, p_choose_comparison, '-o', 'LineWidth', 2);
xlabel('comparison stimulus height');
ylabel('Pr(comparison is taller)');
title('2AFC with within modality MLE cues integration');
grid on;