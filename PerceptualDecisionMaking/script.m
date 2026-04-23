% L1: Simulate a 2AFC bounded accumulation decision-making process for an ambiguous stimulus
% L2: Add several conditions with varying stimulus strength, compute the psychometric function
% L3: Compute the chronometric function
% L4: Simulate the speed-accuracy trade-off commonly seen in human decision-making
% L5: Derive the expected temporal evolution of the variance of the decision-variable, verify
% L6: Derive the expected autocorrelation function of the decision-variable, verify


%L1 - ambiguous so mean ~= 0

dt = 0.001;
T = 2; % max time
time = 0:dt:T;

mu = 0; % ambiguous
sigma = 1;
B = 1; % bound

x = 0;
traj = zeros(size(time));

for t = 1:length(time)
    dx = mu*dt + sigma*sqrt(dt)*randn;
    x = x + dx;
    traj(t) = x;

    if x >= B
        decision = 1;
        RT = time(t);
        break;
    elseif x <= -B
        decision = -1;
        RT = time(t);
        break;
    end
end

plot(time, traj); hold on;
yline(B); yline(-B);
title('Single trial accumulation');

% L2

mus = [-0.2 -0.1 0 0.1 0.2];
nTrials = 200;

accuracy = zeros(size(mus));

for i = 1:length(mus)
    correct = 0;

    for trial = 1:nTrials
        x = 0;

        for t = 1:length(time)
            dx = mus(i)*dt + sigma*sqrt(dt)*randn;
            x = x + dx;

            if x >= B
                if mus(i) > 0, correct = correct + 1; end
                break;
            elseif x <= -B
                if mus(i) < 0, correct = correct + 1; end
                break;
            end
        end
    end

    accuracy(i) = correct / nTrials;
end

plot(mus, accuracy, 'o-');
xlabel('Stimulus strength'); ylabel('P(correct)');
title('Psychometric function');


% L3

RTs = zeros(size(mus));

for i = 1:length(mus)
    rt_list = [];

    for trial = 1:nTrials
        x = 0;

        for t = 1:length(time)
            dx = mus(i)*dt + sigma*sqrt(dt)*randn;
            x = x + dx;

            if abs(x) >= B
                rt_list(end+1) = time(t);
                break;
            end
        end
    end

    RTs(i) = mean(rt_list);
end

plot(mus, RTs, 'o-');
xlabel('Stimulus strength'); ylabel('Reaction time');
title('Chronometric function');

% L4

Bs = [0.5 1 1.5];

for b = 1:length(Bs)
    B = Bs(b);

    accuracy = zeros(size(mus));

    for i = 1:length(mus)
        correct = 0;
    
        for trial = 1:nTrials
            x = 0;
    
            for t = 1:length(time)
                dx = mus(i)*dt + sigma*sqrt(dt)*randn;
                x = x + dx;
    
                if x >= B
                    if mus(i) > 0, correct = correct + 1; end
                    break;
                elseif x <= -B
                    if mus(i) < 0, correct = correct + 1; end
                    break;
                end
            end
        end
    
        accuracy(i) = correct / nTrials;
    end
    RTs = zeros(size(mus));

    for i = 1:length(mus)
        rt_list = [];
    
        for trial = 1:nTrials
            x = 0;
    
            for t = 1:length(time)
                dx = mus(i)*dt + sigma*sqrt(dt)*randn;
                x = x + dx;
    
                if abs(x) >= B
                    rt_list(end+1) = time(t);
                    break;
                end
            end
        end
    
        RTs(i) = mean(rt_list);
    end
    % store results
end


% L5
nTrials = 500;
traj_all = zeros(nTrials, length(time));

for trial = 1:nTrials
    x = 0;

    for t = 1:length(time)
        dx = sigma*sqrt(dt)*randn;
        x = x + dx;
        traj_all(trial,t) = x;
    end
end

var_t = var(traj_all);

plot(time, var_t); hold on;
plot(time, sigma^2*time, '--');
legend('Simulated','Theory');
title('Variance over time');



% L6

t1 = 200; 
t2 = 400;

x1 = traj_all(:, t1);
x2 = traj_all(:, t2);

cov_empirical = mean((x1 - mean(x1)) .* (x2 - mean(x2)));

cov_theory = sigma^2 * min(time(t1), time(t2));