% L1: Simulate spike trains from a homogeneous, inhomogeneous, and modulated Poisson process
% L2: Create a spike raster for each process (100 repeats of a 1-sec neuronal response)
% L3: Compute the variance-to-mean functions using randomly placed and sized counting windows
% L4: Estimate the preferred orientation and gain variability of a real V1 cell (function V1spikes)
% L5: Compute the PSTH for the cell's preferred stimulus using a 10, 25, and 50 ms boxcar filter
% L6: Compute the cell's temporal response modulation (f1/f0) for the preferred stimulus

clear; close all; clc;
rng(0);

T = 1.0;        % seconds
nTr = 50;      % trials  

function tau = HomPoiss(lambda, T)
% Homogeneous Poisson 
    tau = [];
    if lambda <= 0, return; end
    t = 0;
    while true
        t = t + (-log(rand)/lambda);
        if t > T, break; end
        tau(end+1,1) = t; 
    end
end

% Homogeneous Poisson
lambda0 = 20; % spikes per s
spk_h = cell(nTr,1);
for k=1:nTr
    spk_h{k} = HomPoiss(lambda0, T);
end



% % Modulated Poisson process
% lambda2 = @(t) 20 + 10 * sin(2 * pi * 5 * t); 
% spk_mod = cell(nTr, 1);
% for k = 1:nTr
%     spk_mod{k} = InhomPoiss(lambda2, T);
% end 



function tau = InhomPoiss(lambda, T)
% Inhomogeneous Poisson 
    tau = [];
    if lambda <= 0, return; end
    t = 0;
    while true
        t = t + (-log(rand)/lambda(t));
        if t > T, break; end
        tau(end+1,1) = t; 
    end
end




% Inhomogeneous Poisson
lambda1 = @(t) 20 + 10 * sin(2 * pi * t); 
spk_inh = cell(nTr, 1);
for k = 1:nTr
    spk_inh{k} = InhomPoiss(lambda1, T);
end 