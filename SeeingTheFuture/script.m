% L1: Create a 1-D temporal trajectory by summing a slowly changing and a random component
% L2: Create a variable that allows to control the SNR of this trajectory
% L3: Measure and plot the autocorrelation of the trajectory under different SNR values
% L4: Measure the local curvature of the video generate by the function makeVideo
% L5: Process the video with a pair of 2-D filters, compute the curvature of the filter output
% L6: Specify a 2-D 11-point perceptual trajectory, use this to simulate responses in an AXB task


% L1 + L2
nFrames = 11;
stepSize = 15;
imgFile = 'apple.jpeg';

t = 1:nFrames;

slowComp = sin(linspace(0, pi, nFrames));   % slowly varying component
noiseRaw  = randn(1, nFrames);              % random component

snrVals = [0.2 1 5];   % low, medium, high SNR

figure;
for k = 1:numel(snrVals)
    SNR = snrVals(k);

    % signal power / noise power = SNR
    noiseComp = noiseRaw / std(noiseRaw);
    noiseComp = noiseComp * std(slowComp) / sqrt(SNR);

    traj = slowComp + noiseComp;

    subplot(numel(snrVals),1,k);
    plot(t, traj, '-o', 'LineWidth', 2);
    xlabel('Frame');
    ylabel('Value');
    title(sprintf('L1/L2: trajectory, SNR = %.2f', SNR));
    grid on;
end

%% L3: autocorrelation for different SNRs

figure; hold on;
legendStr = cell(numel(snrVals),1);

for k = 1:numel(snrVals)
    SNR = snrVals(k);

    noiseComp = noiseRaw / std(noiseRaw);
    noiseComp = noiseComp * std(slowComp) / sqrt(SNR);

    traj = slowComp + noiseComp;

    [acf, lags] = xcorr(traj - mean(traj), 'coeff');
    plot(lags, acf, 'LineWidth', 2);
    legendStr{k} = sprintf('SNR = %.2f', SNR);
end

xlabel('Lag');
ylabel('Autocorrelation');
title('L3: autocorrelation under different SNRs');
legend(legendStr, 'Location', 'best');
grid on;


%% L4: local curvature of video from makeVideo

makeVideo
video = frameMat;

X = zeros(nFrames, numel(video(:,:,1)));
for i = 1:nFrames
    frame = video(:,:,i);   % frame
    X(i,:) = frame(:);      % vectorize
end

curvPix = localCurvature(X);

figure;
plot(2:nFrames-1, curvPix, '-o', 'LineWidth', 2);
xlabel('Frame index');
ylabel('Local curvature (deg)');
title('L4: local curvature in pixel space');
grid on;

%% show frames
figure;
for i = 1:nFrames
    subplot(3,4,i);
    imagesc(video(:,:,i)); axis image off; colormap gray;
    title(sprintf('F%d', i));
end


%% L5: process with a pair of 2-D filters, compute curvature of filter output

sz = 31;
[x,y] = meshgrid(linspace(-1,1,sz), linspace(-1,1,sz));
sigma = 0.3;
freq = 4;
theta = pi/4;

xr =  x*cos(theta) + y*sin(theta);
yr = -x*sin(theta) + y*cos(theta);

gEnv = exp(-(xr.^2 + yr.^2)/(2*sigma^2));
filt1 = gEnv .* cos(2*pi*freq*xr);   % even
filt2 = gEnv .* sin(2*pi*freq*xr);   % odd

% normalize
filt1 = filt1 / norm(filt1(:));
filt2 = filt2 / norm(filt2(:));

% represent each frame by filter outputs
Y = zeros(nFrames, 2);
for i = 1:nFrames
    frame = video(:,:,i);
    r1 = sum(sum(conv2(frame, filt1, 'valid')));
    r2 = sum(sum(conv2(frame, filt2, 'valid')));
    Y(i,:) = [r1 r2];
end

curvFilt = localCurvature(Y);

figure;
subplot(1,2,1);
plot(Y(:,1), Y(:,2), '-o', 'LineWidth', 2);
xlabel('Filter 1 output');
ylabel('Filter 2 output');
title('L5: 2-D filter trajectory');
grid on;

subplot(1,2,2);
plot(2:nFrames-1, curvFilt, '-o', 'LineWidth', 2);
xlabel('Frame index');
ylabel('Local curvature (deg)');
title('L5: curvature of filter output');
grid on;

%% L6: specify a 2-D 11-point perceptual trajectory and simulate AXB
% slightly curved trajectory
traj2D = [linspace(0,10,nFrames)', 0.8*sin(linspace(0,pi,nFrames))'];

figure;
plot(traj2D(:,1), traj2D(:,2), '-o', 'LineWidth', 2);
xlabel('Dim 1');
ylabel('Dim 2');
title('L6: specified 2-D perceptual trajectory');
grid on;

% pairwise distances
D = squareform(pdist(traj2D));

% simulating AXB responses
nTrialsPerPair = 50;
sigmaNoise = 1.0;   % perceptual noise
pCorrect = zeros(nFrames);

for i = 1:nFrames
    for j = 1:nFrames
        if i == j
            continue;
        end

        d = D(i,j);

        % simple psychometric mapping from distance to AXB accuracy
        % larger distance -> higher correct probability
        p = 0.5 + 0.5 * erf(d / (2*sigmaNoise));

        resp = rand(nTrialsPerPair,1) < p;
        pCorrect(i,j) = mean(resp);
    end
end

figure;
imagesc(pCorrect, [0.5 1]); colorbar;
axis image;
xlabel('Frame j');
ylabel('Frame i');
title('L6: simulated AXB proportion correct');

function c = localCurvature(X)
n = size(X,1);
c = zeros(n-2,1);

for t = 2:n-1
    v1 = X(t,:)   - X(t-1,:);
    v2 = X(t+1,:) - X(t,:);

    if norm(v1) < 1e-12 || norm(v2) < 1e-12
        c(t-1) = 0;
    else
        cosang = dot(v1,v2) / (norm(v1)*norm(v2));
        cosang = max(-1, min(1, cosang));
        c(t-1) = acosd(cosang);
    end
end
end