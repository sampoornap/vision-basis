%% Luminance / Contrast statistics for Natural vs Phase-scrambled images
% L1: local luminance distribution (natural)
% L2: local contrast distribution (natural)
% L3: bivariate luminance-contrast distribution (natural)
% L4: local luminance distribution (phase scrambled)
% L5: local contrast distribution (phase scrambled)
% L6: bivariate luminance-contrast distribution (phase scrambled)


thisPath   = pwd;
natImPath  = fullfile(thisPath, 'NatImages');
windowPath = fullfile(thisPath, 'Window');

cd(windowPath);

[winIdx, map] = imread('Flattop8.tif');

if ~isempty(map)
    window = ind2gray(winIdx, map);
else
    window = im2double(winIdx);
end

window = im2double(window);
window = window(:,:,1);          

% normalize window weights
w = window;
w = w / sum(w(:));

winsize = size(w);               
winsize = winsize(1:2);          
lowR  = floor((winsize(1)-1)/2);
highR = winsize(1) - lowR - 1;

lowC  = floor((winsize(2)-1)/2);
highC = winsize(2) - lowC - 1;

cd(natImPath);

imNames = { ...
    'bridge.jpeg','city.jpeg','daisies.jpeg','face.jpeg', ...
    'highrise.jpeg','lion.jpeg','apple.jpeg'};

images = cell(numel(imNames),1);
for k = 1:numel(imNames)
    I = im2double(imread(imNames{k}));
    if size(I,3) > 1, I = rgb2gray(I); end % using grayscale luminance
    images{k} = I;
end

cd(thisPath);


numPatchesPerImage = 5000;  
useRelativeContrast = false; 

% storage (natural)
L_nat = zeros(numel(images)*numPatchesPerImage,1);
C_nat = zeros(numel(images)*numPatchesPerImage,1);

% storage (phase scrambled)
L_ps  = zeros(numel(images)*numPatchesPerImage,1);
C_ps  = zeros(numel(images)*numPatchesPerImage,1);

%% (keeps amplitude spectrum, randomizes phase)
phaseScramble = @(I) localPhaseScramble(I);

%% sample patches and compute local luminance/contrast
idx = 1;
for k = 1:numel(images)
    I = images{k};
    [imrows, imcols] = size(I);

    % phase scrambled version 
    IPS = phaseScramble(I);

    % (avoid borders)
    rmin = 1 + lowR;
    rmax = imrows - highR;
    cmin = 1 + lowC;
    cmax = imcols - highC;

    for n = 1:numPatchesPerImage
        % picking a random top-left corner so the patch is exactly the window size
        r0 = randi([1, imrows - winsize(1) + 1]);
        c0 = randi([1, imcols - winsize(2) + 1]);
        
        patchNat = I(  r0:r0+winsize(1)-1,  c0:c0+winsize(2)-1 );
        patchPS  = IPS(r0:r0+winsize(1)-1,  c0:c0+winsize(2)-1 );

        % windowed / weighted stats
        [Ln, Cn] = weightedLumContrast(patchNat, w, useRelativeContrast);
        [Lp, Cp] = weightedLumContrast(patchPS,  w, useRelativeContrast);

        L_nat(idx) = Ln;  C_nat(idx) = Cn;
        L_ps(idx)  = Lp;  C_ps(idx)  = Cp;

        idx = idx + 1;
    end
end

%%L1/L2/L3
figure('Name','Natural image patch statistics'); clf;

subplot(2,3,1);
histogram(L_nat, 60, 'Normalization','pdf');
xlabel('Local luminance (weighted mean)'); ylabel('PDF');
title('L1: Natural luminance');

subplot(2,3,2);
histogram(C_nat, 60, 'Normalization','pdf');
xlabel('Local contrast'); ylabel('PDF');
title('L2: Natural contrast');

subplot(2,3,3);
plotBivariate(L_nat, C_nat);
title('L3: Natural bivariate (L,C)');

%%L4/L5/L6
subplot(2,3,4);
histogram(L_ps, 60, 'Normalization','pdf');
xlabel('Local luminance (weighted mean)'); ylabel('PDF');
title('L4: Phase-scrambled luminance');

subplot(2,3,5);
histogram(C_ps, 60, 'Normalization','pdf');
xlabel('Local contrast'); ylabel('PDF');
title('L5: Phase-scrambled contrast');

subplot(2,3,6);
plotBivariate(L_ps, C_ps);
title('L6: Phase-scrambled bivariate (L,C)');


outPath = fullfile(thisPath, 'patch_stats.mat');
save(outPath, 'L_nat','C_nat','L_ps','C_ps','winsize','numPatchesPerImage','useRelativeContrast');
disp(['Saved stats to: ' outPath]);


function [L, C] = weightedLumContrast(patch, w, useRelativeContrast)

    % patch to 2D grayscale 
    if ndims(patch) == 3
        if size(patch,3) >= 3
            patch = 0.2989*patch(:,:,1) + 0.5870*patch(:,:,2) + 0.1140*patch(:,:,3);
        else
            patch = mean(patch,3);
        end
    end

    % w must match patch
    assert(isequal(size(patch), size(w)), 'Patch and window size mismatch.');

    % weighted mean luminance 
    L = sum(patch(:) .* w(:));

    % weighted variance E[x^2] - (E[x])^2
    Ex2  = sum((patch(:).^2) .* w(:));
    varW = max(Ex2 - L.^2, 0);
    stdW = sqrt(varW);

    if useRelativeContrast
        C = stdW / max(L, 1e-6);
    else
        C = stdW;
    end
end

function Iout = localPhaseScramble(Iin)
    % phase-scramble a real image while preserving amplitude spectrum
    
    I = double(Iin);

    F = fft2(I);
    A = abs(F);

    % (keep DC phase = original to preserve mean better)
    phi = angle(F);
    randPhi = (2*pi) * rand(size(I)) - pi;
    randPhi(1,1) = phi(1,1);

    Fscr = A .* exp(1i * randPhi);
    Iscr = real(ifft2(Fscr));

    % mormalizing to [0,1] (so patch stats comparable scale-wise)
    Iscr = Iscr - min(Iscr(:));
    Iscr = Iscr / max(Iscr(:) + eps);

    Iout = Iscr;
end

function plotBivariate(L, C)
    nb = 60;
    edgesL = linspace(min(L), max(L), nb);
    edgesC = linspace(min(C), max(C), nb);
    N = histcounts2(L, C, edgesL, edgesC, 'Normalization','pdf');

    imagesc(edgesL, edgesC, N'); axis xy;
    xlabel('Local luminance'); ylabel('Local contrast');
    colorbar;
end