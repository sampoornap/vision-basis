% Start with a clean slate
clc
clear all
close all

% L1 – L3
% Specify the wavelength spectrum sampling
wavelengthVec = linspace(350, 750, 401);

% Specify the amplitude of the stimulus reflected light
% white light
whiteLightAmp = ones(size(wavelengthVec));

% 500 nm light
nm500LightAmp           = zeros(size(wavelengthVec));
[~, ind500nm]           = min(abs(wavelengthVec - 500));
nm500LightAmp(ind500nm) = 1;

% Get the relative response strength of each cone subtype
[sConeRespWL, mConeRespWL, lConeRespWL]    = getConeResp(wavelengthVec, whiteLightAmp);
[sConeResp500, mConeResp500, lConeResp500] = getConeResp(wavelengthVec, nm500LightAmp);


% Plot response histogram
figure(1)
subplot(2,2,1)
plot(wavelengthVec, whiteLightAmp, 'k-')
hold on, box off, axis square
xlabel('Wavelength (nm)')
ylabel('Light amplitude')

subplot(2,2,2)
bar([1 2 3], [sConeRespWL, mConeRespWL, lConeRespWL])
hold on, box off, axis square
xlabel('Cone type (s, m, l)')
ylabel('Response (a.u.)')

subplot(2,2,3)
plot(wavelengthVec, nm500LightAmp, 'k-')
hold on, box off, axis square
xlabel('Wavelength (nm)')
ylabel('Light amplitude')

subplot(2,2,4)
bar([1 2 3], [sConeResp500, mConeResp500, lConeResp500])
hold on, box off, axis square
xlabel('Cone type (s, m, l)')
ylabel('Response (a.u.)')


% L4
% Load image
natImage = double(imread('apple.jpeg'))/255;

% Filter Constants
imSizeDeg   = 2;
pixSizeDeg  = .01;
filtPrefSf  = 1.5;
filtPrefOri = [0 pi/2];
filtAR      = 2;
filtDO      = 2;

% Specify two Gaussian derivative filters
[filtV] = giveOriFilt(imSizeDeg, pixSizeDeg, filtPrefSf, filtPrefOri(1), filtDO, filtAR);
[filtH] = giveOriFilt(imSizeDeg, pixSizeDeg, filtPrefSf, filtPrefOri(2), filtDO, filtAR);

% Compute filter responses for image center
patchRows    = round(size(natImage, 1)/2) - floor(size(filtV, 1)/2) : round(size(natImage, 1)/2) + floor(size(filtV, 1)/2); 
patchColumns = round(size(natImage, 2)/2) - floor(size(filtV, 2)/2) : round(size(natImage, 2)/2) + floor(size(filtV, 2)/2); 
framePatch   = natImage(patchRows, patchColumns);
respFiltV    = squeeze(sum(sum((framePatch - mean(framePatch(:))).*filtV, 1), 2));
respFiltH    = squeeze(sum(sum((framePatch - mean(framePatch(:))).*filtH, 1), 2));

% Plot results
figure(2)
subplot(2,2,1)
imagesc(filtV), colormap('gray'), hold on, box off, axis off, axis square

subplot(2,2,2)
imagesc(filtH), colormap('gray'), hold on, box off, axis off, axis square

subplot(2,2,3)
imagesc(natImage(patchRows, patchColumns)), colormap('gray'), hold on, box off, axis off, axis square
title('Central image patch')

subplot(2,2,4)
bar([1, 2], [respFiltV, respFiltH])
hold on, box off, axis square
xlabel('Gaussian derivative filter (V, H)')
ylabel('Response (a.u.)')



