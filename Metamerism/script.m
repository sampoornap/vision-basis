% L1: Create a color metamer for white light (use script class7Script)
% L2: Create a color metamer for 550 nm light
% L3: Repeat both exercises for a deuteranope (i.e., a person without M-cones)
% L4: Create a metamer for a natural image processed by a pair of linear filters (line 51–90)

wavelengthVec = linspace(350, 750, 401);
whiteLightAmp = ones(size(wavelengthVec));

% 550 nm light
nm550LightAmp = zeros(size(wavelengthVec));
[~, ind550] = min(abs(wavelengthVec - 550));
nm550LightAmp(ind550) = 1;

% L1
bestParam_white = fminsearch(@(p) giveSSE_coneResp(p, wavelengthVec, whiteLightAmp), params0);

[~, respRef_WL, respCan_WL] = giveSSE_coneResp(bestParam_white, wavelengthVec, whiteLightAmp);

disp('L1: best wavelength index for white light metamer');
disp(bestParam_white);

% L2

bestParam_550 = fminsearch(@(p) giveSSE_coneResp(p, wavelengthVec, nm550LightAmp), params0);

[~, respRef_550, respCan_550] = giveSSE_coneResp(bestParam_550, wavelengthVec, nm550LightAmp);

disp('L2: best wavelength index for 550 nm metamer');
disp(bestParam_550);


% L3 

% first modify mResp to be = 0

bestParam_white_deut = fminsearch(@(p) giveSSE_coneResp(p, wavelengthVec, whiteLightAmp), params0);
bestParam_550_deut   = fminsearch(@(p) giveSSE_coneResp(p, wavelengthVec, nm550LightAmp), params0);

disp('Deuteranope metamer (white):'); disp(bestParam_white_deut);
disp('Deuteranope metamer (550nm):'); disp(bestParam_550_deut);

% L4 - image metamer 

natImage = im2double(imread('apple.jpeg'));
imSizeDeg   = 2;
pixSizeDeg  = 0.01;
filtPrefSf  = 1.5;
filtPrefOri = [0 pi/2];
filtAR      = 2;
filtDO      = 2;


filtV = giveOriFilt(imSizeDeg, pixSizeDeg, filtPrefSf, filtPrefOri(1), filtDO, filtAR);
filtH = giveOriFilt(imSizeDeg, pixSizeDeg, filtPrefSf, filtPrefOri(2), filtDO, filtAR);

% extract center patch
patchRows = round(size(natImage,1)/2) - floor(size(filtV,1)/2) : round(size(natImage,1)/2) + floor(size(filtV,1)/2);

patchCols = round(size(natImage,2)/2) - floor(size(filtV,2)/2) : round(size(natImage,2)/2) + floor(size(filtV,2)/2);

framePatch = natImage(patchRows, patchCols);

% 3 sine components sf, ori, phase, amplitude
params0 = rand(3,4);

% flatten for optimizer
params0_vec = params0(:);

% optimization
bestParams_vec = fminsearch(@(p) giveSSE_OriFilt(reshape(p,3,4), framePatch, pixSizeDeg, filtV, filtH), params0_vec);

bestParams = reshape(bestParams_vec, 3, 4);

% generating metamer
[~, metamer] = giveSSE_OriFilt(bestParams, framePatch, pixSizeDeg, filtV, filtH);



% visualization
figure;
subplot(1,2,1);
imagesc(framePatch); colormap gray; axis image off;
title('Original Patch');

subplot(1,2,2);
imagesc(metamer); colormap gray; axis image off;
title('Metamer (same filter responses)');