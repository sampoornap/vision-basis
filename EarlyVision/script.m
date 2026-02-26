clear; close all; clc;


imSizeDeg  = 2.0;
pixSizeDeg = 0.01;
dOrder     = 2;
aRatio     = 3;

% L1
prefSf  = 2;
prefOri = deg2rad(45);
filt = giveOriFilt(imSizeDeg, pixSizeDeg, prefSf, prefOri, dOrder, aRatio);

figure; imagesc(filt); axis image; colormap gray; colorbar;
title('L1: Example oriented filter');


orisDeg = [0 45 90];
orisRad = deg2rad(orisDeg);

f_oris = cell(3,1);
for i = 1:3
    f_oris{i} = giveOriFilt(imSizeDeg, pixSizeDeg, prefSf, orisRad(i), dOrder, aRatio);
end

figure;
for i = 1:3
    subplot(1,3,i);
    imagesc(f_oris{i}); axis image; colormap gray;
    title(sprintf('L2: ori %d°', orisDeg(i)));
end



sfs = [1 2 4];
prefOri_fixed = deg2rad(45);

f_sfs = cell(3,1);
for i = 1:3
    f_sfs{i} = giveOriFilt(imSizeDeg, pixSizeDeg, sfs(i), prefOri_fixed, dOrder, aRatio);
end

figure;
for i = 1:3
    subplot(1,3,i);
    imagesc(f_sfs{i}); axis image; colormap gray;
    title(sprintf('L3: sf %d', sfs(i)));
end


I = imread('apple.jpeg');

I = im2double(I);

figure;
imagesc(I); axis image off; colormap gray; colorbar;

disp(['Image min intensity: ', num2str(min(I(:)))]);
disp(['Image max intensity: ', num2str(max(I(:)))]);

% creating image where every pixel is gaussian white noise

noiseImg = randn(size(I));
figure;
imagesc(noiseImg); axis image off; colormap gray; colorbar;
title('white noise image');

disp(['Noise mean: ', num2str(mean(noiseImg(:)))]);
disp(['Noise std: ', num2str(std(noiseImg(:)))]);

% normalization 
normf = @(h) h / sqrt(sum(h(:).^2) + 1e-12);


for i = 1:3 % take each neuron
    h = normf(f_oris{i}); % normalize the neuron
    
    % move the receptive field across the image and check how well it matches at every location
    filtered_img   = conv2(I, h, 'same');
    filtered_noise = conv2(noiseImg, h, 'same');

    subplot(3,2,(i-1)*2+1);
    imagesc(filtered_img); axis image off; colormap gray; 
    
end


% measure orientation tuning

h = normf(filt);                 % normalize filter for comparable responses
dims = size(h);                  % stimulus size = filter size
sfFixed = prefSf;                % keep SF fixed at the preferred SF
phase = 0;                       

orisTestDeg = 0:10:170;          % test orientations
orisTestRad = deg2rad(orisTestDeg);

respOri = zeros(size(orisTestRad));

for i = 1:numel(orisTestRad)
    stim = giveSine(dims, pixSizeDeg, sfFixed, orisTestRad(i), phase);
    stim = stim - mean(stim(:));             % remove DC
    respOri(i) = sum(stim(:) .* h(:));       
end

figure;
plot(orisTestDeg, respOri, 'LineWidth', 2);
xlabel('Orientation (deg)');
ylabel('Response');
title('L5: Orientation tuning');
grid on;

% measure spatial frequency tuning
oriFixed = prefOri;              % keep orientation fixed at preferred orientation
phase = 0;

sfsTest = 0.5:0.25:6;            % test SF range (cycles/deg)
respSf = zeros(size(sfsTest));

for i = 1:numel(sfsTest)
    stim = giveSine(dims, pixSizeDeg, sfsTest(i), oriFixed, phase);
    stim = stim - mean(stim(:));             % remove DC
    respSf(i) = sum(stim(:) .* h(:));        
end

figure;
plot(sfsTest, respSf, 'LineWidth', 2);
xlabel('Spatial frequency (cycles/deg)');
ylabel('Response');
title('L6: Spatial frequency tuning');
grid on;
