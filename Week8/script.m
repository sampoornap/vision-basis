% L1: Make of movie of two drifting gratings (sf: 2 c/deg, speed: 3 deg/sec, ori: –45 and +45 deg)
% L2: Make a movie of a plaid composed of the sum of both drifting gratings
% L3: Make a movie of a plaid whose components differ in orientation by an arbitrary amount
% L4: Make a movie of a plaid whose components produce transparent instead of coherent motion
% L5: Make a movie of 3-D luminance noise
% L6: Make a velocity filter in the 3–D Fourier domain and create filtered noise

function [sineWave] = giveSine(dims, pixSizeDeg, sf, ori, phase)

% Compute image
x = linspace(-dims(1)/2, dims(1)/2, dims(1));
y = linspace(-dims(2)/2, dims(2)/2, dims(2));
z = cos(ori) * repmat(y, [dims(1) 1]) + sin(ori) * repmat((1:dims(1))' - (dims(1) + 1)/2, [1 dims(2)]);

sineWave = .5+(.5*(cos(2*pi*sf*pixSizeDeg*z+phase)));
end

sf = 2; % spatial frequency in cycles per degree
speed = 3; % speed in degrees per second
ori1 = -45; % orientation of the first grating
ori2 = 30; % orientation of the second grating
dims = [1024 1024]; 
pixSizeDeg = 0.02;
% t = 1;


fps = 30;
time = 2;
numframes = fps * time;

for f = 1:numframes
    t = (f-1) / fps; 
    phase1 = 2*pi * sf * speed * t;
    phase2 = 1.5*pi * sf * speed * t;


    grating1 = giveSine(dims, pixSizeDeg, sf, deg2rad(ori1), phase1);
    grating2 = giveSine(dims, pixSizeDeg, sf, deg2rad(ori2), phase2);

% combining them 

    img = (grating1 + grating2)/2;
    
    imagesc(img);
    colormap gray;
    axis image off;
    caxis([0 1]);
    drawnow;
end

