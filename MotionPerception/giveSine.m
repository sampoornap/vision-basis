function [sineWave] = giveSine(dims, pixSizeDeg, sf, ori, phase)

% Convert orientation to radians
ori = deg2rad(ori);

% Create coordinate grid
[x, y] = meshgrid( ...
    ((1:dims(2)) - (dims(2)+1)/2), ...
    ((1:dims(1)) - (dims(1)+1)/2));

% Convert pixels → degrees
x = x * pixSizeDeg;
y = y * pixSizeDeg;

% Oriented grating
z = x * cos(ori) + y * sin(ori);

% Generate sine wave
sineWave = 0.5 + 0.5 * cos(2 * pi * sf * z + phase);

end