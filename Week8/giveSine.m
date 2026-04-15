function [sineWave] = giveSine(dims, pixSizeDeg, sf, ori, phase)

% Compute image
x = linspace(-dims(1)/2, dims(1)/2, dims(1));
y = linspace(-dims(2)/2, dims(2)/2, dims(2));
z = cos(ori) * repmat(y, [dims(1) 1]) + sin(ori) * repmat((1:dims(1))' - (dims(1) + 1)/2, [1 dims(2)]);

sineWave = .5+(.5*(cos(2*pi*sf*pixSizeDeg*z+phase)));
end
