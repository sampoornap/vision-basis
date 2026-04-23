clear all
close all
clc

% Set constants
nFrames  = 11;
stepSize = 15; % in pixels

% Load image
natImage       = double(imread('apple.jpeg'))/255;

% Make video
[nRows, nCols] = size(natImage);
shiftVec       = 1:stepSize:nFrames*stepSize;

for iF = 1:nFrames
    frameMat(:,:,iF) = natImage(shiftVec(iF):nRows-shiftVec(nFrames + 1 - iF), shiftVec(iF):nCols-shiftVec(nFrames + 1 - iF));

    % Show video
    imagesc(frameMat(:,:,iF)), colormap('gray'), hold on, box off, axis off, axis square
    pause(.1)
end

