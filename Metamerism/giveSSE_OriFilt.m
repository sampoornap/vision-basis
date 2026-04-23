function [SSE, candPatch] = giveSSE_OriFilt(params, framePatch, pixSizeDeg, filtV, filtH)


% Create a candidate metamer stimulus
for iC = 1:size(params, 1)
    sinMat(:,:,iC) = params(iC,4) * (giveSine(size(framePatch), pixSizeDeg, params(iC,1), params(iC,2), params(iC,3)) -.5) + .5;
end
candPatch = framePatch - mean(framePatch(:)) + sum(sinMat - .5, 3) + .5;

% Keep contrast in check
lumRange  = max(candPatch(:)) - min(candPatch(:));
candPatch = (candPatch - mean(candPatch(:)))/lumRange + .5;

% Get the response strength of the filters for the reference
respFiltVRef = squeeze(sum(sum((framePatch - mean(framePatch(:))).*filtV, 1), 2));
respFiltHRef = squeeze(sum(sum((framePatch - mean(framePatch(:))).*filtH, 1), 2));

% Get the response strength of the filters for the candidate metamer
respFiltVCan = squeeze(sum(sum((candPatch - mean(candPatch(:))).*filtV, 1), 2));
respFiltHCan = squeeze(sum(sum((candPatch - mean(candPatch(:))).*filtH, 1), 2));

% Compute the SSE
SSE = sum(([respFiltVRef respFiltHRef] - [respFiltVCan respFiltHCan]).^2);

