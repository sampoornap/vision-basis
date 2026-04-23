function [SSE, respRelRef, respRelCan] = giveSSE_coneResp(params, wavelengthVec, referenceAmp)

% Create a candidate metamer stimulus
indComp               = max(1, min(length(wavelengthVec), round(params)));
candidateAmp          = zeros(size(wavelengthVec));
candidateAmp(indComp) = 1;

% Get the response strength of each cone subtype for the reference
[sConeRespRef, mConeRespRef, lConeRespRef] = getConeResp(wavelengthVec, referenceAmp);

% Get the response strength of each cone subtype for the candidate metamer
[sConeRespCan, mConeRespCan, lConeRespCan] = getConeResp(wavelengthVec, candidateAmp);

% Compute the relative response distributions
respRelRef = [sConeRespRef, mConeRespRef, lConeRespRef]/lConeRespRef; 
respRelCan = [sConeRespCan, mConeRespCan, lConeRespCan]/lConeRespCan; 

% Compute the SSE
SSE = sum((respRelRef - respRelCan).^2);
