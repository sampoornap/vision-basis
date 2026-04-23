function [sResp, mResp, lResp] = getConeResp(wavelength, profile);

sSel = normpdf(wavelength, 445, 20)/normpdf(445, 445, 20);
mSel = normpdf(wavelength, 545, 35)/normpdf(545, 545, 35);
lSel = normpdf(wavelength, 565, 40)/normpdf(565, 565, 40);

sResp = sum(profile.*sSel);
mResp = 0;
lResp = sum(profile.*lSel);
