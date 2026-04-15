%% SpikeTrainFT
function Y = SpikeTrainFT(Tau, f, T)
% Fourier transfom of an idealized spike train of delta functions
% at times tau (seconds), evaluated at frequencies in f (Hz)
% Division by duration T give units of impulses /second.

if iscell(Tau)
    for i = 1:numel(Tau)
        tau = Tau{i}(:);
        Y(:,i) = sum(exp(-2*pi*sqrt(-1).* tau*f(i)) ,1)./T(i);
    end
else
    tau = Tau(:);
    Y(:,1) = sum(exp(-2*pi*sqrt(-1).* tau*f) ,1)./T;
end
