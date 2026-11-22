function [Lrx,Ltx,Rrx,Rtx] = spatial_correlation_factors(cfg)
%SPATIAL_CORRELATION_FACTORS Build Kronecker correlation matrices/factors.
% Exponential Toeplitz model: R(i,j)=rho^abs(i-j), with rho in [0,0.99].
% The returned lower factors satisfy L*L'=R (MATLAB conjugate transpose).
Nr = cfg.system.nRx;
Nt = cfg.system.nTx;
if ~isfield(cfg.channel,'spatialCorrelation') || ...
        ~isfield(cfg.channel.spatialCorrelation,'enabled') || ...
        ~cfg.channel.spatialCorrelation.enabled
    rhoTx = 0;
    rhoRx = 0;
else
    rhoTx = cfg.channel.spatialCorrelation.txRho;
    rhoRx = cfg.channel.spatialCorrelation.rxRho;
end
validateRho(rhoTx,'txRho');
validateRho(rhoRx,'rxRho');
Rtx = toeplitz(rhoTx.^(0:Nt-1));
Rrx = toeplitz(rhoRx.^(0:Nr-1));
Ltx = chol((Rtx+Rtx')/2,'lower');
Lrx = chol((Rrx+Rrx')/2,'lower');
end

function validateRho(rho,name)
if ~isscalar(rho) || ~isreal(rho) || ~isfinite(rho) || rho < 0 || rho > 0.99
    error('spatial_correlation_factors:InvalidRho', ...
        '%s must be a real scalar in [0,0.99].',name);
end
end
