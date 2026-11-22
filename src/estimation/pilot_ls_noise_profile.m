function noiseByPilot = pilot_ls_noise_profile(thermalNoiseVariance,pilotPower,iciVarianceByPilot)
%PILOT_LS_NOISE_PROFILE Convert received thermal noise plus LS-domain ICI
% variance into the equivalent per-pilot LS observation-noise variance.
%
% Thermal noise is present before pilot division, so its LS variance is
% thermalNoiseVariance/pilotPower. Pilot-innovation ICI estimates are already
% expressed in channel-estimate (LS) units and must NOT be divided by pilot
% power a second time.

if ~isscalar(thermalNoiseVariance) || ~isreal(thermalNoiseVariance) || ...
        ~isfinite(thermalNoiseVariance) || thermalNoiseVariance < 0
    error('pilot_ls_noise_profile:ThermalNoise','Thermal noise variance must be nonnegative and finite.');
end
if ~isscalar(pilotPower) || ~isreal(pilotPower) || ~isfinite(pilotPower) || pilotPower <= 0
    error('pilot_ls_noise_profile:PilotPower','Pilot power must be positive and finite.');
end
if ~isnumeric(iciVarianceByPilot) || ~isreal(iciVarianceByPilot) || ...
        any(~isfinite(iciVarianceByPilot(:))) || any(iciVarianceByPilot(:)<0)
    error('pilot_ls_noise_profile:ICIVariance','ICI variances must be real, finite and nonnegative.');
end
noiseByPilot = thermalNoiseVariance/pilotPower + double(iciVarianceByPilot);
end
