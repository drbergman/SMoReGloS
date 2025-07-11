function sm_data = customSolveSM(sm,p,tt,~,~,~,resample_t)

if isempty(resample_t)
    sm_data = computeTimeSeries(p, tt, [], sm.opts, []);
else
    sm_data = computeTimeSeries(p, resample_t, [], sm.opts, []);
end
