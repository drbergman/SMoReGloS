function out = customSolveSM(~, p, tt, ~, D, condition_on_previous, resample_t)

out = computeTimeSeries(p, tt, D.A, condition_on_previous, resample_t);

