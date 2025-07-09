function cycle_pars = baseCyclePars()

% transition rates from one phase to the next in per day; all based on typical durations of the phases
cycle_pars.g1_to_s = 24/11;
cycle_pars.s_to_g2 = 24/8;
cycle_pars.g2_to_m = 24/4;
cycle_pars.m_to_g1 = 24/1;

% transition from arrested state back into cell cycle in per day
cycle_pars.g1a_to_g1 = 0.06;
cycle_pars.sa_to_g1 = 0.06;
cycle_pars.g2a_to_g1 = 0.06;
cycle_pars.ma_to_g1 = 0.06;
