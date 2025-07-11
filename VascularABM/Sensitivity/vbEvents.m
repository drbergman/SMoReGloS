function [value,isterminal,direction] = vbEvents(~,y)

value = [y-1e-2;1e10-y]; % [close enough to count as 0; big enough to just compute the carrying capacity rather than let it blow up to Inf]
isterminal = true(2,1);
direction = [0;0];
