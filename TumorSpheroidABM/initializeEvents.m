function M = initializeEvents(M)

M.events.times = M.setup.censor_date;
M.events.event_index = Inf;

[M.events.times,order] = sort(M.events.times,1,"ascend");
M.events.event_index = M.events.event_index(order);

M.events.n = length(M.events.times);