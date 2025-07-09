function [M,event_out] = updateTumor(M)

% selects tumor events and then performs them

[M,event_out] = eventSelection(M);

M = performEvents(M,event_out);
