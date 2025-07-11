function out = summarizeSMLHS(in)

% computes the geometric mean of [x_i] where x_i = exp(in)
% geometric mean is used to avoid exponentiating each element of in, some
% of which are expected to be large enough to return Inf upon exponentiation

out = exp(mean(in));