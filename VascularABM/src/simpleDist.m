function out = simpleDist(u,v)

out = sqrt(sum((u-v).^2,2));
