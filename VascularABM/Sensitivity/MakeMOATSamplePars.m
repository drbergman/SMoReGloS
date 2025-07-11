% This script will create the MOAT sample points
clearvars;

nfacs = 4; % number of ABM parameters
npoints = 15; % number of points for an LHS sample

par_names = ["p_{div}";"s_{div}";"r_{mig}";"p_{lim}"];

delta=1/npoints;
coord = (0.5*delta):delta:1;
points = zeros(npoints,nfacs+1,nfacs);
for j=1:nfacs
    points(:,nfacs+1,j) = coord(randperm(length(coord)));
end

for i = 1:npoints
    for j = 1:nfacs
        points(i,j,:) = points(i,nfacs+1,:);
        if points(i,nfacs+1,j)<=0.5
            points(i,j,j) = points(i,nfacs+1,j) + 0.5;
        else
            points(i,j,j) = points(i,nfacs+1,j) - 0.5;
        end
    end
end

points = reshape(points,[],nfacs);

D = makeMOATDistributions(par_names);

for i = 1:size(points,1)
    for j = 1:nfacs
        points(i,j) = icdf(D(par_names(j)),points(i,j));
    end
end

points(:,4) = round(points(:,4));

save("data/MOATLHSSample.mat","points")