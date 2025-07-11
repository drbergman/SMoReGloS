function out = smModelPalette()

model_type = ["exponential";"logistic";"von_bertalanffy"];
colors = lines(3);
out = containers.Map();
for i = 1:length(model_type)
    out(model_type(i)) = colors(i,:);
end