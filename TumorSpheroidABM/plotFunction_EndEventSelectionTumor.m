function plotFunction_EndEventSelectionTumor(M, prob_matrix)

%% tumor probabilities
order = [3,2,1]; % first, separate by fast vs slow, then engaged/non-engaged
xs = linspace(0,1,size(prob_matrix,1));
ys = sortrows(prob_matrix(:,order));
for i = 1:size(prob_matrix,2)
    set(M.fig.tum_prob_bar(i),'XData',xs,'YData',ys(:,i))
end
if length(xs)>1
    M.fig.ax(M.fig.tum_prob_ind).XLim = [0,1] + .5*[-1,1]*xs(2);
end
