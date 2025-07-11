% plots histograms of the sensitivities of the parameters in both models

clearvars;

addpath("../../../src/myfunctions/") % replace with path (rel or abs) to myfunctions

[~,~,~] = mkdir("figures/fig");
[~,~,~] = mkdir("figures/png");

direct_method_color = 0.1*[1 1 1];
sm_model_palette = containers.Map("sm",lines(1));

ABM = load("data/MOATDirect.mat","par_names","order","mu_star");
ABM_SM = load("data/MOATIndirect.mat","mu_star");

C = convertParNames(ABM.par_names(ABM.order));
y = [ABM.mu_star;ABM_SM.mu_star];

f=figureOnRight;
hold on;
rectangle('Position',[1.5, 0, 2, 350], 'FaceColor', 'yellow', "EdgeColor", "none")
ax = gca;
b = bar(1:length(C),y);
xticks(ax, 1:length(C))
xticklabels(ax, C)
b(1).FaceColor = direct_method_color;
b(2).FaceColor = sm_model_palette("sm");
legend({"Direct","Indirect"},"location","northeast","FontSize",8)

%% 
f.Units = "inches";
f.Position(3:4) = [2.5,1];
set(ax,'FontSize',8)

%% set margins
margin = struct("left",.16,"right",.01,"top",.05,"bottom",.33);
spacing = struct("horizontal",.01,"vertical",.09);
uniformAxisSpacing(ax,margin,spacing);

%%
savefig("figures/fig/CompareABMSensitivityMethods")
print("-r300","figures/png/CompareABMSensitivityMethods","-dpng")

%%
y = y./sum(y,2);
figure;
bar(C,y)
ylabel("\mu^* / <\mu^*>")
legend({"Sensitivity on ABM directly","Sensitivity on ABM via SM"},"location","northeast","FontSize",16)
set(gca,'FontSize',16)

%%
savefig("figures/fig/CompareABMSensitivityMethods_NormalizedMu")
print("-r300","figures/png/CompareABMSensitivityMethods_NormalizedMu","-dpng")

%%
f=figureOnRight;
ax = gca;
s = ["Direct","Indirect"];
B = bar(1:2,y,"stacked");
B(1).BarWidth = 0.8;
B(2).BarWidth = 0.8;
ylabel("\mu^* / <\mu^*>")
xticklabels(s)
legend(flip(B),flip(C),"location","bestoutside","FontSize",8)
set(ax,"FontSize",8)

%% 
f.Units = "inches";
f.Position(3:4) = [2.5,1.5];
set(ax,'FontSize',8)

%%
xlim([0.4 2.6])

%% set margins
margin = struct("left",.17,"right",.4,"top",.04,"bottom",.1);
spacing = struct("horizontal",.01,"vertical",.09);
uniformAxisSpacing(ax,margin,spacing);

%%
savefig("figures/fig/CompareABMSensitivityMethods_NormalizedMu_Stacked")
print("-r300","figures/png/CompareABMSensitivityMethods_NormalizedMu_Stacked","-dpng")
