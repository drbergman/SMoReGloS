% plots histograms of the sensitivities of the parameters in both models

clearvars;

addpath("../../src/myfunctions/") % replace with path (rel or abs) to myfunctions
addpath("../ODEFitting/")

save_figs = true;
reprint = true;
file_types = ["fig","png"];
resolution = "-r300";

direct_method_color = 0.1*[1 1 1];
sm_model_palette = smModelPalette();

% model_type = "exponential";
% model_type = "logistic";
model_type = ["exponential","logistic","von_bertalanffy"];
endpoints = ["final_size", "AUC", "time_to_half"];

for endpoint = endpoints
    legend_entries = ["Direct","Exponential","Logistic","Von Bertalanffy"];

    if endpoint == "time_to_half"
        model_type(model_type=="von_bertalanffy") = [];
        legend_entries(legend_entries=="Von Bertalanffy") = [];
    end

    abm_sensitivty_file = "data/GlobalSensitivityEFASTDirect"+ "_" + endpoint + ".mat";
    ABM = load(sprintf(abm_sensitivty_file,endpoint));
    for i = 1:length(model_type)
        sm_sensitivty_file = sprintf("data/GlobalSensitivityEFASTIndirect_%s_%s.mat",model_type(i), endpoint);
        ABM_SM(i) = load(sm_sensitivty_file,"display_par_names","S1","ST");
    end

    line_width = 1;

    %% get indirect order to match direct order
    y_all = zeros(size(ABM.S1,1),1+length(model_type),2);
    y_all(:,1,1) = ABM.S1;
    y_S1 = zeros(size(ABM.S1,1),1+length(model_type));
    [y_S1(:,1),order_abm] = sort(ABM.S1,"descend");
    C = categorical(ABM.display_par_names(order_abm),ABM.display_par_names(order_abm));
    for i = 1:length(model_type)
        y_S1(:,i+1) = ABM_SM(i).S1(order_abm);
        y_all(:,i+1,1) = ABM_SM(i).S1;
    end


    %% S1
    f=figureOnRight("Name",sprintf("CompareeFAST_S1_%s",endpoint));
    hold on;
    b=bar(y_S1,"LineWidth",line_width);
    b(1).FaceColor = direct_method_color;
    for i = 1:length(model_type)
        b(i+1).FaceColor = sm_model_palette(model_type(i));
    end
    xticks(1:size(y_S1,1))
    xticklabels(C)
    set(gca,'FontSize',16)
    ylim([0 1])
    yticks([0 1])
    saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,resolution=resolution)

    %% get indirect order to match direct order
    y_ST = zeros(size(ABM.ST,1),1+length(model_type));
    y_all(:,1,2) = ABM.ST;
    [y_ST(:,1),order_abm] = sort(ABM.ST,"descend");
    C = categorical(ABM.display_par_names(order_abm),ABM.display_par_names(order_abm));
    for i = 1:length(model_type)
        y_ST(:,i+1) = ABM_SM(i).ST(order_abm);
        y_all(:,i+1,2) = ABM_SM(i).ST;
    end

    %% ST
    f=figureOnRight("Name",sprintf("CompareeFAST_ST_%s",endpoint));
    hold on;
    b=bar(y_ST,"LineWidth",line_width);
    b(1).FaceColor = direct_method_color;
    for i = 1:length(model_type)
        b(i+1).FaceColor = sm_model_palette(model_type(i));
    end
    xticks(1:size(y_ST,1))
    xticklabels(C)
    [ngroups,nbars] = size(y_ST);
    set(gca,'FontSize',16)
    ylim([0 1])
    yticks([0 1])
    saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,resolution=resolution)

    %% S1 and ST stacked
    Y = y_all;
    [~,order_all] = sort(Y(:,1,2),"descend");
    Y = Y(order_all,:,:);
    Y(:,:,2) = Y(:,:,2) - Y(:,:,1); % get difference for correct stacking of these bars
    f=figureOnRight("Name",sprintf("CompareeFAST_Both_%s",endpoint));
    hold on;
    b  =cell(4,3);
    dw = 1/(length(model_type)+2);
    width = 1/(length(model_type)+3);
    s1_transform = @(x) nthroot(x,1);
    sT_transform = @(x) nthroot(x,1);
    sT_alpha = 0.5;
    top_bar_line_style = "None";

    bottom_bar_line_width = 0.5;
    bottom_bar_line_style = "-";
    C = categorical(ABM.display_par_names(order_all),ABM.display_par_names(order_all));
    for par_ind = 1:4
        b{par_ind,1}=bar(par_ind-dw,squeeze(Y(par_ind,1,:)),"stacked");
        b{par_ind,1}(1).BarWidth = width;
        b{par_ind,1}(2).BarWidth = width;
        b{par_ind,1}(1).FaceColor = s1_transform(direct_method_color);
        b{par_ind,1}(2).FaceColor = sT_transform(direct_method_color);
        b{par_ind,1}(2).FaceAlpha = sT_alpha;

        b{par_ind,1}(1).LineWidth = bottom_bar_line_width;
        b{par_ind,1}(1).LineStyle = bottom_bar_line_style;
        b{par_ind,1}(2).LineStyle = top_bar_line_style;

        for model_index = 1:length(model_type)
            b{par_ind,model_index+1}=bar(par_ind+dw*model_index - dw,squeeze(Y(par_ind,model_index+1,:)),"stacked");
            b{par_ind,model_index+1}(1).BarWidth = width;
            b{par_ind,model_index+1}(2).BarWidth = width;
            b{par_ind,model_index+1}(1).FaceColor = s1_transform(sm_model_palette(model_type(model_index)));
            b{par_ind,model_index+1}(2).FaceColor = sT_transform(sm_model_palette(model_type(model_index)));
            b{par_ind,model_index+1}(2).FaceAlpha = sT_alpha;

            b{par_ind,model_index+1}(1).LineWidth = bottom_bar_line_width;
            b{par_ind,model_index+1}(1).LineStyle = bottom_bar_line_style;

            b{par_ind,model_index+1}(2).LineStyle = top_bar_line_style;
        end
    end

    xticks(1:size(y_ST,1))
    xticklabels(C)
    [ngroups,nbars] = size(y_ST);
    set(gca,'FontSize',16)
    xlim([1-2*dw, par_ind+dw*model_index])
    ylim([0 1])
    yticks([0 1])
    saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,resolution=resolution)


    %% normalized mu
    y_S1_sum = sum(y_S1,1);
    y_S1_normalized = y_S1./y_S1_sum;

    %% stacked histogram
    f = figureOnRight("Name",sprintf("CompareeFAST_S1_Methods_NormalizedMu_Stacked_All_%s",endpoint));
    B = bar(categorical(legend_entries,legend_entries),y_S1_normalized',"stacked");
    % ylabel("\mu^* / <\mu^*>")
    legend(flip(B),flip(C),"location","bestoutside","FontSize",16)
    abm_par_palette = abmParameterPalette();
    for i = 1:numel(B)
        B(i).FaceColor = abm_par_palette(B(i).DisplayName);
    end
    set(gca,"FontSize",16)
    ylim([0 1])
    yticks([0 1])

    %% save figs
    saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,resolution=resolution)

    %% normalized mu
    y_ST_sum = sum(y_ST,1);
    y_ST_normalized = y_ST./y_ST_sum;

    %% stacked histogram
    f = figureOnRight("Name",sprintf("CompareeFAST_ST_Methods_NormalizedMu_Stacked_All_%s",endpoint));
    B = bar(categorical(legend_entries,legend_entries),y_ST_normalized',"stacked");
    % ylabel("\mu^* / <\mu^*>")
    legend(flip(B),flip(C),"location","bestoutside","FontSize",16)
    abm_par_palette = abmParameterPalette();
    for i = 1:numel(B)
        B(i).FaceColor = abm_par_palette(B(i).DisplayName);
    end
    set(gca,"FontSize",16)
    ylim([0 1])
    yticks([0 1])

    %% save figs
    saveFigures(f,save_figs=save_figs,reprint=reprint,file_types=file_types,resolution=resolution)

end
%% reset path
rmpath("../ODEFitting/")

