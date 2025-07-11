function C = convertParNames(par_names_ordered)

C = strings(1,length(par_names_ordered));

for i = 1:length(par_names_ordered)
    switch par_names_ordered(i)
        case "carrying_capacity"
            C(i) = "K_A";
        case "g1_to_s"
            C(i) = "\rho_{G1\rightarrowS}";
        case "s_to_g2"
            C(i) = "\rho_{S\rightarrowG2}";
        case "g2_to_m"
            C(i) = "\rho_{G2\rightarrowM}";
        case "m_to_g1"
            C(i) = "\rho_{M\rightarrowG1}";
        case "move_rate_microns"
            C(i) = "s";
        case "occmax_2d"
            C(i) = "T_{con}";
    end
end
