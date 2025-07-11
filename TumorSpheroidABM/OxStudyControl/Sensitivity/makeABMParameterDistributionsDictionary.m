function D = makeABMParameterDistributionsDictionary()

% defines distributions on each of the varied abm parameters
D = dictionary();

D("carrying_capacity") = makedist("Uniform",500,1500);
D("occmax_2d") = makedist("Uniform",4,7); % cannot make a discrete uniform distribution object, so the draw from this will be floored to produce on integer
D("move_rate_microns") = makedist("Uniform",0,20);
D("g1_to_s") = makedist("Uniform",24/11 * 0.9,24/11 * 1.1);
D("s_to_g2") = makedist("Uniform",2.7,3.3);
D("g2_to_m") = makedist("Uniform",5.4,6.6);
D("m_to_g1") = makedist("Uniform",21.6,26.4);
