function identifiability_index_palette = identifiabilityIndexPalette()

identifiability_index_palette = containers.Map();
% zero_color = [2,204,254]/255;
% one_color = [231,82,135]/255;
% two_color = [120,81,169]/255;
% zero_color = [65,105,225]/255;
% one_color = [255,0,255]/255;
% two_color = [255,215,0]/255;
% zero_color = [152,255,152]/255;
% one_color = [255,218,185]/255;
% two_color = [200,162,200]/255;

% M = [ 64, 224, 208 ;
%    255, 255, 85 ;
%    255, 105, 180 ]/255;
% M = [ 33, 33, 33 ;
%    25, 25, 112 ;
%    139, 0, 0 ]/255;
% M = [
%     10, 20, 40;   % Shadowy Slate
%     0, 30, 15;    % Deep Forest Green
%     20, 0, 60      % Midnight Blue
% ]/255;

M = [
    0, 129, 109;   
    142, 10, 38;     % Maroon (University of Minnesota Duluth)
    0, 46, 98       % Sable (Johns Hopkins University)
]/255;

% identifiability_index_palette(0) = zero_color;
identifiability_index_palette("0") = M(1,:);

% identifiability_index_palette(1) = one_color;
identifiability_index_palette("1") = M(2,:);

% identifiability_index_palette(2) = two_color;
identifiability_index_palette("2") = M(3,:);