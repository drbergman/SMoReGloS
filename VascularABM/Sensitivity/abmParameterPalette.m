function abm_par_palette = abmParameterPalette()

% M = [
%     0, 31, 63;     % Dark Blue (#001F3F)
%     46, 204, 113;   % Emerald Green (#2ECC71)
%     231, 76, 60;    % Crimson Red (#E74C3C)
%     255, 215, 0     % Golden Yellow (#FFD700)
% ] / 255;

M = lines(7);
M = M(4:end,:);
abm_par_names = ["p_{lim}","s_{div}","p_{div}","r_{mig}"];

abm_par_palette = containers.Map();
for i = 1:length(abm_par_names)
    abm_par_palette(abm_par_names(i)) = M(i,:);
end