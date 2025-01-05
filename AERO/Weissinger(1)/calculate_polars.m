% Calcolare le polari per C_L e C_D
function [cl, cd] = calculate_polars(gamma, b, c, AR)
    % Calcolo di C_L e C_D
    cl = 2 * pi * gamma / (b * c); % Coefficiente di portanza
    cd = (cl.^2) / (pi * AR); % Coefficiente di resistenza indotto
end