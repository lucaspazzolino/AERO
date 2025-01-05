% Calcolare le polari per C_L e C_D
function [cl, cd] = calculate_polars_tail(gamma, b, c, AR, S_tail, AR_tail)
    % Calcolo di C_L e C_D considerando anche il piano di coda
    cl = 2 * pi * gamma / (b * c); % Coefficiente di portanza
    cd = (cl.^2) / (pi * AR); % Coefficiente di resistenza indotto

    % Aggiungere resistenza indotta del piano di coda
    cd_tail = (cl.^2) / (pi * AR_tail); % Resistenza indotta del piano di coda
    cd = cd + cd_tail; % Somma della resistenza
end
