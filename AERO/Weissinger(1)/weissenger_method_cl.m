% Metodo di Weissenger per calcolare la circolazione e il C_L/alpha
function [gamma, cl_alpha] = weissenger_method_cl(b, c, S, S_tail, alpha_range)
    % Numero di segmenti
    n = 50; 
    dx = b / n; % Passo in x
    x = linspace(-b/2, b/2, n); % Coordinate sui segmenti

    % Nuova pendenza della curva C_L / alpha considerando i piani di coda
    cl_alpha = (2 * pi * b^2 / S) / (1 + (b^2 / S) * (S_tail / S));

    % Calcolare la circolazione su ciascun segmento
    gamma = zeros(length(alpha_range), n); % Cambiato per avere una matrice
    
    for j = 1:length(alpha_range)
        for i = 1:n
            % Calcolo della circolazione per ogni angolo di attacco
            gamma(j, i) = cl_alpha * alpha_range(j) * dx; % Gamma dipende da alpha
        end
    end
end