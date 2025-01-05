% Metodo di Weissenger per calcolare la circolazione e il C_L/alpha
function [gamma, cl_alpha] = weissenger_method(b, c, S, alpha_range)
    % Numero di segmenti
    n = 50; 
    dx = b / n; % Passo in x
    x = linspace(-b/2, b/2, n); % Coordinate sui segmenti
    cl_alpha = 2 * pi * b / S; % Pendenza della curva C_L / alpha per un'ala ellittica

    % Calcolare la circolazione su ciascun segmento
    gamma = zeros(length(alpha_range), n); % Cambiato per avere una matrice
    
    for j = 1:length(alpha_range)
        for i = 1:n
            % Calcolo della circolazione per ogni angolo di attacco
            gamma(j, i) = cl_alpha * alpha_range(j) * dx; % Gamma dipende da alpha
        end
    end
end
