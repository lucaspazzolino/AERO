clc
close all
clear 

% Carica i dati del profilo (assunto calcolato in precedenza)
data = readmatrix('rae2822_coordinates.csv');
x = data(:, 1);        % Coordinate x
y_sup = data(:, 2);    % Coordinate y dell'estradosso
y_inf = data(:, 3);    % Coordinate y dell'intradosso

% Calcolo della linea media
y_mean = (y_sup + y_inf) / 2;

% Lunghezza della corda (L)
L = 1; % Normalizzato a 1 (modifica se necessario)

% Calcolo della derivata della linea media usando differenze finite
dy_dx = gradient(y_mean, x); % Approssimazione della derivata

% Preparazione per l'integrale
N = 1000; % Numero di punti per la discretizzazione
s = linspace(-L/2, L/2, N); % Variabile di integrazione

% Calcolo dell'integrando in corrispondenza dei punti discretizzati
% Usando interpolazione per allineare dy_dx al dominio di s
dy_dx_interp = interp1(x - L/2, dy_dx, s, 'linear', 'extrap');
integrand = dy_dx_interp ./ sqrt((L^2 / 4) - s.^2);

% Applicazione del metodo dei trapezi
delta_s = s(2) - s(1); % Passo di integrazione
integral_value = sum((integrand(1:end-1) + integrand(2:end)) / 2 * delta_s);

% Output del risultato
fprintf('Il valore del richiesto e: %.6f\n', integral_value);
