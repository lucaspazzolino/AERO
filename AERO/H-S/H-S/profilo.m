clc
close all
clear

% Script MATLAB per calcolare la derivata prima della linea media

% Nomi dei file contenenti i dati
dorso_file = 'dorso.dat';
ventre_file = 'ventre.dat';

% Controllo se i file esistono
if ~isfile(dorso_file)
    error('Il file %s non esiste!', dorso_file);
end

if ~isfile(ventre_file)
    error('Il file %s non esiste!', ventre_file);
end

% Leggere i dati dal file del dorso
dorso_data = load(dorso_file);
x_dorso = dorso_data(:, 1);
y_dorso = dorso_data(:, 2);

% Leggere i dati dal file del ventre
ventre_data = load(ventre_file);
x_ventre = ventre_data(:, 1);
y_ventre = ventre_data(:, 2);

% Verifica che le coordinate 'x' del dorso e del ventre siano allineate
if any(abs(x_dorso - x_ventre) > 1e-3)  % Controllo della differenza tra x
    warning('Le coordinate x del dorso e del ventre non sono allineate. Interpolazione consigliata.');
    
    % Interpolazione dei dati per allineare le coordinate x
    x_allineato = linspace(min(x_dorso), max(x_dorso), max(length(x_dorso), length(x_ventre)));
    y_dorso_interp = interp1(x_dorso, y_dorso, x_allineato, 'linear', 'extrap');
    y_ventre_interp = interp1(x_ventre, y_ventre, x_allineato, 'linear', 'extrap');
    x_profile = x_allineato;
else
    % Se sono allineati, si possono unire direttamente i dati
    x_profile = x_dorso;
    y_dorso_interp = y_dorso;
    y_ventre_interp = y_ventre;
end



% Calcolare la linea media
y_media = (y_dorso_interp + y_ventre_interp) / 2;

% Visualizzare il profilo alare
figure;
hold on;
plot(x_profile, y_dorso_interp, 'b-', 'LineWidth', 1.5); % Dorso
plot(x_profile, y_ventre_interp, 'r-', 'LineWidth', 1.5); % Ventre
plot(x_profile, y_media, 'k--', 'LineWidth', 2); % Linea media
grid on;
xlabel('x (Posizione orizzontale)');
ylabel('y (Posizione verticale)');
title('Profilo Alare con Linea Media');
legend('Dorso', 'Ventre', 'Linea Media');
axis equal;
hold off;




% Opzione 1: Utilizzare la Spline (interpolazione cubica)
spline_fit = spline(x_profile, y_media);

% Calcolare la derivata prima della spline
spline_deriv = ppDer(spline_fit);

% Opzione 2: Utilizzare Polyfit (approssimazione polinomiale)
grado = 5; % Grado del polinomio da adattare
poly_coeff = polyfit(x_profile, y_media, grado);

% Calcolare la derivata prima del polinomio
poly_deriv_coeff = poly_coeff(1:end-1) .* (grado:-1:1);  % Derivata del polinomio
poly_deriv = @(x) polyval(poly_deriv_coeff, x);  % Funzione della derivata

% Calcolo della derivata prima usando le differenze finite
h = mean(diff(x_profile));  % Passo medio tra i punti
diff_deriv = (y_media(3:end) - y_media(1:end-2)) / (2 * h);  % Differenze finite centrali
x_diff = x_profile(2:end-1);  % Coordinate corrispondenti alle derivate

% Visualizzare la derivata prima usando la Spline
figure;
fplot(@(x) ppval(spline_deriv, x), [min(x_profile), max(x_profile)], 'g-', 'LineWidth', 2);
xlabel('x');
ylabel('y'' (Derivata prima)');
title('Derivata Prima della Linea Media (Spline)');

% Visualizzare la derivata prima usando il polinomio
figure;
fplot(poly_deriv, [min(x_profile), max(x_profile)], 'm-', 'LineWidth', 2);
xlabel('x');
ylabel('y'' (Derivata prima)');
title('Derivata Prima della Linea Media (Polinomiale)');

% Visualizzare la derivata prima usando le differenze finite
figure;
plot(x_diff, diff_deriv, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('y'' (Derivata prima)');
title('Derivata Prima della Linea Media (Differenze Finite)');


% Parametri
L = 10; % Lunghezza totale del profilo (modificabile)

% Creazione dei punti s (evitando i bordi s = ±L/2)
s = linspace(-L/2 + 1e-6, L/2 - 1e-6, 1000);  % Piccolo offset per evitare la singolarità

% Funzione per la derivata prima della spline
deriv_spline = @(x) ppval(spline_deriv, x);

% Funzione per la derivata prima del polinomio
deriv_polynomial = @(x) poly_deriv(x);

% Funzione per la derivata prima tramite differenze finite
deriv_diff_fin = @(x) interp1(x_diff, diff_deriv, x, 'linear', 'extrap');

% Funzione da integrare per la spline
integrand_spline = @(s) (1/pi) * deriv_spline(s) ./ sqrt((L^2 / 4) - s.^2);

% Funzione da integrare per il polinomio
integrand_polynomial = @(s) (1/pi) * deriv_polynomial(s) ./ sqrt((L^2 / 4) - s.^2);

% Funzione da integrare per le differenze finite
integrand_diff_fin = @(s) (1/pi) * deriv_diff_fin(s) ./ sqrt((L^2 / 4) - s.^2);

% Metodo dei trapezi per calcolare l'integrale
function integral = trapezoidal_rule(f, x)
    % Calcolare il passo
    h = x(2) - x(1);  
    % Somma ai bordi
    sum_f = f(x(1)) + f(x(end));  
    % Somma agli interni
    sum_f = sum_f + 2 * sum(f(x(2:end-1)));  
    % Applicare la formula dei trapezi
    integral = (h / 2) * sum_f;  
end

% Calcolare l'integrale numerico per la spline (Metodo dei Trapezi)
integral_spline_trapezoidal = trapezoidal_rule(integrand_spline, s);

% Calcolare l'integrale numerico per il polinomio (Metodo dei Trapezi)
integral_polynomial_trapezoidal = trapezoidal_rule(integrand_polynomial, s);

% Calcolare l'integrale numerico per le differenze finite (Metodo dei Trapezi)
integral_diff_fin_trapezoidal = trapezoidal_rule(integrand_diff_fin, s);

% Metodo di Integrazione Adattativa con 'integral' di MATLAB (gestisce singolarità)
integral_spline_integral = integral(integrand_spline, -L/2, L/2, 'RelTol', 1e-6);
integral_polynomial_integral = integral(integrand_polynomial, -L/2, L/2, 'RelTol', 1e-6);
integral_diff_fin_integral = integral(integrand_diff_fin, -L/2, L/2, 'RelTol', 1e-6);

% Visualizzare i risultati
fprintf('Integrale usando la spline (Metodo dei Trapezi): %.4f\n', integral_spline_trapezoidal);
fprintf('Integrale usando il polinomio (Metodo dei Trapezi): %.4f\n', integral_polynomial_trapezoidal);
fprintf('Integrale usando le differenze finite (Metodo dei Trapezi): %.4f\n', integral_diff_fin_trapezoidal);

fprintf('Integrale usando la spline (Metodo integrale): %.4f\n', integral_spline_integral);
fprintf('Integrale usando il polinomio (Metodo integrale): %.4f\n', integral_polynomial_integral);
fprintf('Integrale usando le differenze finite (Metodo integrale): %.4f\n', integral_diff_fin_integral);
