clc
close all
clear

% Script MATLAB per visualizzare il profilo alare unendo dorso e ventre

% Nomi dei file contenenti i dati
dorso_file = 'dorso_naca0012.dat';
ventre_file = 'ventre_NACA0012.dat';

% Leggere i dati dal file del dorso
dorso_data = load(dorso_file);
x_dorso = dorso_data(:, 1);
y_dorso = dorso_data(:, 2);

% Leggere i dati dal file del ventre
ventre_data = load(ventre_file);
x_ventre = ventre_data(:, 1);
y_ventre = ventre_data(:, 2);

% Unire i punti del dorso e del ventre
x_profile = [x_dorso; flipud(x_ventre)];
y_profile = [y_dorso; flipud(y_ventre)];

% Visualizzare il profilo alare
figure;
plot(x_profile, y_profile, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('y');
title('Profilo Alare');
axis equal;

% Calcolare la linea media (corretto)
y_midline = (y_dorso + y_ventre) / 2

% Visualizzare la linea media
hold on;
plot(x_dorso, y_midline, 'r-', 'LineWidth', 1.5);
legend('Profilo Alare', 'Linea Media');
hold off;

% Calcolare la derivata prima della linea media
dx = diff(x_dorso);  % Calcolare la differenza nei punti x
dy = diff(y_midline);  % Differenza nei valori y della linea media
deriv_diff_fin = dy ./ dx;  % Derivata con differenze finite

% Calcolare la derivata prima tramite spline
spline_deriv = spline(x_dorso, y_midline);
deriv_spline = @(x) ppval(fnder(spline_deriv), x);

% Calcolare la derivata prima tramite polinomio (fit polinomiale)
poly_coeffs = polyfit(x_dorso, y_midline, 3);
deriv_polynomial = @(x) polyval(polyder(poly_coeffs), x);

% Visualizzare le derivate
figure;
plot(x_dorso(2:end), deriv_diff_fin, 'b-', 'LineWidth', 1.5);
hold on;
fplot(deriv_spline, [min(x_dorso), max(x_dorso)], 'r--', 'LineWidth', 1.5);
fplot(deriv_polynomial, [min(x_dorso), max(x_dorso)], 'g-.', 'LineWidth', 1.5);
xlabel('x');
ylabel('Derivata della Linea Media');
legend('Derivata con Differenze Finite', 'Derivata con Spline', 'Derivata con Polinomio');
grid on;
title('Derivata Prima della Linea Media');
hold off;

% Calcolare l'integrale della derivata prima (come esempio)
L = max(x_dorso) - min(x_dorso);  % Lunghezza totale del profilo
integrand_spline = @(s) (1/pi) * deriv_spline(s) ./ sqrt((L^2 / 4) - s.^2);
integrand_polynomial = @(s) (1/pi) * deriv_polynomial(s) ./ sqrt((L^2 / 4) - s.^2);
integrand_diff_fin = @(s) (1/pi) * interp1(x_dorso(1:end-1), deriv_diff_fin, s, 'linear', 'extrap') ./ sqrt((L^2 / 4) - s.^2);

% Metodo dei trapezi per il calcolo dell'integrale
function integral = trapezoidal_rule(f, x)
    h = x(2) - x(1);  % Calcolare il passo
    sum_f = f(x(1)) + f(x(end));  % Somma ai bordi
    sum_f = sum_f + 2 * sum(f(x(2:end-1)));  % Somma agli interni
    integral = (h / 2) * sum_f;  % Applicare la formula dei trapezi
end

% Calcolare l'integrale numerico per ciascun metodo
s = linspace(min(x_dorso), max(x_dorso), 1000);  % Definire i punti per l'integrazione
integral_spline_trapezoidal = trapezoidal_rule(integrand_spline, s);
integral_polynomial_trapezoidal = trapezoidal_rule(integrand_polynomial, s);
integral_diff_fin_trapezoidal = trapezoidal_rule(integrand_diff_fin, s);

% Metodo di integrazione adattiva con 'integral' di MATLAB
integral_spline_integral = integral(integrand_spline, min(x_dorso), max(x_dorso), 'RelTol', 1e-6);
integral_polynomial_integral = integral(integrand_polynomial, min(x_dorso), max(x_dorso), 'RelTol', 1e-6);
integral_diff_fin_integral = integral(integrand_diff_fin, min(x_dorso), max(x_dorso), 'RelTol', 1e-6);

% Calcolare l'errore relativo per i vari metodi
true_value = integral_spline_integral;  % Usare il valore dell'integrale adattivo come riferimento

% Errore relativo per il metodo dei trapezi
error_spline_trapezoidal = abs(integral_spline_trapezoidal - true_value) / true_value;
error_polynomial_trapezoidal = abs(integral_polynomial_trapezoidal - true_value) / true_value;
error_diff_fin_trapezoidal = abs(integral_diff_fin_trapezoidal - true_value) / true_value;

% Errore relativo per il metodo adattivo
error_spline_integral = abs(integral_spline_integral - true_value) / true_value;
error_polynomial_integral = abs(integral_polynomial_integral - true_value) / true_value;
error_diff_fin_integral = abs(integral_diff_fin_integral - true_value) / true_value;

% Visualizzare i risultati
fprintf('Integrale usando la spline (Metodo dei Trapezi): %.4f\n', integral_spline_trapezoidal);
fprintf('Integrale usando il polinomio (Metodo dei Trapezi): %.4f\n', integral_polynomial_trapezoidal);
fprintf('Integrale usando le differenze finite (Metodo dei Trapezi): %.4f\n', integral_diff_fin_trapezoidal);

fprintf('Integrale usando la spline (Metodo integrale): %.4f\n', integral_spline_integral);
fprintf('Integrale usando il polinomio (Metodo integrale): %.4f\n', integral_polynomial_integral);
fprintf('Integrale usando le differenze finite (Metodo integrale): %.4f\n', integral_diff_fin_integral);

fprintf('\nErrori relativi (Metodo dei Trapezi):\n');
fprintf('Spline: %.6f\n', error_spline_trapezoidal);
fprintf('Polinomio: %.6f\n', error_polynomial_trapezoidal);
fprintf('Differenze finite: %.6f\n', error_diff_fin_trapezoidal);

fprintf('\nErrori relativi (Metodo integrale):\n');
fprintf('Spline: %.6f\n', error_spline_integral);
fprintf('Polinomio: %.6f\n', error_polynomial_integral);
fprintf('Differenze finite: %.6f\n', error_diff_fin_integral);
