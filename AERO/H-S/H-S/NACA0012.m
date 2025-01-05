clc
close all
clear

% 1. Parametri del profilo NACA 0012
n_points = 200;  % Numero di punti da utilizzare per il profilo
x = linspace(0, 1, n_points);  % Punti lungo la corda (0 a 1)

% 2. Formula per calcolare la coordinata y del profilo NACA 0012
y = 0.12 * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x.^2 + 0.2843 * x.^3 - 0.1036 * x.^4);

% 3. Profilo alare completo (dorso e ventre)
x_profile = [x, fliplr(x)];  % Unire il bordo di attacco e di uscita
y_profile = [y, -fliplr(y)];  % Il ventre è il riflesso del dorso

% Visualizzare il profilo alare
figure;
plot(x_profile, y_profile, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('x');
ylabel('y');
title('Profilo Alare NACA 0012');
axis equal;

% 4. Calcolare la linea media
% La linea media è la media tra il dorso e il ventre, che è simmetrico per il NACA 0012
y_midline = (y + (-fliplr(y))) / 2;  % Media tra dorso e ventre

% Visualizzare la linea media
hold on;
plot(x, y_midline, 'r-', 'LineWidth', 1.5);
legend('Profilo Alare', 'Linea Media');
hold off;

% 5. Derivata prima della linea media usando tre metodi
% a) Derivata tramite differenze finite
dx = diff(x);  % Calcolare la differenza nei punti x
dy = diff(y_midline);  % Differenza nei valori y della linea media
deriv_diff_fin = dy ./ dx;

% b) Derivata tramite spline
spline_deriv = spline(x, y_midline);
deriv_spline = @(x) ppval(fnder(spline_deriv), x);

% c) Derivata tramite polinomio (fit polinomiale)
poly_coeffs = polyfit(x, y_midline, 3);
deriv_polynomial = @(x) polyval(polyder(poly_coeffs), x);

% 6. Definire la funzione dell'integrale
L = max(x) - min(x);  % Lunghezza totale del profilo (assumiamo che x sia simmetrico)
integrand_spline = @(s) (1/pi) * deriv_spline(s) ./ sqrt((L^2 / 4) - s.^2);
integrand_polynomial = @(s) (1/pi) * deriv_polynomial(s) ./ sqrt((L^2 / 4) - s.^2);
integrand_diff_fin = @(s) (1/pi) * interp1(x(1:end-1), deriv_diff_fin, s, 'linear', 'extrap') ./ sqrt((L^2 / 4) - s.^2);

% 7. Metodo dei trapezi per il calcolo dell'integrale
function integral = trapezoidal_rule(f, x)
    h = x(2) - x(1);  % Calcolare il passo
    sum_f = f(x(1)) + f(x(end));  % Somma ai bordi
    sum_f = sum_f + 2 * sum(f(x(2:end-1)));  % Somma agli interni
    integral = (h / 2) * sum_f;  % Applicare la formula dei trapezi
end

% 8. Calcolare l'integrale numerico per ciascun metodo
% Metodo dei trapezi
s = linspace(min(x), max(x), 1000);  % Definire i punti per l'integrazione
integral_spline_trapezoidal = trapezoidal_rule(integrand_spline, s);
integral_polynomial_trapezoidal = trapezoidal_rule(integrand_polynomial, s);
integral_diff_fin_trapezoidal = trapezoidal_rule(integrand_diff_fin, s);

% Metodo di integrazione adattiva con 'integral' di MATLAB
integral_spline_integral = integral(integrand_spline, min(x), max(x), 'RelTol', 1e-6);
integral_polynomial_integral = integral(integrand_polynomial, min(x), max(x), 'RelTol', 1e-6);
integral_diff_fin_integral = integral(integrand_diff_fin, min(x), max(x), 'RelTol', 1e-6);

% 9. Calcolare l'errore relativo per i vari metodi
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
