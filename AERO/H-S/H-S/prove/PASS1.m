clc 
clear 
close all 

% Definizione della lunghezza della corda (normalizzata)
L = 1; % La lunghezza è normalizzata

% Numero di punti per la discretizzazione
N = 1000; 

% Generazione delle coordinate x (corda normalizzata)
x = linspace(0, L, N)';

% Parametri caratteristici del profilo RAE2822 (modifica se necessario)
t = 0.12; % Spessore massimo relativo
m = 0.02; % Curvatura massima relativa
p = 0.4;  % Posizione della curvatura massima relativa

% Calcolo delle coordinate della camber line (linea di curvatura)
yc = zeros(size(x));
dyc_dx = zeros(size(x));

for i = 1:N
    if x(i) < p * L
        yc(i) = m * (x(i) / (p * L)^2) * (2 * p * L - x(i));
        dyc_dx(i) = 2 * m * (p * L - x(i)) / (p * L)^2;
    else
        yc(i) = m * ((L - x(i)) / ((1 - p) * L)^2) * (1 + x(i) - 2 * p * L);
        dyc_dx(i) = 2 * m * (p * L - x(i)) / ((1 - p) * L)^2;
    end
end

% Calcolo dello spessore simmetrico
yt = 5 * t * ( ...
      0.2969 * sqrt(x / L) ...
    - 0.1260 * (x / L) ...
    - 0.3516 * (x / L).^2 ...
    + 0.2843 * (x / L).^3 ...
    - 0.1015 * (x / L).^4);

% Calcolo delle coordinate dell'estradosso e dell'intradosso
theta = atan(dyc_dx);
y_sup = yc + yt .* cos(theta);
y_inf = yc - yt .* cos(theta);

% Salvataggio dei dati in un file CSV
data = [x, y_sup, y_inf];
writematrix(data, 'rae2822_coordinates.csv');
disp('Coordinate del profilo RAE2822 salvate in "rae2822_coordinates.csv".');

% Visualizzazione del profilo
figure;
hold on;
plot(x, y_sup, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Estradosso');
plot(x, y_inf, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Intradosso');
title('Profilo RAE2822');
xlabel('Corda (x)');
ylabel('Altezza (y)');
axis equal;
grid on;
legend('Location', 'best');
hold off;
