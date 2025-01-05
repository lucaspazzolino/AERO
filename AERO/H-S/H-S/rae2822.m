clc
close all
clear 

% Script MATLAB per visualizzare il profilo alare unendo dorso e ventre

% Nomi dei file contenenti i dati
dorso_file = 'dorso.dat';
ventre_file = 'ventre.dat';

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