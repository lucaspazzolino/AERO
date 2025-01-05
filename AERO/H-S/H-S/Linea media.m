

 


% Carica i dati del profilo RAE2822
% Assumiamo che il file 'rae2822_coordinates.csv' abbia tre colonne:
% x, y_sup (estradosso), y_inf (intradosso)
data = readmatrix('rae2822_coordinates.csv');

% Estrai le colonne dal file
x = data(:, 1);        % Coordinate x
y_sup = data(:, 2);    % Coordinate y dell'estradosso
y_inf = data(:, 3);    % Coordinate y dell'intradosso

% Calcolo della linea media
y_mean = (y_sup + y_inf) / 2;

% Plot del profilo e della linea media
figure;
hold on;
plot(x, y_sup, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Estradosso');
plot(x, y_inf, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Intradosso');
plot(x, y_mean, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Linea Media');
hold off;

% Aggiunta delle etichette e della leggenda
title('Linea Media del Profilo RAE2822');
xlabel('Corda (x)');
ylabel('Altezza (y)');
legend('Location', 'best');
grid on;
axis equal;

% Carica i dati del profilo RAE2822
% Assumiamo che il file 'rae2822_coordinates.csv' abbia tre colonne:
% x, y_sup (estradosso), y_inf (intradosso)
data = readmatrix('rae2822_coordinates.csv');

% Estrai le colonne dal file
x = data(:, 1);        % Coordinate x
y_sup = data(:, 2);    % Coordinate y dell'estradosso
y_inf = data(:, 3);    % Coordinate y dell'intradosso

% Calcolo della linea media
y_mean = (y_sup + y_inf) / 2;

% Plot del profilo e della linea media
figure;
hold on;
plot(x, y_sup, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Estradosso');
plot(x, y_inf, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Intradosso');
plot(x, y_mean, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Linea Media');
hold off;

% Aggiunta delle etichette e della leggenda
title('Linea Media del Profilo RAE2822');
xlabel('Corda (x)');
ylabel('Altezza (y)');
legend('Location', 'best');
grid on;
axis equal;
