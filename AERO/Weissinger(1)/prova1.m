clc
clear all 
close

% Parametri dell'ala del Cessna 172
b_cessna = 11.03; % apertura alare in metri
c_cessna = 1.2;   % corda media in metri
S_cessna = 16.2;  % superficie alare in m^2
AR_cessna = b_cessna^2 / S_cessna; % Aspetto dell'ala Cessna

% Parametri dell'ala di un altro aereo (esempio: Piper Cub)
b_piper = 8.5;   % apertura alare in metri
c_piper = 1.1;   % corda media in metri
S_piper = 11.9;  % superficie alare in m^2
AR_piper = b_piper^2 / S_piper; % Aspetto dell'ala Piper

% Intervallo di angoli di attacco
alpha_range = -5:1:15; % angoli di attacco in gradi
alpha = deg2rad(alpha_range); % angoli di attacco in radianti



% Calcolare la distribuzione della portanza per il Cessna 172
[gamma_cessna, cl_alpha_cessna] = weissenger_method(b_cessna, c_cessna, S_cessna, alpha);

% Calcolare la distribuzione della portanza per il Piper Cub
[gamma_piper, cl_alpha_piper] = weissenger_method(b_piper, c_piper, S_piper, alpha);

% Media della circolazione per ciascun angolo di attacco
gamma_cessna_avg = mean(gamma_cessna, 2); % Media lungo i segmenti
gamma_piper_avg = mean(gamma_piper, 2); % Media lungo i segmenti

% Comparazione della distribuzione della circolazione e della pendenza di C_L/alpha
figure;
subplot(2,1,1);
plot(alpha_range, gamma_cessna_avg, 'b', alpha_range, gamma_piper_avg, 'r');
xlabel('Angolo di attacco (°)');
ylabel('Circolazione Media (Gamma)');
legend('Cessna 172', 'Piper Cub');
title('Distribuzione della Circolazione (media)');

subplot(2,1,2);
plot(alpha_range, cl_alpha_cessna*ones(size(alpha_range)), 'b', alpha_range, cl_alpha_piper*ones(size(alpha_range)), 'r');
xlabel('Angolo di attacco (°)');
ylabel('C_L/\alpha');
legend('Cessna 172', 'Piper Cub');
title('Pendenza della Curva C_L/\alpha');


% Calcolare le polari per il Cessna 172
[cl_cessna, cd_cessna] = calculate_polars(gamma_cessna_avg, b_cessna, c_cessna, AR_cessna);

% Calcolare le polari per il Piper Cub
[cl_piper, cd_piper] = calculate_polars(gamma_piper_avg, b_piper, c_piper, AR_piper);

% Tracciare le polari
figure;
plot(cd_cessna, cl_cessna, 'b', cd_piper, cl_piper, 'r');
xlabel('C_D');
ylabel('C_L');
legend('Cessna 172', 'Piper Cub');
title('Polari C_L vs C_D');