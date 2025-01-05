clc
clear all
close all

% Definisci i percorsi di file
bash_script_path = '/home/luca-spazzolino/script/istruzioni';  % Percorso completo del file bash
output_file = '/home/luca-spazzolino/Scaricati/AERO/H-S/H-S/results_naca0012.txt';  % File di output di XFOIL

% Esegui lo script bash
status = system(['bash ', bash_script_path]);

% Verifica se lo script è stato eseguito correttamente
if status ~= 0
    error('Errore nell''esecuzione dello script bash');
end

% Attendi un po' per assicurarsi che XFOIL abbia scritto i risultati
pause(3);  % Pausa per permettere a XFOIL di completare

% Leggi il file di output di XFOIL
fid_output = fopen(output_file, 'r');
if fid_output == -1
    error('Impossibile aprire il file di output di XFOIL');
end

% Cerca i valori di CL e CM nel file di output
cl = NaN;
cm = NaN;

% Nome del file
filename = 'results_naca0012.txt';

% Apertura e lettura del file
fid = fopen(filename, 'r');
if fid == -1
    error('Impossibile aprire il file %s', filename);
end


% Leggere il contenuto del file
data = textscan(fid, '%f %f %f %f %f %f %f %f %f', ...
                'HeaderLines', 12); % Skip delle righe di intestazione
fclose(fid);

% Estrarre i dati
alpha = data{1}; % Angolo di attacco (alpha)
CL = data{2};    % Coefficiente di portanza
CM = data{5};    % Coefficiente di momento

% Mostrare i risultati
disp('Risultati estratti:');
disp(table(alpha, CL, CM));

% Trova valori per alpha = 2.0
alpha_target = 2.0;
idx = find(alpha == alpha_target, 1);
if ~isempty(idx)
    fprintf('Per alpha = %.1f°:\n', alpha_target);
    fprintf('  CL = %.4f\n', CL(idx));
    fprintf('  CM = %.4f\n', CM(idx));
else
    fprintf('Valore di alpha = %.1f° non trovato.\n', alpha_target);
end


