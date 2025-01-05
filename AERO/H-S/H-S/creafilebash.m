clc
clear
close

% Percorso del file di input e di output
input_file_path = '/home/luca-spazzolino/codici/naca0012_input.txt';
output_file_path = '/home/luca-spazzolino/codici/xfoil_output.txt';

% Numero di Reynolds e angolo di attacco
Re = 1e6;
alpha = 2.0;

% Crea il file di input per XFOIL
fileID = fopen(input_file_path, 'w');
fprintf(fileID, 'NACA 0012\n');                  % Specifica il profilo NACA
fprintf(fileID, 'OPER\n');                       % Entra in modalità operativa
fprintf(fileID, 'Visc %.0f\n', Re);              % Imposta il numero di Reynolds
fprintf(fileID, 'Alfa %.2f\n', alpha);           % Imposta l'angolo di attacco
fprintf(fileID, 'PACC\n');                       % Attiva la scrittura dei risultati
fprintf(fileID, '%s\n\n', output_file_path);     % Specifica il file di output
fclose(fileID);
