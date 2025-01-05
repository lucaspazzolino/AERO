clc
clear all
close all 

% Definisci il percorso dello script Bash
bash_script_path = '/home/luca-spazzolino/Scaricati/AERO/H-S/H-S/clcm';

% Esegui lo script Bash
[status, cmdout] = system(bash_script_path);

% Verifica se lo script è stato eseguito correttamente
if status == 0
    disp('XFOIL eseguito correttamente');
    
    % Leggi il file di output che contiene C_L e C_M
    output_file = 'cl_cm_output.txt';
    fileID = fopen(output_file, 'r');
    results = textscan(fileID, '%s %f', 'Delimiter', ':');
    fclose(fileID);
    
    % Estrai i risultati
    CL = results{2}(1);
    CM = results{2}(2);
    
    % Visualizza i risultati
    fprintf('Coefficienti ottenuti da XFOIL:\n');
    fprintf('C_L = %.4f\n', CL);
    fprintf('C_M = %.4f\n', CM);
else
    disp('Errore durante l esecuzione dello script Bash');
    disp(cmdout);
end
