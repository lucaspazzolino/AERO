clc
clear
close all

% Definisci il percorso di XFOIL e la directory di lavoro
xfoil_path = '/usr/bin/xfoil';  % Modifica con il percorso corretto di XFOIL
input_file = 'clcm';   % File di input per XFOIL (modifica se necessario)
output_file = 'cl_cm_output.txt'; % Nome del file di output

% Crea il file di input XFOIL (se non esiste già)
fid_input = fopen(input_file, 'w');
fprintf(fid_input, 'NACA 0012\n'); % Profilo NACA 0012
fprintf(fid_input, 'Visc 500000.0\n'); % Numero di Reynolds
fprintf(fid_input, 'Alfa 2.0\n'); % Angolo di attacco di 2 gradi
fprintf(fid_input, 'Pacc\n'); % Output di un file
fprintf(fid_input, 'File %s\n', output_file); % Nome del file di output
fprintf(fid_input, 'Quit\n'); % Comando per uscire
fclose(fid_input);

% Esegui XFOIL
system(['cd ' fileparts(xfoil_path) ' && ./xfoil < ' input_file]);

% Attendi che XFOIL finisca e che il file di output venga scritto
pause(2);

% Leggi il file di output di XFOIL
fid_output = fopen(output_file, 'r');
if fid_output == -1
    error('Impossibile aprire il file di output di XFOIL');
end

% Cerca i valori di Cl e Cm nel file di output
cl = NaN;
cm = NaN;

while ~feof(fid_output)
    line = fgetl(fid_output);
    
    % Cerca i valori di Cl e Cm nelle righe corrispondenti
    if contains(line, 'CL')
        cl_values = sscanf(line, ' CL = %f');
        if ~isempty(cl_values)
            cl = cl_values(1);
        end
    end
    
    if contains(line, 'CM')
        cm_values = sscanf(line, ' CM = %f');
        if ~isempty(cm_values)
            cm = cm_values(1);
        end
    end
end

fclose(fid_output);

% Mostra i risultati
disp(['Cl: ', num2str(cl)]);
disp(['Cm: ', num2str(cm)]);

% Controlla se i valori sono NaN
if isnan(cl) || isnan(cm)
    error('I valori di Cl o Cm non sono stati trovati nel file di output di XFOIL.');
end
