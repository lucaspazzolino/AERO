clc
clear 
close all

% Impostazioni iniziali
% Specifica il percorso di XFoil
xfoil_path = '/usr/local/XFoil/xfoil';
% Profilo alare RAE2822 (assicurati di avere il file RAE2822.dat nel percorso giusto)
airfoil_file = 'rae2822_coordinates.csv';  % Nome del file con le coordinate del profilo

% Angolo di attacco
alpha = 2;  % Angolo di attacco in gradi

% Impostazioni di flusso
Re = 6e6;  % Numero di Reynolds (esempio)
mach = 0.1;  % Numero di Mach (esempio, per flusso subsonico)

% Crea il comando per eseguire XFoil
cmd = [xfoil_path, ' < ', airfoil_file, ' <<EOF\n'];

% Comando XFoil
cmd = [cmd, 'PLOP\n', 'G\n', 'OPER\n', 'VISC\n', num2str(Re), '\n', 'MACH\n', num2str(mach), '\n', 'ALFA\n', num2str(alpha), '\n', 'EOF\n'];

% Esegui XFoil utilizzando system per invocare il comando
[status, cmdout] = system(cmd);

% Analizza i risultati di XFoil
% Salva i risultati di Cp, CL, e CM in un file temporaneo
output_file = 'xfoil_output.txt';
if status == 0
    % Leggi il file di output di XFoil
    output_data = load(output_file);
    % Estrai e visualizza i dati di Cp, CL, e CM
    Cp = output_data(:, 2);  % Colonna per Cp
    x = output_data(:, 1);    % Posizione sulla corda
    CL = output_data(end, 3);  % CL finale
    CM = output_data(end, 4);  % CM finale
    
    % Plot dei risultati
    figure;
    subplot(1, 2, 1);
    plot(x, Cp);
    title('Distribuzione Cp');
    xlabel('Posizione sulla corda');
    ylabel('Cp');
    
    subplot(1, 2, 2);
    bar([CL, CM]);
    title('CL e CM');
    ylabel('Coefficiente');
    set(gca, 'xticklabel', {'CL', 'CM'});
else
    error('Errore nell''esecuzione di XFoil');
end
