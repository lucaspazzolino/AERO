function [xfoil] = xfoil2matlab(AoA, Codiceprofilo)

file_to_open = 'data/Cp_import.txt';
NACA = ["NACA", Codiceprofilo];

n_angles = length(AoA);

xfoil(n_angles) = struct();

for i = 1:n_angles

    polar_file = ['data/polar_', num2str(i), '.txt'];
    cp_file = ['data/cp_', num2str(i), '.txt'];

    fileID = fopen(file_to_open,'w');

    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', NACA);
    fprintf(fileID, '%s\n', 'PPAR');
    fprintf(fileID, '%s\n', 'N 103');
    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', 'OPER');
    fprintf(fileID, '%s\n', 'PACC');
    fprintf(fileID, '%s\n', polar_file);
    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', ['ALFA ', num2str(AoA(i))]);
    fprintf(fileID, '%s\n', ['CPWR ', cp_file]);
    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', ' ');
    fprintf(fileID, '%s\n', 'QUIT');

    fclose(fileID);

    system('xfoil < data/Cp_import.txt');

    data = importdata(cp_file, ' ', 1);
    xfoil(i).cp = data.data;
    data = importdata(polar_file, ' ', 12);
    xfoil(i).polar = data.data;

    system('rm data/*');

end
