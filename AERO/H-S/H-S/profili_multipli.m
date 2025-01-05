clc
close all
clear 

addpath mat_functions


%% Input

U_inf = 1; 
alpha = 0; 
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);


% VALIDAZIONE 2
% Numero di corpi da analizzare 
NCorpi = 2;
CodiceProfilo = cell(NCorpi, 1);
CodiceProfilo{1} = '23012'; % Modificare per aereo cisterna
CodiceProfilo{2} = '6716'; % Modificare per aereo dietro
Chord = [5,2];
NPannelli = [102,102];

LE_X_Position = [0,30];
LE_Y_Position = [0,-10];

flap_xpos = [0,0]; 
flap_ypos = [0,0];
flap_alpha = [0,0];
theta_prof = [7.1193, 7.9873]; % Modificare per impostare gli angoli di incidenza

% % VALIDAZIONE 3
% % Numero di corpi da analizzare (effetto suolo)
% NCorpi = 4;
% CodiceProfilo = cell(NCorpi, 1);
% CodiceProfilo{1} = '0012';
% CodiceProfilo{2} = '0012'; 
% CodiceProfilo{3} = '23012';
% CodiceProfilo{4} = '23012';
% Chord = [1,1,1,1];
% NPannelli = [101,101,101,101];
% 
% LE_X_Position = [0,2,0,2];
% LE_Y_Position = [0,0,2,2];
% 
% flap_xpos = [0,0,0,0]; 
% flap_ypos = [0,0,0,0];
% flap_alpha = [0,0,0,0];
% theta_prof = [10,25,-10,-5]; % rotazione profilo in senso orario


%% Creazione profilo

for i=1:NCorpi
    [x,y]=createProfile(CodiceProfilo{i},NPannelli(i),Chord(i),flap_xpos(i),flap_ypos(i),flap_alpha(i),theta_prof(i));
    Corpi{i}.x=x;
    Corpi{i}.y=y;
end

%% Creazione di una struttura di pannelli

Centro = cell(NCorpi, 1);
Normale = cell(NCorpi, 1);
Tangente = cell(NCorpi, 1);
Estremo_1 = cell(NCorpi, 1);
Estremo_2 = cell(NCorpi, 1);
alpha = cell(NCorpi, 1);
lunghezza = cell(NCorpi, 1);
L2G_TransfMatrix = cell(NCorpi, 1);
G2L_TransfMatrix = cell(NCorpi, 1);


figure;
legend_string = cell(NCorpi, 1);

for i = 1:NCorpi
    [Centro{i}, Normale{i}, Tangente{i}, Estremo_1{i}, Estremo_2{i}, alpha{i}, lunghezza{i}, L2G_TransfMatrix{i}, G2L_TransfMatrix{i}] = CreaStrutturaPannelli(Corpi{i}, LE_X_Position(i), LE_Y_Position(i));
    
    plot(Centro{i}(:, 1), Centro{i}(:, 2), '-', 'LineWidth', 2);
    legend_string{i} = strcat("Corpo ", num2str(i), ' - NACA ', CodiceProfilo{i});
    hold on;
end

axis equal;
title('Coordinate profili');
xlabel('Z [m]', 'Interpreter', 'latex');
ylabel('Y [m]', 'Interpreter', 'latex');
grid on;
grid minor;
legend(legend_string);

%% Inizializzazione matrici e vettori


NCols = sum(NPannelli) + NCorpi;
NRows = NCols;
matrixA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%%

Uv = zeros(2*sum(NPannelli),sum(NPannelli));
Us = zeros(2*sum(NPannelli),sum(NPannelli));

A1 = 0;

for Corpo_i = 1:NCorpi
    for i = 1:NPannelli(Corpo_i)
        A2 = 0;

        Centro_qui = Centro{Corpo_i}(i, :)';
        Tangente_qui = Tangente{Corpo_i}(i, :)'; 
        Normale_qui = Normale{Corpo_i}(i, :)'; 

        A1 = A1+1;
        u1 = 2*A1-1;
        u2 = u1+1;

        for Corpo_j = 1:NCorpi
            for j = 1:NPannelli(Corpo_j)

                Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';             
                Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
                L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

                U_Sorgente = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                U_Vortice = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                A2 = A2 +1;
                matrixA(A1,A2) = dot(U_Sorgente',Normale_qui);
                Uv(u1,A2) = U_Vortice(1,1);
                Uv(u2,A2) = U_Vortice(2,1);
                Us(u1,A2) = U_Sorgente(1,1);
                Us(u2,A2) = U_Sorgente(2,1);
            end
        end
    end
end
%% Completamento matrice
Vec_norm = zeros (2*NCorpi,2);

A4 = 1;
UJS = zeros(2,1);
UJSN = zeros(2,1);
VN1 = zeros(2,1);
VNN = zeros(2,1);
A5 = 0;
A6 = 0;

for Corpo_i = 1:NCorpi
    for i = 1:NPannelli(Corpo_i)
        A5 = A5+1;
        B5 = 0;
        UV = zeros(2,1);
        Normale_qui = Normale{Corpo_i}(i, :)'; 
        % calcolo parte b
        for Corpo_j = 1:NCorpi
            for j = 1:NPannelli(Corpo_j)
                B5 = B5+1;
                UV(1,1) = Uv(2*A5-1,B5);
                UV(2,1) = Uv(2*A5,B5);
                matrixA(A5,sum(NPannelli)+Corpo_j) = matrixA(A5,sum(NPannelli)+Corpo_j)+dot(UV',Normale_qui);
            end
        end       
    end

    % quantità per calcolo c  
    A4N = A4+2*NPannelli(Corpo_i)-2;
    A7 = 0;
    for Corpo_j = 1:NCorpi
        for j = 1:NPannelli(Corpo_j)
            A7 = A7+1;
            UJS = Us(A4:A4+1,A7);
            UJSN = Us(A4N:A4N+1,A7);
            VT1 = Tangente{Corpo_j}(1, :);
            VTN = Tangente{Corpo_j}(end, :);
            matrixA((sum(NPannelli)+Corpo_i),A7) = dot(UJS,VT1)+dot(UJSN,VTN);
        end
    end
    %end
    A4 = A4+2*NPannelli(Corpo_i);
    UR1 = zeros(2,1);
    URN = zeros(2,1);

    % calcolo d
    for i = 1:NPannelli(Corpo_i)
        A8 = 0;
        A6 = A6+1;
        VT1 = Tangente{Corpo_i}(1, :);
        VTN = Tangente{Corpo_i}(end, :);
        for Corpo_j = 1:NCorpi
            for j = 1:NPannelli(Corpo_j)
                A8 = A8+1;
                if i == 1  
                    UR1(1,1) = Uv(A6*2-1,A8);
                    UR1(2,1) = Uv(A6*2,A8);
                    matrixA(sum(NPannelli)+Corpo_i,sum(NPannelli)+Corpo_j) = matrixA(sum(NPannelli)+Corpo_i,sum(NPannelli)+Corpo_j)+dot(UR1',VT1);
                end
                if i == NPannelli(Corpo_i)
                    URN(1,1) = Uv(A6*2-1,A8);
                    URN(2,1) = Uv(A6*2,A8);
                    matrixA(sum(NPannelli)+Corpo_i,sum(NPannelli)+Corpo_j) = matrixA(sum(NPannelli)+Corpo_i,sum(NPannelli)+Corpo_j)+dot(URN',VTN);
                end
           end
        end
    end
end

%% Termine noto
A7 = 0;
for Corpo_i = 1:NCorpi
    for i =1:NPannelli(Corpo_i)
        % BSi
        A7 = A7+1;
        Normale_qui = Normale{Corpo_i}(i, :)'; 
        TermineNoto(A7) = -dot(U_inf,Normale_qui);
    end
    % BV
    T1 = Tangente{Corpo_i}(1, :);
    TN = Tangente{Corpo_i}(end, :);
    TT = T1+TN;
    TermineNoto(sum(NPannelli)+Corpo_i) = -dot(U_inf,TT);
end

%% Risoluzione sistema lineare
Soluzione = linsolve(matrixA,TermineNoto);

A8 = 0;
for Corpo_i = 1:NCorpi
    for i = 1:NPannelli(Corpo_i)
        A8 = A8+1;
        sigma_mia{Corpo_i}(i) = Soluzione(A8);
    end
    gamma_mia(Corpo_i) = Soluzione(sum(NPannelli)+Corpo_i);
end


%% Calcolo del cp e della velocità sui pannelli

U_Pannelli = cell(NCorpi, 1);
Ut_Pannelli = cell(NCorpi, 1);
Un_Pannelli = cell(NCorpi, 1);
Cp = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
    
    U_Pannelli{Corpo_i} = zeros(NPannelli(Corpo_i),2);
    Ut_Pannelli{Corpo_i} = zeros(NPannelli(Corpo_i),1);
    Un_Pannelli{Corpo_i} = zeros(NPannelli(Corpo_i),1);
    
end

for Corpo_i = 1:NCorpi
    for i =1:NPannelli(Corpo_i)

        U_Pannelli{Corpo_i}(i, :) = U_inf'; 
        Centro_qui = Centro{Corpo_i}(i, :)';
        Tangente_qui = Tangente{Corpo_i}(i, :)'; 
        Normale_qui = Normale{Corpo_i}(i, :)'; 
    
        for Corpo_j = 1:NCorpi
            for j = 1:NPannelli(Corpo_j)

                Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';             
                Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
                L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

                U_Sorgente = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                U_Vortice = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                U_Pannelli{Corpo_i}(i, :) = U_Pannelli{Corpo_i}(i, :) + sigma_mia{Corpo_j}(j) .* U_Sorgente' + gamma_mia(Corpo_j) .* U_Vortice';

            end
        end
        
        Ut_Pannelli{Corpo_i}(i) = dot(U_Pannelli{Corpo_i}(i, :)', Tangente_qui);
        Un_Pannelli{Corpo_i}(i) = dot(U_Pannelli{Corpo_i}(i, :)', Normale_qui);
    end
    
    Cp{Corpo_i} = 1-Ut_Pannelli{Corpo_i}.^2/norm(U_inf)^2;
end


Cl = cell(NCorpi, 1);
for Corpo_i = 1:NCorpi
    Cl_qui = 0;
    for i =1:NPannelli(Corpo_i)

        Normale_qui = Normale{Corpo_i}(i, :)';
        lunghezza_qui = lunghezza{Corpo_i}(i);

        Cl_qui = Cl_qui + (-Cp{Corpo_i}(i)*( lunghezza_qui.*dot(Normale_qui, U_inf_normal )));
    end
    
    Cl{Corpo_i} = Cl_qui / Chord(Corpo_i);
end

%% Grafici

figure;
legend_string = cell(2 * NCorpi, 1);

for Corpo_i = 1:NCorpi
    % Plot Cp curve
    plot(Centro{Corpo_i}(:, 1), -Cp{Corpo_i}, '-', 'LineWidth', 2);
    legend_string{Corpo_i} = strcat('$C_P$ NACA ', CodiceProfilo{Corpo_i});
    hold on;

    % Plot airfoil shape
    plot(Corpi{Corpo_i}.x + LE_X_Position(Corpo_i), Corpi{Corpo_i}.y + LE_Y_Position(Corpo_i), '-', 'LineWidth', 2);
    legend_string{NCorpi + Corpo_i} = strcat('Profilo NACA ', CodiceProfilo{Corpo_i});
end
axis equal
grid on;
grid minor;
title('Struttura profilo e $C_P$');
legend(legend_string);


for Corpo_i = 1:NCorpi
    half = floor(NPannelli(Corpo_i) / 2) + 1;
    Cp_dorso = Cp{Corpo_i};
    Cp_ventre = Cp{Corpo_i}(half:end);
    figure();
    plot(Centro{Corpo_i}(1:half, 1) - min(Centro{Corpo_i}(1:half, 1)), -Cp_dorso(1:half), 'b-', 'LineWidth', 2);
    hold on;
    plot(Centro{Corpo_i}(half:end, 1) - min(Centro{Corpo_i}(half:end, 1)), -Cp_ventre, 'r-', 'LineWidth', 2);
    %Plot airfoil shape
    plot(Corpi{Corpo_i}.x, Corpi{Corpo_i}.y, '-', 'LineWidth', 2);
    legend('C_P Ventre', 'C_P Dorso', ['Profilo NACA ', CodiceProfilo{Corpo_i}]);

    axis equal
    grid on;
    grid minor;
    title(['C_P e profilo NACA ', CodiceProfilo{Corpo_i}])
    hold off;
end



for Corpo_i = 1:NCorpi
    half = floor(NPannelli(Corpo_i)/2)+1;
    Cp_dorso = Cp{Corpo_i};
    figure()
    legend_string = cell(NCorpi, 1);
    plot(Centro{Corpo_i}(1:half, 1) - min(Centro{Corpo_i}(1:half, 1)), -Cp{Corpo_i}(1:half), 'b-', 'LineWidth', 2)
    hold on
    plot(Centro{Corpo_i}(half:end, 1) - min(Centro{Corpo_i}(half:end, 1)), -Cp{Corpo_i}(half:end), 'r-', 'LineWidth', 2)
    legend('Ventre','Dorso')
    hold on
    % % INIZIO: AGGIUNTA PER VALIDAZIONE H-S XFOIL
    % if NCorpi == 1
    %     xfoil = xfoil2matlab(theta_prof, CodiceProfilo{1});
    %     plot(xfoil.cp(:,1), -xfoil.cp(:,2), 'ko', 'MarkerSize', 10);
    %     legend('Ventre','Dorso','xfoil')
    %     xlabel('$\frac{x}{c}$', 'Interpreter','latex');
    %     ylabel('$C_P$', 'Interpreter', 'latex');
    % end
    % % FINE: AGGIUNTA PER VALIDAZIONE H-S XFOIL
    hold off

end

hold off
grid on
grid minor
title("Cp Corpo " + num2str(Corpo_i), 'interpreter','latex')



figure;
legend_string = cell(2 * NCorpi, 1); % Adjusted to accommodate both types of plots

for Corpo_i = 1:NCorpi
    % Plot Cp curve
    plot(Centro{Corpo_i}(:, 1) - min(Centro{Corpo_i}(:, 1)), -Cp{Corpo_i}, '-', 'LineWidth', 2);
    legend_string{Corpo_i} = strcat('$C_P$ NACA ', CodiceProfilo{Corpo_i});
    hold on;

    % Plot airfoil shape
    plot(Corpi{Corpo_i}.x, Corpi{Corpo_i}.y, '-', 'LineWidth', 2);
    legend_string{NCorpi + Corpo_i} = strcat('Profilo NACA ', CodiceProfilo{Corpo_i});
end

axis equal;
grid on;
legend(legend_string, 'interpreter', 'latex');
title('Profili sovrapposti', 'interpreter', 'latex');



%% Altri grafici
% true = creazione attiva
% false = creazione non attiva

ifSaveFigures=false;

if ifSaveFigures 

    xMin = 2000;
    xMax = -2000;
    yMin = 2000;
    yMax = -2000;
    
    for Corpo_i = 1:NCorpi
        
        xMin = min(xMin, min(Centro{Corpo_i}(:, 1)));
        xMax = max(xMax, max(Centro{Corpo_i}(:, 1)));
        
        yMin = min(yMin, min(Centro{Corpo_i}(:, 2)));
        yMax = max(yMax, max(Centro{Corpo_i}(:, 2)));
    end
    
    xMin = xMin - 1;
    xMax = xMax + 1;
    yMin = yMin - 1;
    yMax = yMax + 1;
    



    
    Nx = 400;
    Ny = 400;
    
    
    x = linspace(xMin,xMax, Nx);
    y = linspace(yMin, yMax, Ny);

    [X, Y] = meshgrid(x, y);

    isIn = zeros(Nx, Ny);
    
    t = cputime;
    for (i = 1:Nx)
        for (j = 1:Ny)
            for( Corpo_i = 1:NCorpi)
                Boundary = [ Corpi{Corpo_i}.x+LE_X_Position(Corpo_i) Corpi{Corpo_i}.y+LE_Y_Position(Corpo_i)];
                if (inpolygon(X(i, j), Y(i, j), Boundary(:, 1), Boundary(:, 2)))
                    isIn(i, j) = 1;
                end
            end
        end
    end
    cputime - t


    U_Mesh = zeros(Nx, Ny);
    V_Mesh = zeros(Nx, Ny);
    U_Mesh_Mag = zeros(Nx, Ny);
    Cp_Mesh = zeros(Nx, Ny);

    t = cputime;

        for PointIndex_i = 1:Nx
    %         PointIndex_i
            for PointIndex_j = 1:Ny

                if(~isIn(PointIndex_i, PointIndex_j))

                    U = U_inf'; 
                    Centro_qui = [X(PointIndex_i, PointIndex_j); Y(PointIndex_i, PointIndex_j)];

                    for Corpo_j = 1:NCorpi
                        for j = 1:NPannelli(Corpo_j)

                            Estremo_1_qui = Estremo_1{Corpo_j}(j, :)';             
                            Estremo_2_qui = Estremo_2{Corpo_j}(j, :)';
                            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix{Corpo_j}(j, :, :));
                            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix{Corpo_j}(j, :, :));

                            U_Sorgente = ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
                            U_Vortice = ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);

                            U = U + sigma_mia{Corpo_j}(j) .* U_Sorgente' + gamma_mia(Corpo_j) .* U_Vortice';

                        end
                    end

                    U_Mesh(PointIndex_i, PointIndex_j) = U(1);
                    V_Mesh(PointIndex_i, PointIndex_j) = U(2);
                    U_Mesh_Mag(PointIndex_i, PointIndex_j) = norm(U);
                    Cp_Mesh(PointIndex_i, PointIndex_j) = 1 - norm(U) / norm(U_inf);

                else
                    U_Mesh(PointIndex_i, PointIndex_j) = NaN;
                    V_Mesh(PointIndex_i, PointIndex_j) = NaN;
                    U_Mesh_Mag(PointIndex_i, PointIndex_j) = NaN;
                    Cp_Mesh(PointIndex_i, PointIndex_j) = NaN;
                end
            end
        end


    SavingNameStart = "./Figures/test_";


    UMag_fig = figure;
    contourf(X, Y, U_Mesh_Mag, 100,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x+LE_X_Position(Corpo_j), Corpi{Corpo_j}.y+LE_Y_Position(Corpo_j), '-k');
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
%     legend("$U_{MAG}$", 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_UMag.eps');
    saveas(UMag_fig, SavingName);
    
    
    U_fig = figure;
    contourf(X, Y, U_Mesh, 100,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x+LE_X_Position(Corpo_j), Corpi{Corpo_j}.y+LE_Y_Position(Corpo_j), '-k');
    end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_U.eps');
    saveas(U_fig, SavingName);
    
    
    V_fig = figure;
    contourf(X, Y, V_Mesh, 1000,'LineStyle','None');
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x+LE_X_Position(Corpo_j), Corpi{Corpo_j}.y+LE_Y_Position(Corpo_j), '-k');
     end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $V$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_V.eps');
    saveas(V_fig, SavingName);


    Cp_fig = figure;
    contourf(X, Y, Cp_Mesh, 100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar
    hold on
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x+LE_X_Position(Corpo_j), Corpi{Corpo_j}.y+LE_Y_Position(Corpo_j), '-k');
   end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_Cp.eps');
    saveas(Cp_fig, SavingName);
    
    
    Streamlines_fig = figure;
    contourf(X, Y, U_Mesh_Mag, 100,'LineStyle','None');
    colormap(flipud(hot));
    colorbar
    hold on
    streamslice(X, Y, U_Mesh, V_Mesh, 10);
    for Corpo_j = 1:NCorpi
        plot(Corpi{Corpo_j}.x+LE_X_Position(Corpo_j), Corpi{Corpo_j}.y+LE_Y_Position(Corpo_j), '-k');
   end
    axis equal
    xlabel('x', 'interpreter', 'latex')
    ylabel('y', 'interpreter', 'latex')
    title('Contorno di $U_{MAG}$ e linee di corrente', 'interpreter', 'latex')
    SavingName = strcat(SavingNameStart, '_Streamlines.eps');
    saveas(Streamlines_fig, SavingName);
    
end

