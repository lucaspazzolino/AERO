clc
close all
clear 

addpath mat_functions


%% Input

NACA = {'0006'; '0012'; '23012'; '23045'}; % Specificare i profili da analizzare

alpha_i = 0; % angolo iniziale
alpha_f = 10; % angolo finale
alpha_step = 1; % incremento



%%
v = [alpha_i:alpha_step:alpha_f];

for n = 1:length(NACA)

    for k = 1:length(v)
        theta_prof = v(k);
        
        U_inf = 1; 
        alpha = 0;   
        U_inf_x = U_inf * cos(deg2rad(alpha));
        U_inf_y = U_inf * sin(deg2rad(alpha));
        
        U_inf = [U_inf_x; U_inf_y];
        U_inf_normal = [-U_inf(2); U_inf(1)];
        U_inf_normal = U_inf_normal ./ norm(U_inf_normal);
        
        % VALIDAZIONE 1
        NCorpi = 1;
        CodiceProfilo = cell(NCorpi, 1);
        CodiceProfilo{1} = NACA{n};
        Chord = [1];
        NPannelli = [102];
        
        LE_X_Position = [0];
        LE_Y_Position = [0];
        
        flap_xpos = [0]; 
        flap_ypos = [0];
        flap_alpha = [0];
        
        
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

        legend_string = cell(NCorpi, 1);
        
        for i = 1:NCorpi
            [Centro{i}, Normale{i}, Tangente{i}, Estremo_1{i}, Estremo_2{i}, alpha{i}, lunghezza{i}, L2G_TransfMatrix{i}, G2L_TransfMatrix{i}] = CreaStrutturaPannelli(Corpi{i}, LE_X_Position(i), LE_Y_Position(i));

        end
        
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
       

        xfoil = xfoil2matlab(theta_prof, CodiceProfilo{1});
        
        alpha_vec(k) = theta_prof;
        if NCorpi == 1
            CL_XFOIL(k) = xfoil(1).polar(2);
            CL_HS(k) = Cl{1};
            err_a_cl(k) = abs(CL_HS(k) - CL_XFOIL(k));
            err_rel_cl(k) = err_a_cl(k)/CL_XFOIL(k);
        end

    end

    figure
    plot(alpha_vec, CL_HS, 'LineWidth', 2);
    hold on;
    plot(alpha_vec, CL_XFOIL, 'ko', 'MarkerSize', 10);
    hold off;
    title(['NACA', CodiceProfilo]);
    legend('Hess-Smith', 'xfoil');
    xlabel('$\alpha$ [deg]', 'interpreter', 'latex');
    ylabel('$C_{L}$', 'interpreter', 'latex');
    grid on;
    grid minor;

    figure
    plot(alpha_vec, err_rel_cl, 'r', 'LineWidth', 2);
    title('errore relativo del Cl su alpha');

    err_rel(:, n) = err_rel_cl;

end

err_rel_perc = err_rel .* 100;


figure;
hold on;

for n = 1:length(NACA)
    plot(alpha_vec, err_rel_perc(:, n), 'LineWidth', 2);
    legend_entries{n} = ['NACA ' num2str(NACA{n})];
end

title('Errore relativo percentuale')
xlabel('$\alpha$ [deg]', 'Interpreter', 'latex');
ylabel('$E_{rel_{\%}}$', 'Interpreter', 'latex');
grid on;
grid minor;

legend(legend_entries);
