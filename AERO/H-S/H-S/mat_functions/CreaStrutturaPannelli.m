function [Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo_Input,X_LE,Y_LE)

NPannelli = length(Corpo_Input.x)-1;

Centro = zeros(NPannelli, 2);
Normale = zeros(NPannelli, 2);
Tangente = zeros(NPannelli, 2);
Estremo_1 = zeros(NPannelli, 2);
Estremo_2 = zeros(NPannelli, 2);
alpha = zeros(NPannelli, 1);
L2G_TransfMatrix = zeros(NPannelli, 2, 2);
G2L_TransfMatrix = zeros(NPannelli, 2, 2);



for i = 1:NPannelli
    
    Centro(i, 1) = X_LE+1/2*(Corpo_Input.x(i)+Corpo_Input.x(i+1));
    Centro(i, 2) = Y_LE+1/2*(Corpo_Input.y(i)+Corpo_Input.y(i+1));
    
    Estremo_1(i, 1) = X_LE+Corpo_Input.x(i);
    Estremo_1(i, 2) = Y_LE+Corpo_Input.y(i);
    
    Estremo_2(i, 1) = X_LE+Corpo_Input.x(i+1);
    Estremo_2(i, 2) = Y_LE+Corpo_Input.y(i+1);
    
    dy = Estremo_2(i,2)-Estremo_1(i,2);
    dx = Estremo_2(i,1)-Estremo_1(i,1);
    
    % compute angle and matrices
    angle = atan2(dy, dx);
    
    L2G_TransfMatrix(i, :, :) = [cos(angle) ,  -sin(angle);
                                 sin(angle),  cos(angle)];
                             
    G2L_TransfMatrix(i, :, :) = [cos(angle) ,  sin(angle);
                                 -sin(angle),  cos(angle)];
                             
    Normale(i, 1) = -sin(angle);
    Normale(i, 2) = cos(angle);
    
    Tangente(i, 1) = cos(angle);
    Tangente(i, 2) = sin(angle);
    
    lunghezza(i) = norm(Estremo_2(i, :) - Estremo_1(i, :));
    
end


