% Aggiungere l'effetto dei piani di coda (approssimazione semplificata)
function [cl_tail, cd_tail] = tail_effects(S_tail, S, cl, cd)
    % Effetto dei piani di coda sull'aerodinamica (approccio semplificato)
    AR_tail = 6; % Aspetto del piano di coda (approssimato)
    cl_tail = cl * (1 - S_tail / S); % Effetto sulla portanza
    cd_tail = cd * (1 + S_tail / S); % Effetto sulla resistenza
end