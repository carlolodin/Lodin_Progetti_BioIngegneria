function [num_basi, error] = conta_basi(s)
% conta le basi di un semento di DNA e plotta i risultati in un grafico a barre

% definizione errori e array basi
n = length(s);
error = 0; % se =1 è presente un errore
num_basi = zeros(1,4); % A C G T

% controllo errore lunghezza
if n>1000
    error = 1;
    num_basi = [NaN, NaN, NaN, NaN];
    fprintf('Errore analizzando la stringa %s: la stringa fornita è più lunga di 1000 caratteri\n', inputname(1));
    return;
end

% conta delle basi
for i=1:n
    base = s(i);
    switch base
        case 'A'
            num_basi(1) = num_basi(1) + 1;
        case 'C'
            num_basi(2) = num_basi(2) + 1;
        case 'G'
            num_basi(3) = num_basi(3) + 1;
        case 'T'
            num_basi(4) = num_basi(4) + 1;
        case 'U'
            error = 1;
            num_basi = [NaN, NaN, NaN, NaN];
            fprintf('Errore analizzando la stringa %s: la stringa fornita potrebbe rappresentare RNA\n', inputname(1));
            return;
        otherwise
            error = 1;
            num_basi = [NaN, NaN, NaN, NaN];
            fprintf('Errore analizzando la stringa %s: la stringa fornita non rappresenta DNA\n', inputname(1));
            return;
    end
end



% plot del grafico a barre
X = categorical({'Adenina','Citosina','Guanina','Timina'});
Y = num_basi;
figure();
b = bar(X,Y);
b.FaceColor = 'flat';
b.CData(1,:) = [0.6350 0.0780 0.1840]; % personalizzazione dei colori
b.CData(2,:) = [0.9290 0.6940 0.1250];
b.CData(3,:) = [0.4660 0.6740 0.1880];
b.CData(4,:) = [0.3010 0.7450 0.9330];
title("Conteggio delle basi di DNA di " + inputname(1));

end