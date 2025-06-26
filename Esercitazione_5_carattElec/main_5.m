%%Preparazione dati
dati = readmatrix("dati_caratterizzazione elettrica.xlsx", "Range", "A3:Z203");

tf = cell(4,3); % FdT, righe = tipologia, colonne = campione

moduli = cell(4,3); 
fasi = cell(4,3);

freq = dati(:,1);
puls = dati(:,2);

for i = 1:4
    for j = 1:3
        moduli{i,j} = dati(:, 3 + (i-1)*6+ (j-1)*2);
        fasi{i,j} = dati(:, 4 + (i-1)*6+ (j-1)*2);
    end 
end

%% FdT e fitting

fase_fit = cell(4,3); % Modelli fittat. Una riga per tipologia, una colonna per campione
mod_fit = cell(4,3);

% determinazione delle funzioni di trasferimento
for i = 1:4  
    for j = 1:3
        imm = sqrt(-1);
        
        d = frd(moduli{i, j}.*exp(deg2rad(fasi{i, j})*imm), puls); 
        tf{i, j} = tfest(d, 2, 1);  % 2 poli e 1 zero

        resp = squeeze(freqresp(tf{i, j}, puls));
        mod_fit{i, j} = (abs(resp));
        fase_fit{i, j} = rad2deg((angle(resp)));

    end
end

%% punto 1

titoli.tipo = {'P-G 1:0.1', 'P-G 1:0.2', 'P-G 1:0.3', 'P-G 1:0.4'};
titoli.cam = {'campione A', 'campione B', 'campione C'};

for i = 1:4
    figure();
    sgtitle(titoli.tipo{i})
    for j = 1:3   
        subplot(2, 3, j)
        semilogx(puls, moduli{i,j}/10^3)
        hold on
        semilogx(puls, mod_fit{i, j}/10^3, ':', 'LineWidth', 1)
        xlim([0, 10^8])
        title(['Modulo ',' ', titoli.cam{j}]);
        legend('dati sperimentali', 'fitting')
        xlabel('\omega [rad/s]')
        ylabel('[k\Omega]')
        hold off

        subplot(2, 3, j+3)
        semilogx(puls, fasi{i,j});
        hold on
        semilogx(puls, fase_fit{i, j}, ':', 'LineWidth', 1);
        xlim([0, 10^8])
        title(['Fase', ' ', titoli.cam{j}]);
        legend('dati sperimentali', 'fitting')
        xlabel('\omega [rad/s]')
        ylabel('[deg]')
        hold off
    end

end

%% punto 2
% Per ciascuna tipologia calcoli il valor medio di modulo e fase 
% punto per punto e ne tracci il grafico

moduli_medi = zeros(4, length(dati));
fasi_medie = zeros(4, length(dati));
moduli_medi_fit = zeros(4, length(dati));
fasi_medie_fit = zeros(4, length(dati));

for i = 1:4
    moduli_medi(i,:) = mean(cat(2, moduli{i, :}), 2);
    fasi_medie(i,:) = mean(cat(2, fasi{i, :}), 2);
    moduli_medi_fit(i,:) = mean(cat(2, mod_fit{i, :}), 2);
    fasi_medie_fit(i,:) = mean(cat(2, fase_fit{i, :}), 2);
end

for i = 1:4
    figure()
    sgtitle(['Media dei campioni per', ' ', titoli.tipo{i}])

    subplot(2, 1, 1)
    semilogx(puls, moduli_medi(i, :)/10^3);
    xlabel('\omega [rad/s]')
    ylabel('[k\Omega]')
    hold on
    semilogx(puls, moduli_medi_fit(i, :)/10^3, ':', 'LineWidth', 1)
    hold off
    legend('dati sperimentali', 'fitting')
    title('Modulo medio')

    subplot(2, 1, 2)
    semilogx(puls, fasi_medie(i, :));
    xlabel('\omega [rad/s]')
    ylabel('[deg]')
    hold on
    semilogx(puls, fasi_medie_fit(i, :), ':', 'LineWidth', 1)
    hold off
    legend('dati sperimentali', 'fitting')
    title('Fase media')
end

%% punto 3
% Riporti su un grafico il valor medio dei moduli medi per ciascuna 
% tipologia di campione.

figure;
semilogx(puls, mean(moduli_medi, 1), 'LineWidth', 1.5)
title('Moduli medi')
hold on
for i = 1:4
    semilogx(puls, moduli_medi_fit(i,:), 'LineWidth', 1.5)
end
hold off
xlabel('\omega [rad/s]')
ylabel('[k\Omega]')
legend('Media delle medie dei dati originali', 'P-G 1:0.1', 'P-G 1:0.2', ...
    'P-G 1:0.3','P-G 1:0.4')

%% punto 4
% Per ciascun campione calcoli l’errore quadratico medio tra i moduli di
% fitting e dati sperimentali

MSE = zeros(4, 3); % righe = tipologie, colonne = campioni 

for i = 1:4
    for j = 1:3
    mse = sum((mod_fit{i, j} - moduli{i, j}).^2)/length(moduli{i, j});
    MSE(i,j) = mse;
    end
end

%% punto 5

circuito1 = imread('circuito1.jpg');
circuito2 = imread('circuito2.jpg');

% Mostra l'immagine
figure
imshow(circuito1);
title('modello circuitale 1')
figure
imshow(circuito2);
title('modello circuitale 2')

%% punto 6

% primo circuito
R11 = zeros(4,3);
R21 = zeros(4,3);
C11 = zeros(4,3);
C21 = zeros(4,3);

for i = 1:4
    for j = 1:3
        den2 = tf{i,j}.Denominator(2);
        den3 = tf{i,j}.Denominator(3);
        num1 = tf{i,j}.Numerator(1);
        num2 = tf{i,j}.Numerator(2);
        
        R11(i, j) = num2/den3;
        C11(i, j) = 1/num1;
        C21(i, j) = den2/num2 - 1/num1 - den3*num1/((num2)^2);
        R21(i, j) = num1/(num2*C21(i, j));
    end
end

R11_m = abs(mean(R11, 2));
R21_m = abs(mean(R21, 2));
C11_m = abs(mean(C11, 2));
C21_m = abs(mean(C21, 2));

R11_std = std(R11, 1, 2);
R21_std = std(R21, 1, 2);
C11_std = std(C11, 1, 2);
C21_std = std(C21, 1, 2);

% secondo circuito
R12 = zeros(4,3);
R22 = zeros(4,3);
C12 = zeros(4,3);
C22 = zeros(4,3);

for i = 1:4
    for j = 1:3
        den2 = tf{i,j}.Denominator(2);
        den3 = tf{i,j}.Denominator(3);
        num1 = tf{i,j}.Numerator(1);
        num2 = tf{i,j}.Numerator(2);

        C12(i, j) = (den2*num2-num1*den3)/((num2)^2);
        R12(i, j) = (den2/den3 - num1/num2)/C12(i, j);
        C22(i, j) = C12(i, j)/(C12(i, j)*num1 - 1);
        R22(i, j) = (C12(i, j)*num1 - 1)/(C12(i, j)^2*num2);
    end
end

R12_m = mean(R12, 2);
R22_m = mean(R22, 2);
C12_m = mean(C12, 2);
C22_m = mean(C22, 2);

R12_std = std(R12, 1, 2);
R22_std = std(R22, 1, 2);
C12_std = std(C12, 1, 2);
C22_std = std(C22, 1, 2);

%% punto 7

figure
sgtitle('parametri circuitali medi al variare della concentrazione')
x = [0.1, 0.2, 0.3, 0.4];
subplot(2,2,1)
errorbar(x, R11_m/1000, R11_std/1000, '--o')
xlim([0, 0.5])
hold on
errorbar(x, R21_m/1000, R21_std/1000, '--o')
hold off 
legend('R1', 'R2')
xticks(x)
xticklabels(titoli.tipo)
title('resistenze primo circuito')
xlabel('concentrazione di gelatina')
ylabel('[k\Omega]')

subplot(2,2,2)
errorbar(x, C11_m*10^9, C11_std*10^9, '--o')
xlim([0, 0.5])
hold on
errorbar(x, C21_m*10^9, C21_std*10^9, '--o')
hold off 
legend('C1', 'C2')
xticks(x)
xticklabels(titoli.tipo)
title('capacità primo circuito')
xlabel('concentrazione di gelatina')
ylabel('[nF]')

subplot(2,2,3)
errorbar(x, R12_m/1000, R12_std/1000, '--o')
xlim([0, 0.5])
hold on
errorbar(x, R22_m/1000, R22_std/1000,'--o')
hold off 
legend('R1', 'R2')
xticks(x)
xticklabels(titoli.tipo)
title('resistenze secondo circuito')
xlabel('concentrazione di gelatina')
ylabel('[k\Omega]')

subplot(2,2,4)
errorbar(x, C12_m*10^9, C12_std*10^9,'--o')
xlim([0, 0.5])
hold on
errorbar(x, C22_m*10^9, C22_std*10^9,'--o')
hold off 
legend('C1', 'C2')
xticks(x)
xticklabels(titoli.tipo)
title('capacità secondo circuito')
xlabel('concentrazione di gelatina')
ylabel('[nF]')



