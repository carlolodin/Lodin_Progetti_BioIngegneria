%Esercitazione 4
% Sono forniti i dati di Forza (N) e spostamento (mm) di 3 tipologie di campioni a sezione
%rettangolare a base di cheratina e gelatina (3x)
% Lunghezza, larghezza e spessore iniziale sono noti
% Per ciascuna tipologia di campione tracciare i diagrammi sforzo(kPa) deformazione (%)
%dopo aver ‘ripulito’ i dati grezzi
% Determinare visivamente il tratto lineare e calcolare il modulo elastico
% Confronto regressione lineare (stima R2 ed errore massimo) vs coefficiente angolare
%estremi
% Determinare sforzo a snervamento, sforzo massimo, deformazione a massimo sforzo,
%sforzo a rottura, massima elongazione
% Calcolare resilienza e tenacità a rottura (trapz vs formula)
% Riportare i dati su un grafico a barre come media±deviazione standard

data = cell(3, 1); % 1:1, 1:2, 2:1
data{1} = readmatrix("Campioni1_1.xlsx");
data{2} = readmatrix("Campioni1_2.xlsx");
data{3} = readmatrix("Campioni2_1.xlsx");

dimensioni = cell(3, 1);
dimensioni{1} = [11.7, 1.1, 20.3; 9.7, 1.1, 22.5; 10.7, 1, 22.2]; % base1, base2, lunghezza in mm, le righe sono i campioni
dimensioni{2} = [9.7, 1, 22.8; 9.9, 1, 17.6; 10, 1, 18.2]; % base1, base2, lunghezza in mm
dimensioni{3} = [9.3, 1.1, 22.3; 7.6, 1.2, 23.8; 9.2, 1.1, 24.8]; % base1, base2, lunghezza in mm

val = zeros(3, 3, 9); % sigmasnerv, epssnerv, sigmamax, epsmax, sigmarot, epsrot, resilienza, tenacita, E
stat = zeros(3, 2, 9); % righe = tipologie, colonne = media/std, terza dim = caratteristica

young_modulus = zeros(3,3);
sigma_snerv = zeros(3,3);
eps_snerv = zeros(3,3);
sigma_max = zeros(3,3);
eps_max_stress = zeros(3,3);
sigma_break = zeros(3,3);
eps_break = zeros(3,3);
resilienza = zeros(3,3);
tenacita = zeros(3,3);

for tipologia = 1:3
    for campione = 1:3
        % Estrazione dei dati
        forza = data{tipologia}(:, 1+(campione-1)*2);
        forza = forza - 9.6*forza(1)/10;
        spostamento = data{tipologia}(:, 2+(campione-1)*2);
        
        % Pulizia dei dati
        forza(forza < 0) = 0; % Forza negativa non ha senso
        spostamento(spostamento < 0) = 0; % Spostamento negativo non ha senso
        
        % Calcolo della sezione del campione
        sezione = dimensioni{tipologia}(campione, 1) * dimensioni{tipologia}(campione, 2); % mm^2
        
        % Calcolo dello sforzo (kPa)
        sforzo = (forza / sezione) * 1000; % Convertito in kPa
        sforzo = smoothdata(sforzo, 'loess', 300);
        
        % Calcolo della deformazione (%)
        deformazione = (spostamento / dimensioni{tipologia}(campione, 3)) * 100; % Convertito in %
        
        valid_idx = ~isnan(sforzo) & ~isnan(deformazione);
        sforzo = sforzo(valid_idx);
        deformazione = deformazione(valid_idx);
        
        % Regressione lineare per calcolare il modulo elastico
        p = polyfit(deformazione(1:floor(length(deformazione)/4)), sforzo(1:floor(length(deformazione)/4)), 1);
        young_modulus(tipologia, campione) = p(1);
        % Definizione della retta: y = p(1)*(x - 0.2)
        retta = p(1) * (deformazione - 0.2);

        % Trova la differenza tra retta e sforzo
        diff = sforzo - retta;
        
        % Trova dove cambia segno la differenza (intersezione)
        idx = find(diff(1:end-1).*diff(2:end) < 0, 1);
        sigma_snerv(tipologia, campione) = forza(idx);
        eps_snerv(tipologia, campione) = deformazione(idx);


        % Visualizza il punto di intersezione
        figure;
        hold on;
        grid on;
        plot(deformazione(idx), sforzo(idx), 'ro', 'MarkerFaceColor', 'r');
        plot(deformazione, retta);
        text(deformazione(idx) + 0.1, sforzo(idx)+ .1, sprintf('(%.2f, %.2f)', deformazione(idx), sforzo(idx)), 'FontSize', 12);

        %Trova massimo sforzo
        [max_val, idx_max] = max(sforzo);
        sigma_max(tipologia, campione) = sforzo(idx_max);
        eps_max_stress(tipologia, campione) = deformazione(idx_max);
        plot(eps_max_stress(tipologia, campione), sigma_max(tipologia, campione), 'ro', 'MarkerFaceColor', 'r');
        text(eps_max_stress(tipologia, campione) + 0.1, sigma_max(tipologia, campione)+ .1, ...
        sprintf('(%.2f, %.2f)', eps_max_stress(tipologia, campione), sigma_max(tipologia, campione)), 'FontSize', 12);

        %Trova sforzo e deformazione a rottura
        sigma_break(tipologia, campione) = sforzo(end);
        eps_break(tipologia, campione) = deformazione(end);
        plot(deformazione(end), sforzo(end), 'ro', 'MarkerFaceColor', 'r');
        text(deformazione(end) + 0.1, sforzo(end)+ .1, ...
        sprintf('(%.2f, %.2f)', deformazione(end), sforzo(end)), 'FontSize', 12);

        % Tracciamento del diagramma sforzo-deformazione
        plot(deformazione, sforzo);
        xlabel('Deformazione (%)');
        ylabel('Sforzo (kPa)');
        title(['Campione ' num2str(campione) ' - Tipologia ' num2str(tipologia)]);

        %Calcola resilienza e tenacità
        resilienza(tipologia, campione) = trapz(deformazione(1:idx), sforzo(1:idx));
        tenacita(tipologia, campione) = trapz(deformazione, sforzo);

        text(3, 1, ['Resilienza: ', num2str(resilienza(tipologia, campione)), ...
            'J/m^3  Tenacità: ', num2str(tenacita(tipologia, campione)),'J/m^3'], ...
            'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    end
end

media = [ ...
    mean(young_modulus, 2), ...
    mean(sigma_snerv, 2), ...
    mean(eps_snerv, 2), ...
    mean(sigma_max, 2), ...
    mean(eps_max_stress, 2), ...
    mean(sigma_break, 2), ...
    mean(eps_break, 2), ...
    mean(resilienza, 2), ...
    mean(tenacita, 2) ...
];

devstd = [ ...
    std(young_modulus, 0, 2), ...
    std(sigma_snerv, 0, 2), ...
    std(eps_snerv, 0, 2), ...
    std(sigma_max, 0, 2), ...
    std(eps_max_stress, 0, 2), ...
    std(sigma_break, 0, 2), ...
    std(eps_break, 0, 2), ...
    std(resilienza, 0, 2), ...
    std(tenacita, 0, 2) ...
];

nomi = { ...
    'Modulo Young', ...
    'Sforzo snerv.', ...
    'Deform. snerv.', ...
    'Sforzo max', ...
    'Deform. max', ...
    'Sforzo rottura', ...
    'Deform. rottura', ...
    'Resilienza', ...
    'Tenacità' ...
};

figure;
for i = 1:9
    subplot(3,3,i)
    bar(1:3, media(:,i), 'FaceColor', [0.2 0.6 0.8]);
    hold on
    errorbar(1:3, media(:,i), devstd(:,i), 'k.', 'LineWidth', 1.5, 'CapSize', 12);
    hold off
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Tipo 1', 'Tipo 2', 'Tipo 3'});
    title(nomi{i});
    grid on
end
sgtitle('Media e deviazione standard delle proprietà meccaniche')
