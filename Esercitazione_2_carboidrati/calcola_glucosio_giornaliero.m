function calcola_glucosio_giornaliero(kcal_min, kcal_max, num_soggetti)

    %% Costanti del problema
    
    pm_glucosio = 180;              % g/mol
    energia_gluc_ossidasi = 3000;   % KJ/mol
    fraz_energ_carb = .55;        
    fraz_colazione = .3;           
    fraz_pranzo = .5;            
    fraz_cena = .2;
    
    %% Verifica degli input
    if ~isscalar(kcal_min) || ~isscalar(kcal_max) || ~isscalar(num_soggetti) || mod(num_soggetti,1) ~= 0 || ~isnumeric(kcal_min) || ~isnumeric(kcal_max) || ~isnumeric(num_soggetti)
        error('kcal_min e kcal_max devono essere scalari compresi tra 2000 e 3000, num_soggetti deve essere un intero >= 1.');
    end
    if kcal_min < 2000 || kcal_min > 3000
        error('I valori di kcal_min devono essere compresi tra 2000 e 3000.');
    end
    if kcal_max < 2000 || kcal_max > 3000
        error('I valori di kcal_max devono essere compresi tra 2000 e 3000.');
    end
    if num_soggetti <= 0
        error('Il numero di soggetti deve essere maggiore di zero.');
    end
    
    %% Inizializzazione e calcolo
    grammi_glucosio = zeros(num_soggetti, 3); % Inizializza la matrice dei grammi di glucosio
    
    %Inizializza le energie giornaliere da carboidrati in un vettore colonna
    energia_giorn_carb = randi([kcal_min kcal_max], num_soggetti, 1)*fraz_energ_carb*4.184;  % Energia in KJ da carboidrati (kcal * 4.184)
    
    grammi_glucosio_tot = energia_giorn_carb*pm_glucosio/energia_gluc_ossidasi;  % Calcola i grammi di glucosio giornalieri
    grammi_glucosio(:, 1) = grammi_glucosio_tot * fraz_colazione;  % Colazione
    grammi_glucosio(:, 2) = grammi_glucosio_tot * fraz_pranzo;     % Pranzo
    grammi_glucosio(:, 3) = grammi_glucosio_tot *  fraz_cena;      % Cena
    
    %% Plot grafici
    figure;

    subplot(3,1,1);
    plot(1:num_soggetti, grammi_glucosio(:,1), 'o-b');
    xlabel('Numero soggetto');
    ylabel('Glucosio Colazione (g)');
    title('Glucosio assunto a Colazione');

    subplot(3,1,2);
    plot(1:num_soggetti, grammi_glucosio(:,2), 's-r');
    xlabel('Numero soggetto');
    ylabel('Glucosio Pranzo (g)');
    title('Glucosio assunto a Pranzo');

    subplot(3,1,3);
    plot(1:num_soggetti, grammi_glucosio(:,3), 'd-g');
    xlabel('Numero soggetto');
    ylabel('Glucosio Cena (g)');
    title('Glucosio assunto a Cena');

    % Calcola media e range per ciascun pasto
    media_glucosio = mean(grammi_glucosio, 1);
    min_glucosio = min(grammi_glucosio, [], 1);
    max_glucosio = max(grammi_glucosio, [], 1);

    figure;
    bar(media_glucosio);
    hold on;
    errorbar(1:3, media_glucosio, media_glucosio - min_glucosio, max_glucosio - media_glucosio, 'k.', 'CapSize', 15);
    hold off;
    set(gca, 'XTickLabel', {'Colazione', 'Pranzo', 'Cena'});
    ylabel('Glucosio medio (g)');
    title('Media e range di glucosio assunto per pasto');
    legend('Media', 'Range');
end
