function calcola_lipidi(peso_min, peso_max, num_pazienti, minL, maxL)

    %% Costanti del problema
    
    dens_lipidi = 4e18;             % lipidi/m^2
    dens_cellule = 5.314e14;        % cell/m^3
    PM_lipidi = 500;                % g/mol
    diam_cellula = 20e-6;           % m
    dens_tessuto = 1000;            % kg/m^3
    N_Avogadro = 6.022e23;

    %% Inizializzazione e calcolo
    
    pesi = randi([peso_min peso_max], num_pazienti, 1);    %Inizializza i pesi dei pazienti in un vettore colonna
    volumi = pesi ./ dens_tessuto;                      % Calcola il volume del tessuto per paziente (m^3)
    num_cellule = volumi * dens_cellule;                 % Calcola il numero di cellule per paziente
    
    superficie_cellula = 4 * pi * (diam_cellula/2)^2;    % Calcola la superficie di una cellula (m^2)

    lipidi_tot = dens_lipidi * superficie_cellula .* num_cellule;  % numero di lipidi per paziente
    lipidi_grammi = (lipidi_tot / N_Avogadro) * PM_lipidi;            % grammi di lipidi per paziente

    %% Plot istogramma

    edges = minL:10:maxL;  % intervalli con passo 10
    histogram(lipidi_grammi, edges)
    xlabel('Grammi di lipidi')
    ylabel('Numero di pazienti')
    title('Distribuzione grammi di lipidi per paziente')
    grid on;

    % Verifica se alcuni valori non sono inclusi nel grafico
    if min(lipidi_grammi) < minL || max(lipidi_grammi) > maxL
        warning('Alcuni valori di lipidi non sono riportati sul grafico. Intervallo calcolato: [%.2f, %.2f] g', valoreMin, valoreMax);
    end

end
