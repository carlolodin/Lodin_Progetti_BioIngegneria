function calcola_lipidi(pesoMin, pesoMax, numPazienti, minL, maxL)

    %% Costanti del problema
    
    dens_lipidi = 4e18;             % lipidi/m^2
    dens_cellule = 5.314e14;        % cell/m^3
    PM_lipidi = 500;                % g/mol
    diam_cellula = 20e-6;           % m
    dens_tessuto = 1000;            % kg/m^3
    N_Avogadro = 6.022e23;

    %% Inizializzazione e calcolo
    
    pesi = randi([pesoMin pesoMax], numPazienti, 1);    %Inizializza i pesi dei pazienti in un vettore colonna
    volumi = pesi ./ dens_tessuto;                      % Calcola il volume del tessuto per paziente (m^3)
    numCellule = volumi * dens_cellule;                 % Calcola il numero di cellule per paziente
    
    superficieCellula = 4 * pi * (diam_cellula/2)^2;    % Calcola la superficie di una cellula (m^2)

    lipidiTot = dens_lipidi * superficieCellula .* numCellule;  % numero di lipidi per paziente
    lipidiGr = (lipidiTot / N_Avogadro) * PM_lipidi;            % grammi di lipidi per paziente

    %% Plot istogramma

    edges = minL:10:maxL;  % intervalli con passo 10
    histogram(lipidiGr, edges)
    xlabel('Grammi di lipidi')
    ylabel('Numero di pazienti')
    title('Distribuzione grammi di lipidi per paziente')
    grid on;

    % Verifica se alcuni valori non sono inclusi nel grafico
    if min(lipidiGr) < minL || max(lipidiGr) > maxL
        warning('Alcuni valori di lipidi non sono riportati sul grafico. Intervallo calcolato: [%.2f, %.2f] g', valoreMin, valoreMax);
    end

end
