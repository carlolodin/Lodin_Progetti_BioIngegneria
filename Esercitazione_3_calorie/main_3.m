data = readmatrix("pesi.xlsx", "Range", "A2:F51"); %Import to a 50x6 matrix
save("Esercitazione_3_calorie/workspace.mat")

load("Esercitazione_3_calorie/workspace.mat"); %Load the matrix from the file

fasce_eta = [30, 30, 14, 14, 5, 5]';
sessi = ['M', 'F', 'M', 'F', 'M', 'F']';
percent = [.15, .30, .55]; %Proteine, Grassi, Carboidrati
media_gruppi = zeros(6, 3);
risultati = zeros(6, 50, 3); %Gruppo, Soggetto, Nutriente (1=Proteine, 2=Grassi, 3=Carboidrati)
kcal = zeros(50);

figure;
for gruppo = 1:6
    eta = fasce_eta(gruppo);
    sesso = sessi(gruppo);
    pesi = data(:, gruppo);

    for soggetto = 1:50
        kcal = calcola_kcal(eta, sesso, pesi(soggetto));
        risultati(gruppo, soggetto, :) = kcal * percent; % Proteine
    end
    media_gruppi(gruppo, :) = mean(risultati(gruppo, :, :), 2); % Calcola la media per il gruppo

    % Subplot per gruppo
    subplot(2, 3, gruppo);
    bar(squeeze(risultati(gruppo, :, :)), 'stacked');
    title(['Gruppo ', num2str(gruppo) ' (' sesso ', età mediana=' num2str(eta) ')']);
    xlabel('Pazienti'); ylabel('kCal');
    legend({'Proteine', 'Lipidi', 'Carboidrati'});

end

figure;
bar(media_gruppi, 'stacked');
set(gca, 'XTickLabel', {'M18-60', 'F18-60', 'M10-18', 'F10-18', 'M3-10', 'F3-10'});
xlabel('Gruppi');
ylabel('kCal');
title('Media apporto calorico per nutrienti');
legend({'Proteine', 'Lipidi', 'Carboidrati'});
