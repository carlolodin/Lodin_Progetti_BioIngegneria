load('SamplEEG.mat')
data = SamplEEG(1, :); % Prendo il primo canale
fs = sampling_rate; % Frequenza di campionamento
DATA = fftshift(fft(data)); % Trasformata di Fourier del segnale

% IIR filter evaluation for 1-8 Hz bandpass filter

%% Parametri globali
banda = [8 16];                   % Banda di interesse
N = length(SamplEEG(1, :));     % Numero di campioni
T = 1/sampling_rate;            % Periodo di campionamento
freqs = linspace(-(N/2)/(N*T), ((N/2)-1)/(N*T), size(SamplEEG, 2)); % Asse delle frequenze
t = linspace(0, (N-1)*T, N); % Asse del tempo

%% FIltri IIR
%% Butterworth method
ordini = [2 10 14 18 26 60 100 150 200]; % Ordini da testare - MIGLIORE ORDINE 60
butter_filters = cell(length(ordini),1); % Cell array per i filtri equiripple

% Progettazione dei filtri
for k = 1:length(ordini)
    n = ordini(k);
    tic;
    disp(['Progettazione filtro butterworth di ordine ' num2str(ordini(k))]);
    butter_filters{k} = designfilt("bandpassiir", FilterOrder= ordini(k), HalfPowerFrequency1=banda(1)-.1, ...
        HalfPowerFrequency2=banda(2)+.1, SampleRate=sampling_rate, DesignMethod='butter');   
    toc;
end

%Verifica stabilità dei filtri
for k = 1:length(ordini)
    if ~isstable(butter_filters{k})
        warning(['IL filtro Butterworth di ordine ' num2str(ordini(k)) ' NON è stabile.']);
    end
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;


for k = 1:length(ordini)
    [h, f] = freqz(butter_filters{k}, length(freqs), sampling_rate);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end

xlim([0 18]);   % Zoom sul range di interesse
ylim([-50 5]); % Zoom sull'asse y
plot([16 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;

title('Modulo della risposta in frequenza dei filtri IIR (Cheby1 Method)');
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (dB)');
legend show;

%GRAFICI DELLO SPETTRO DEL SEGNALE ORIGINALE E FILTRATO
figure;
subplot(5,2,1)  
hold on;
plot(freqs, fftshift(abs(fft(double(SamplEEG(1, :))))/N)); % TDF del segnale originale
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (V e-6)');
title("Segnale originale");
xlim([-18, 18]);
plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');

for i = 1:length(ordini)   
    y = filter(butter_filters{i}, data);
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18, 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri IIR di ordine crescente tramite il metodo Butterworth. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 5000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

%% Cheby1 method
ordini = [2 10 14 18 26 60 100 150 200]; % Ordini da testare - MIGLIORE ORDINE 26
cheby1_filters = cell(length(ordini),1); % Cell array per i filtri equiripple

% Progettazione dei filtri
for k = 1:length(ordini)
    n = ordini(k);
    tic;
    disp(['Progettazione filtro Cheby1 di ordine ' num2str(ordini(k))]);
    cheby1_filters{k} = designfilt("bandpassiir", FilterOrder= ordini(k), PassbandFrequency1=banda(1), ...
        PassbandFrequency2=banda(2), SampleRate=sampling_rate, DesignMethod='cheby1', PassBandRipple=0.1);   
    toc;
end

%Verifica stabilità dei filtri
for k = 1:length(ordini)
    if ~isstable(cheby1_filters{k})
        warning(['IL filtro Cheby1 di ordine ' num2str(ordini(k)) ' NON è stabile.']);
    end
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;


for k = 1:length(ordini)
    [h, f] = freqz(cheby1_filters{k}, length(freqs), sampling_rate);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end

xlim([0 18]);   % Zoom sul range di interesse
ylim([-50 5]); % Zoom sull'asse y
plot([16 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;

title('Modulo della risposta in frequenza dei filtri IIR (Butterworth Method)');
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (dB)');
legend show;

%GRAFICI DELLO SPETTRO DEL SEGNALE ORIGINALE E FILTRATO
figure;
subplot(5,2,1)  
hold on;
plot(freqs, fftshift(abs(fft(double(SamplEEG(1, :))))/N)); % TDF del segnale originale
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (V e-6)');
title("Segnale originale");
xlim([-18, 18]);
plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');

for i = 1:length(ordini)   
    y = filter(cheby1_filters{i}, data);
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18, 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri IIR di ordine crescente tramite il metodo Cheby1. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 5000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

%% Cheby2 method
ordini = [2 10 14 18 26 60 100 150 200]; % Ordini da testare - MIGLIORE ORDINE 60
cheby2_filters = cell(length(ordini),1); % Cell array per i filtri equiripple

% Progettazione dei filtri
for k = 1:length(ordini)
    n = ordini(k);
    tic;
    disp(['Progettazione filtro cheby2 di ordine ' num2str(ordini(k))]);
    cheby2_filters{k} = designfilt("bandpassiir", FilterOrder= ordini(k), StopbandFrequency1=banda(1)-.1, ...
        StopbandFrequency2=banda(2)+.1, SampleRate=sampling_rate, DesignMethod='cheby2', StopBandAttenuation=6);   
    toc;
end

%Verifica stabilità dei filtri
for k = 1:length(ordini)
    if ~isstable(cheby2_filters{k})
        warning(['Il filtro Cheby2 di ordine ' num2str(ordini(k)) ' NON è stabile.']);
    end
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;


for k = 1:length(ordini)
    [h, f] = freqz(cheby2_filters{k}, length(freqs), sampling_rate);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end

xlim([0 18]);   % Zoom sul range di interesse
ylim([-50 5]); % Zoom sull'asse y
plot([16 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;

title('Modulo della risposta in frequenza dei filtri IIR (Cheby2 Method)');
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (dB)');
legend show;

%GRAFICI DELLO SPETTRO DEL SEGNALE ORIGINALE E FILTRATO
figure;
subplot(5,2,1)  
hold on;
plot(freqs, fftshift(abs(fft(double(SamplEEG(1, :))))/N)); % TDF del segnale originale
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (V e-6)');
title("Segnale originale");
xlim([-18, 18]);
plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');

for i = 1:length(ordini)   
    y = filter(butter_filters{i}, data);
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18, 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri IIR di ordine crescente tramite il metodo Cheby1. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 5000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);
        
%% Elliptic method
ordini = [2 6 10 14 18 26 50]; % Ordini da testare - MIGLIORE ORDINE 14
elliptic_filters = cell(length(ordini),1); % Cell array per i filtri equiripple

% Progettazione dei filtri
for k = 1:length(ordini)
    n = ordini(k);
    tic;
    disp(['Progettazione filtro elliptic di ordine ' num2str(ordini(k))]);
    elliptic_filters{k} = designfilt("bandpassiir", FilterOrder= ordini(k), PassbandFrequency1=banda(1), ...
        PassbandFrequency2=banda(2), SampleRate=sampling_rate, DesignMethod='ellip', StopbandAttenuation1=20, PassbandRipple=0.1, ...
         StopbandAttenuation2=6);   
    toc;
end

%Verifica stabilità dei filtri
for k = 1:length(ordini)
    if ~isstable(elliptic_filters{k})
        warning(['IL filtro Elliptic di ordine ' num2str(ordini(k)) ' NON è stabile.']);
    end
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;


for k = 1:length(ordini)
    [h, f] = freqz(elliptic_filters{k}, length(freqs), sampling_rate);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end

xlim([0 18]);   % Zoom sul range di interesse
ylim([-50 5]); % Zoom sull'asse y
plot([16 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;

title('Modulo della risposta in frequenza dei filtri IIR (Elliptic Method)');
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (dB)');
legend show;

%GRAFICI DELLO SPETTRO DEL SEGNALE ORIGINALE E FILTRATO
figure;
subplot(5,2,1)  
hold on;
plot(freqs, fftshift(abs(fft(double(SamplEEG(1, :))))/N)); % TDF del segnale originale
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (V e-6)');
title("Segnale originale");
xlim([-18, 18]);
plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');

for i = 1:length(ordini)   
    y = filter(elliptic_filters{i}, data);
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18, 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri IIR di ordine crescente tramite il metodo Elliptic. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 5000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

%% Confronto tra i metodi IIR
figure;
hold on;
plot(t, SamplEEG(1, :), 'DisplayName', 'Segnale originale'); % Segnale originale
plot(t, filter(butter_filters{6}, data), 'DisplayName', 'IIR Butter Order 60'); % Filtro FIR con metodo della finestra
plot(t, filter(cheby1_filters{5}, data), 'DisplayName', 'IIR Cheby1 Order 26'); % Filtro FIR con metodo equiripple
plot(t, filter(cheby2_filters{6}, data), 'DisplayName', 'IIR Cheby2 Order 60'); % Filtro FIR con metodo least-squares
plot(t, filter(elliptic_filters{4}, data), 'DisplayName', 'IIR Elliptic Order 14'); % Filtro FIR con metodo least-squares
xlabel('Tempo (s)');
ylabel('Ampiezza (V e-6)');
title('Confronto tra i metodi FIR');
grid on;
legend show;
xlim([0 3]); % Zoom sul range di interesse

%% Risultati FIR
% Tutti i filtri mostrano comportamenti peggiori rispetto ai filtri FIR,
% Tra questi, il filtro cheby2 di ordine 60 mostra un buon compromesso tra prestazioni e dimensioni,
% ed ha la distorsione di fase minore ma non raggiunge le performance del filtro FIR di ordine 3000.
