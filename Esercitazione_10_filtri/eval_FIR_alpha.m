load('SamplEEG.mat')

% FIR filter evaluation for 1-8 Hz bandpass filter

%% Parametri globali
banda = [8 16];                   % Banda di interesse
N = length(SamplEEG(1, :));     % Numero di campioni
T = 1/sampling_rate;            % Periodo di campionamento
freqs = linspace(-(N/2)/(N*T), ((N/2)-1)/(N*T), size(SamplEEG, 2)); % Asse delle frequenze
t = linspace(0, (N-1)*T, N); % Asse del tempo

%% FILTRI FIR
%% Metodo della finestra
ordini = [5 50 100 500 1000 3000 5000 50000]; % Ordini da testare
window_filters = cell(length(ordini),1);        % Cell array per i filri con metodo della finestra

% Progettazione dei filti
for k = 1:length(ordini)
    n = ordini(k);
    window_filters{k} = fir1(n, banda/(sampling_rate/2), 'bandpass'); % N.B. la frequenza è normalizzata 
                                                                     % rispetto alla frequenza di Nyquist
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;

responses_window = cell(length(ordini),1);  %Cell array per le risposte in frequenza

for k = 1:length(ordini)
    [h, f] = freqz(window_filters{k}, 1, length(freqs), sampling_rate);
    responses_window{k} = abs(h);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end

xlim([0 18]);   % Zoom sul range di interesse
plot([16 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;

title('Modulo della risposta in frequenza dei filtri FIR (Window Method)');
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (dB)');
legend show;

%GRAFICI DELLO SPETTRO DEL SEGNALE ORIGINALE E FILTRATO
figure;
subplot(5,2,1)  
plot(freqs, fftshift(abs(fft(double(SamplEEG(1, :))))/N)); % TDF del segnale originale
xlabel('Frequenza (Hz)');
ylabel('Ampiezza (V e-6)');
title("Segnale originale");
xlim([-18 18]);

for i = 1:length(ordini)   
    y = conv(double(SamplEEG(1, :)), window_filters{i}, 'same'); %Convoluzione con il filtro mantendo la dimensione del segnale originale
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri FIR di ordine crescente tramite il metodo della finestra. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 5000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

%% METODO EQUIRIPPLE %%%%%%%%%%%%%%%%%%%

ordini = [5 50 100 500 1000 3000 4000]; % Ordini da testare
equiripple_filters = cell(length(ordini),1); % Cell array per i filtri equiripple

% Progettazione dei filtri
for k = 1:length(ordini)
    n = ordini(k);
    tic;
    disp(['Progettazione filtro equiripple di ordine ' num2str(ordini(k))]);
    equiripple_filters{k} = designfilt('bandpassfir', 'FilterOrder', ordini(k), ...
        'PassbandFrequency1', banda(1), 'PassbandFrequency2', ...
        banda(2), 'StopbandFrequency1',banda(1)-.1,  'Stopbandfrequency2', banda(2) +.1 , ...
        'SampleRate', sampling_rate, 'DesignMethod', 'equiripple');   
    toc;
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;

responses_equiripple = cell(length(ordini),1);  %Cell array per le risposte in frequenza

for k = 1:length(ordini)
    [h, f] = freqz(equiripple_filters{k}.Coefficients, 1, length(freqs), sampling_rate);
    responses_equiripple{k} = abs(h);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end

xlim([0 18]);   % Zoom sul range di interesse
plot([1 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;

title('Modulo della risposta in frequenza dei filtri FIR (Equiripple Method)');
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
xlim([-18 18]);
plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');

for i = 1:length(ordini)   
    y = conv(double(SamplEEG(1, :)), equiripple_filters{i}.Coefficients, 'same'); %Convoluzione con il filtro mantendo la dimensione del segnale originale
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri FIR di ordine crescente tramite il metodo della finestra. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 5000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

%% Metodo least-squares
ordini = [5 50 100 500 1000 3000 4000]; % Ordini da testare
least_squares_filters = cell(length(ordini),1); % Cell array per i filtri least-squares 

% Progettazione dei filtri
for k = 1:length(ordini)
    n = ordini(k);
    tic;
    disp(['Progettazione filtro least-squares di ordine ' num2str(ordini(k))]);
    least_squares_filters{k} = designfilt('bandpassfir', 'FilterOrder', ordini(k), ...
        'PassbandFrequency1', banda(1), 'PassbandFrequency2', ...
        banda(2), 'StopbandFrequency1',banda(1)-.1,  'Stopbandfrequency2', banda(2) +.1 , ...
        'SampleRate', sampling_rate, 'DesignMethod', 'ls');   
    toc;
end

%GRAFICI DELLA RISPOSTA IN FREQUENZA
figure;
hold on;
responses_least_squares = cell(length(ordini),1);  %Cell array per le risposte in frequenza
for k = 1:length(ordini)
    [h, f] = freqz(least_squares_filters{k}.Coefficients, 1, length(freqs), sampling_rate);
    responses_least_squares{k} = abs(h);
    plot(f, 20*log10(abs(h)), 'DisplayName', ['Order ' num2str(ordini(k))]);
end
xlim([0 18]);   % Zoom sul range di interesse
plot([16 16], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '16 Hz');  %Plot linee guida della banda
plot([8 8], ylim, 'k--', 'LineWidth', 1.2, 'DisplayName', '8 Hz');
plot(xlim, [-3 -3], 'r--', 'LineWidth', 1.2, 'DisplayName', '-3 dB');
grid on;
title('Modulo della risposta in frequenza dei filtri FIR (Least-squares Method)');
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
xlim([-18 18]);
plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '1 Hz'); %Linee guida per la banda
plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
for i = 1:length(ordini)   
    y = conv(double(SamplEEG(1, :)), least_squares_filters{i}.Coefficients, 'same'); %Convoluzione con il filtro mantendo la dimensione del segnale originale
    Y = fftshift(fft(y)/N);   % TDF del segnale filtrato
    
    subplot(5,2, i+1)               % Plot modulo della risposta in frequenza
    hold on;
    plot(freqs, abs(Y))
    xlim([-18 18]);
    ylim([-.1, 4]);
    plot([16 16], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '16 Hz'); %Linee guida per la banda
    plot([8 8], ylim, 'r--', 'LineWidth', .8, 'DisplayName', '8 Hz');
    title(strcat("ordine:" , int2str(ordini(i))))
    xlabel('Frequenza (Hz)');
    ylabel('Ampiezza (V e-6)');
    grid on;
    
end
annotation('textbox', [0.1, 0.01, 0.8, 0.08], ...
            'String', 'Confronto tra la trasformata di Fourier del segnale originale e dei segnali filtrati con filtri FIR di ordine crescente tramite il metodo least-squares. Tra questi, il miglior compromesso tra performance e dimensioni è rappresentato dal filtro di ordine 3000', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 10);

%% Confronto tra i metodi FIR
figure;
hold on;
plot(t, SamplEEG(1, :), 'DisplayName', 'Segnale originale'); % Segnale originale
plot(t, conv(double(SamplEEG(1, :)), window_filters{6}, 'same'), 'DisplayName', 'FIR Window Order 5000'); % Filtro FIR con metodo della finestra
plot(t, conv(double(SamplEEG(1, :)), equiripple_filters{6}.Coefficients, 'same'), 'DisplayName', 'FIR Equiripple Order 3000'); % Filtro FIR con metodo equiripple
plot(t, conv(double(SamplEEG(1, :)), least_squares_filters{6}.Coefficients, 'same'), 'DisplayName', 'FIR Least-squares Order 3000'); % Filtro FIR con metodo least-squares
xlabel('Tempo (s)');
ylabel('Ampiezza (V e-6)');
title('Confronto tra i metodi FIR');
grid on;
legend show;
xlim([0 3]); % Zoom sul range di interesse

%% Risultati FIR
% Nel dominio del tempo, finestra e least-squares mostrano andamenti simili. 
% Tra i due, prediligo least-squares per il suo ordine più basso.

