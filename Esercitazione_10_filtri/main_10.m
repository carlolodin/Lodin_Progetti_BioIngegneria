% FIR: Least squares ordine 3000
% IIR: Cheby2 ordine 60

load('SamplEEG.mat')
data = SamplEEG(1, :); % Prendo il primo canale
fs = sampling_rate; % Frequenza di campionamento
DATA = fftshift(fft(data, 17999)); % Trasformata di Fourier del segnale

% IIR filter evaluation for 1-8 Hz bandpass filter

%% Parametri globali
banda_1 = [1 8];                   % Banda di interesse
banda_alpha = [8 16];
N = length(SamplEEG(1, :));     % Numero di campioni
T = 1/sampling_rate;            % Periodo di campionamento
freqs = linspace(-(N/2)/(N*T), ((N/2)-1)/(N*T), size(SamplEEG, 2)); % Asse delle frequenze
t = linspace(0, (N-1)*T, N); % Asse del tempo

%% Progettazione dei filtri
iir_filt_1 = designfilt("bandpassiir", FilterOrder = 60, StopbandFrequency1=banda_1(1)-.1, ...
        StopbandFrequency2=banda_1(2)+.1, SampleRate=sampling_rate, DesignMethod='cheby2', StopBandAttenuation=6);
iir_filt_alpha = designfilt("bandpassiir", FilterOrder = 60, StopbandFrequency1=banda_alpha(1)-.1, ...
        StopbandFrequency2=banda_alpha(2)+.1, SampleRate=sampling_rate, DesignMethod='cheby2', StopBandAttenuation=6);

fir_filt_1 = designfilt('bandpassfir', 'FilterOrder', 3000, ...
        'PassbandFrequency1', banda_1(1), 'PassbandFrequency2', ...
        banda_1(2), 'StopbandFrequency1',banda_1(1)-.1,  'Stopbandfrequency2', banda_1(2) +.1 , ...
        'SampleRate', sampling_rate, 'DesignMethod', 'ls');
fir_filt_alpha = designfilt('bandpassfir', 'FilterOrder', 3000, ...
        'PassbandFrequency1', banda_alpha(1), 'PassbandFrequency2', ...
        banda_alpha(2), 'StopbandFrequency1',banda_alpha(1)-.1,  'Stopbandfrequency2', banda_alpha(2) +.1 , ...
        'SampleRate', sampling_rate, 'DesignMethod', 'ls');

%% Applicazione dei fritri
data_iir_filt_1 = filter(iir_filt_1, data);
data_iir_filt_alpha = filter(iir_filt_alpha, data);

data_fir_filt_1_conv = conv(fir_filt_1.Coefficients, data);
data_fir_filt_alpha_conv = conv(fir_filt_alpha.Coefficients, data);

FIR_FILT_1 = fftshift(fft(fir_filt_1.Coefficients, 17999));
FIR_FILT_ALPHA = fftshift(fft(fir_filt_alpha.Coefficients, 17999));

DATA_FIR_FILT_1_PROD = DATA .* FIR_FILT_1;
DATA_FIR_FILT_ALPHA_PROD = DATA .* FIR_FILT_ALPHA;

data_fir_filt_1_prod = ifft(ifftshift(DATA_FIR_FILT_1_PROD), 'symmetric');
data_fir_filt_alpha_prod = ifft(ifftshift(DATA_FIR_FILT_ALPHA_PROD), 'symmetric');

figure;
hold on;
plot(t, data, 'k', 'DisplayName', 'Original Signal');
plot(t, data_iir_filt_1, 'r', 'DisplayName', 'IIR 1-8 Hz');
plot(t, data_fir_filt_1_conv(1:15000), 'b', 'DisplayName', 'FIR 1-8 Hz (Conv)');
plot(t, data_fir_filt_1_prod(1:15000), 'g', 'DisplayName', 'FIR 1-8 Hz (Prod)');
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal Filtering Comparison (1-8 Hz)');
legend('show');
grid on;

figure;
hold on;
plot(t, data, 'k', 'DisplayName', 'Original Signal');
plot(t, data_iir_filt_alpha, 'r', 'DisplayName', 'IIR 8-16 Hz');
plot(t, data_fir_filt_alpha_conv(1:15000), 'b', 'DisplayName', 'FIR 8-16 Hz (Conv)');
plot(t, data_fir_filt_alpha_prod(1:15000), 'g', 'DisplayName', 'FIR 8-16 Hz (Prod)');
xlabel('Time (s)'); 
ylabel('Amplitude');
title('Signal Filtering Comparison (8-16 Hz)');
legend('show');
grid on;

%% Spettro dei segnali filtrati
DATA_15000 = fftshift(fft(data, 15000)/15000);
DATA_IIR_FILT_1 = fftshift(fft(data_iir_filt_1, 15000)/15000);
DATA_IIR_FILT_ALPHA = fftshift(fft(data_iir_filt_alpha, 15000)/15000);
DATA_FIR_FILT_1_CONV = fftshift(fft(data_fir_filt_1_conv(1:15000), 15000)/15000);
DATA_FIR_FILT_ALPHA_CONV = fftshift(fft(data_fir_filt_alpha_conv(1:15000), 15000)/15000);
DATA_FIR_FILT_1_PROD = fftshift(fft(data_fir_filt_1_prod(1:15000), 15000)/15000);
DATA_FIR_FILT_ALPHA_PROD = fftshift(fft(data_fir_filt_alpha_prod(1:15000), 15000)/15000);

figure;
hold on;
plot(freqs, abs(DATA_15000), 'k', 'DisplayName', 'Original Signal');
plot(freqs, abs(DATA_IIR_FILT_1), 'r', 'DisplayName', 'IIR 1-8 Hz');
plot(freqs, abs(DATA_FIR_FILT_1_CONV), 'b', 'DisplayName', 'FIR 1-8 Hz (Conv)');
plot(freqs, abs(DATA_FIR_FILT_1_PROD), 'g', 'DisplayName', 'FIR 1-8 Hz (Prod)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Filtered Signals (1-8 Hz)');
legend('show');
xlim([-10 10])

figure;
hold on;
plot(freqs, abs(DATA_15000), 'k', 'DisplayName', 'Original Signal');
plot(freqs, abs(DATA_IIR_FILT_ALPHA), 'r', 'DisplayName', 'IIR 8-16 Hz');
plot(freqs, abs(DATA_FIR_FILT_ALPHA_CONV), 'b', 'DisplayName', 'FIR 8-16 Hz (Conv)');
plot(freqs, abs(DATA_FIR_FILT_ALPHA_PROD), 'g', 'DisplayName', 'FIR 8-16 Hz (Prod)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Spectrum of Filtered Signals (1-8 Hz)');
legend('show');
xlim([-18 18])

%% Applicazione filtri a tutti i canali
data_iir_filt_alpha_all = zeros(size(SamplEEG));
data_fir_filt_alpha_all = zeros(size(SamplEEG));

for i = 1:size(SamplEEG, 1)
    data_iir_filt_alpha_all(i, :) = filter(iir_filt_alpha, SamplEEG(i, :));
    
    data_fir_filt_alpha_all(i, :) = conv(SamplEEG(i, :), fir_filt_alpha.Coefficients, 'same');
end

% Calcolo della potenza delle righe di data_fir_filt_alpha_all
potenza_fir_filt_alpha = sum(data_fir_filt_alpha_all.^2, 2)./15000*sampling_rate;
potenza_iir_filt_alpha = sum(data_iir_filt_alpha_all.^2, 2)./15000*sampling_rate;

comparison = table(Chans_names(:), potenza_fir_filt_alpha(:), potenza_iir_filt_alpha(:), ...
    'VariableNames', {'Channel', 'FIR_Alpha_Power', 'IIR_Alpha_Power'});
comparison = sortrows(comparison, 'FIR_Alpha_Power', 'descend');

catChan = categorical(comparison.Channel);
catChan = reordercats(catChan, string(comparison.Channel)); % Mantieni l'ordine della tabella

figure;
bar(catChan, [comparison.FIR_Alpha_Power, comparison.IIR_Alpha_Power]);
xlabel('Channel');
ylabel('Power');
legend('FIR Alpha Power', 'IIR Alpha Power');
title('Alpha Band Power by Channel');
grid on;



