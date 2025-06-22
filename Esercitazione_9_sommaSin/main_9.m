%% Esercitazione 9 - Somma di Sinusoidi

%% Parametri delle funzioni

A1 = 3;     % Ampiezze
A2 = 4;
f = @(x)A1*sin(2*pi*2*x)+A2*sin(2*pi*2.1*x); % Funzione da analizzare
Tc = 1/10; % Periodo di campionamento

%% Studio con N Fissato - T variabile con zero padding per arrivare a N1 campioni
T1 = [3, 5, 10, 15, 20, 30];
N1 = 300;
Y1 = zeros(length(T1), N1);

figure;
hold on;

frequencies = linspace(-5, ((N1/2)-1)/(N1*Tc), N1); % Frequenze per il grafico
for i = 1:length(T1)
    
    x = 0:Tc:T1(i); 
    y = f(x);
    Y1(i, :) = fftshift(fft(y, N1)/N1);
    subplot(3, 2, i);
    plot(frequencies, abs(Y1(i, :)));
    title(['T = ' num2str(T1(i))]);
    
end

%% Studio con T Fissato
T2 = 15;
N2 = [30, 50, 100, 150, 200, 300];
Y2 = cell(1, 6);

figure;
hold on;

x = 0:Tc:T2-1; % Common x for all N2
y = f(x);

for i = 1:length(N2)
    frequencies = linspace(-5, ((N2(i)/2)-1)/(Tc*N2(i)), N2(i));
    Y2{i} = fftshift(fft(y, N2(i))/N2(i));
    
    subplot(3, 2, i);
    plot(frequencies, abs(Y2{i}));
    title(['N = ' num2str(N2(i))]);
end




