fir = double(out.simout.Data);
FIR = fftshift(fft(fir)/length(fir));

iir = double(out.simout1.Data);
IIR = fftshift(fft(iir)/length(iir));

f = linspace(-500/(0.001*1001), 500/(0.001*1001), length(fir));

figure;
subplot(2,1,1);
hold on;
plot(f, abs(FIR), 'DisplayName', 'FIR');
plot(f, abs(IIR), 'DisplayName', 'IIR');
legend show;
xlabel('frequenza (Hz)')
ylabel('Ampiezza')

subplot(2,1,2);
hold on;
plot(f, angle(FIR), 'DisplayName', 'FIR');
plot(f, angle(IIR), 'DisplayName', 'IIR');
legend show;
xlabel('frequenza (Hz)')
ylabel('Fase (Rad)')