load("SampleHEP.mat")
T = t(2) - t(1); % Sampling period

%Coorrect an error in R_peaks
for i=1:length(R_peaks)
    if R_peaks(i) == floor(58.778/T)
        break;
    end
end
R_peaks(i) = floor((58.778/T) - ((58.778 - 58.634)/T));

%Plot ECG signal with R-peaks
figure;
hold on;
plot(t, ECG);
scatter(R_peaks*T, ECG(R_peaks), 'rx');
title('ECG Signal with R-peaks');
ylabel('ECG Amplitude (Ve-6)');
xlabel('Time (seconds)');

% compute HRV
HRV = zeros(1, length(R_peaks) - 1);
for i = 1:length(R_peaks) - 1
    HRV(i) = (R_peaks(i + 1) - R_peaks(i)) * T; % Convert indices to time
end

average_BPM = 60/mean(HRV);

% Plot HRV
figure;
hold on;
plot(HRV, 'o-');
title('Heart Rate Variability (HRV)');
xlabel('# Peak');
ylabel('Time (seconds)');
text(0.5, 0.9, ['Average BPM: ' num2str(average_BPM)], 'Units', 'normalized', 'HorizontalAlignment', 'center');
grid on;

HEP = zeros(length(0:T:1));
for hep_idx = 1:length(HEP)
    for peak_idx=1:length(R_peaks)
        idx = R_peaks(peak_idx) + hep_idx - 1;
        if idx >= 0 && idx <= length(C4)
            HEP(hep_idx) = HEP(hep_idx) + C4(idx);
        end
    end
    HEP(hep_idx) = HEP(hep_idx) / length(R_peaks);
end

figure;
hold on;
plot(0:T:1, HEP);
title('Heartbeat Evoked Potential (HEP)');
xlabel('Time (seconds)');
ylabel('Amplitude (Ve-6)');