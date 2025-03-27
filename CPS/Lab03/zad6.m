clearvars; close all; clc;

% Wczytanie dźwięków
[canary, Fs1] = audioread('canary.wav');
[motor, Fs2] = audioread('motor.wav');

% Dopasowanie częstotliwości próbkowania
if Fs1 ~= Fs2
    targetFs = min(Fs1, Fs2); % Wybór niższej częstotliwości
    canary = resample(canary, targetFs, Fs1);
    motor = resample(motor, targetFs, Fs2);
    Fs = targetFs;
else
    Fs = Fs1;
end

% Dopasowanie długości sygnałów
len = max(length(canary), length(motor));
canary = [canary; zeros(len - length(canary), 1)];
motor = [motor; zeros(len - length(motor), 1)];

% Oś czasu
time = (0:len-1) / Fs;

% Wykresy oryginalnych sygnałów
figure;
subplot(2,1,1);
plot(time, canary);
title('Oryginalny sygnał - Śpiew ptaka'); xlabel('Czas (s)'); ylabel('Amplituda');
subplot(2,1,2);
plot(time, motor);
title('Oryginalny sygnał - Warkot silnika'); xlabel('Czas (s)'); ylabel('Amplituda');

% Obliczenie DFT
CanaryDFT = fft(canary);
MotorDFT = fft(motor);

% Oś częstotliwości
freqs = (-len/2:len/2-1) * (Fs / len);

% Wykresy widm
figure;
subplot(2,1,1);
plot(freqs, abs(fftshift(CanaryDFT)));
title('Widmo DFT - Śpiew ptaka'); xlabel('Częstotliwość (Hz)'); ylabel('Amplituda');
subplot(2,1,2);
plot(freqs, abs(fftshift(MotorDFT)));
title('Widmo DFT - Warkot silnika'); xlabel('Częstotliwość (Hz)'); ylabel('Amplituda');

% Suma sygnałów
sum_signal = canary + motor;
SumDFT = fft(sum_signal);

% Wykres widma sumy
figure;
plot(freqs, abs(fftshift(SumDFT)));
title('Widmo DFT - Suma sygnałów'); xlabel('Częstotliwość (Hz)'); ylabel('Amplituda');

% Filtracja - usunięcie niskich częstotliwości (warkot silnika)
thresh = 1000; % Próg częstotliwości do usunięcia

SumDFT_shifted = fftshift(SumDFT);
idx = (freqs > -thresh) & (freqs < thresh);
SumDFT_shifted(idx) = 0;
SumDFT_shifted(:,2)=0;
SumDFT_filtered = ifftshift(SumDFT_shifted);

% IDFT
filtered_signal = real(ifft(SumDFT_filtered));

% Odsłuch
sound(filtered_signal, Fs);

% Wykres przefiltrowanego sygnału
figure(4);
plot(time, filtered_signal); 
xlabel('Czas (s)'); ylabel('Amplituda');
title('Przefiltrowany sygnał');

% Zapis do pliku
audiowrite('filtered_signal.wav', filtered_signal, Fs);
