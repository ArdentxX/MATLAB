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

% Obliczenie DFT
CanaryDFT = fft(canary);
MotorDFT = fft(motor);

% Oś częstotliwości
freqs = linspace(-Fs/2, Fs/2, len);

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
thresh = 500; % Próg częstotliwości do usunięcia

SumDFT_filtered = SumDFT;
idx = abs(freqs) < thresh;
SumDFT_filtered(idx) = 0;

% IDFT
filtered_signal = ifft(SumDFT_filtered, 'symmetric');

% Odsłuch
sound(filtered_signal, Fs);

% Wykres przefiltrowanego sygnału
figure;
plot(filtered_signal);
title('Przefiltrowany sygnał'); xlabel('Próbki'); ylabel('Amplituda');

% Zapis do pliku
audiowrite('filtered_signal.wav', filtered_signal, Fs);
