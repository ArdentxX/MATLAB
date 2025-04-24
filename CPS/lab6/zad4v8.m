clc; clear; close all;
% ZADANIE 4 - Filtrowanie dźwięków rzeczywistych

% Wczytanie dwóch nagrań
[x1, fs] = audioread('mowa.wav');
[x2, ~] = audioread('ptak.wav');

% Dopasowanie długości sygnałów
L = min(length(x1), length(x2));
x1 = x1(1:L);
x2 = x2(1:L);

% Mieszanie sygnałów
x_mix = x1 + x2;

% Analiza częstotliwościowa (FFT) 
X1 = abs(fft(x1));
X2 = abs(fft(x2));
XM = abs(fft(x_mix));

% Rysowanie widm FFT
figure;
plot(linspace(0, fs/2, L/2), X1(1:L/2));
title('Widmo FFT sygnału mowa.wav');
xlabel('Częstotliwość (Hz)');

figure;
plot(linspace(0, fs/2, L/2), X2(1:L/2));
title('Widmo FFT sygnału ptak.wav');
xlabel('Częstotliwość (Hz)');

figure;
plot(linspace(0, fs/2, L/2), XM(1:L/2));
title('Widmo FFT mieszaniny');
xlabel('Częstotliwość (Hz)');

% Poprawione spektrogramy
figure;
pspectrum(x1, fs, 'spectrogram', 'FrequencyLimits', [0 8000], 'TimeResolution', 0.1);
title('Spektrogram sygnału mowa.wav');

figure;
pspectrum(x2, fs, 'spectrogram', 'FrequencyLimits', [0 8000], 'TimeResolution', 0.1);
title('Spektrogram sygnału ptak.wav');

figure;
pspectrum(x_mix, fs, 'spectrogram', 'FrequencyLimits', [0 8000], 'TimeResolution', 0.1);
title('Spektrogram mieszaniny');

% Projekt filtru IIR - pasmowo-przepustowy lub LP dla mowy (do 3 kHz)
[b, a] = butter(4, 2*3000/fs);  % LP 4 rzedu

% Filtracja sygnalu mieszanego
x_filtered = filter(b, a, x_mix);

% Odsłuch i zapis
x_filtered = x_filtered - mean(x_filtered);
x_filtered = x_filtered / max(abs(x_filtered));
soundsc(x_filtered, fs);

audiowrite('mowa_odfiltrowana.wav', x_filtered, fs);

figure;
plot(linspace(0, fs/2, L/2), x_filtered(1:L/2));
title('Widmo FFT mieszaniny po filtrze');
xlabel('Częstotliwość (Hz)');

figure;
pspectrum(x_filtered, fs, 'spectrogram', 'FrequencyLimits', [0 8000], 'TimeResolution', 0.1);
title('Spektrogram - po filtracji');

% Zera i bieguny filtru
figure;
zplane(b, a);
title('Zera i bieguny filtru');

% Charakterystyka czestotliwosciowa
fvtool(b, a);
