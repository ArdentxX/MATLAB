% odbiornik FM: P. Swiatkiewicz, T. Twardowski, T. Zielinski, J. Bułat

clear all; close all;

fs = 3.2e6;         % sampling frequency
N  = 32e6;         % number of samples (IQ)
fc = 0.50e6;        % central frequency of MF station

bwSERV = 80e3;     % bandwidth of an FM service (bandwidth ~= sampling frequency!)
bwAUDIO = 16e3;     % bandwidth of an FM audio (bandwidth == 1/2 * sampling frequency!)

f = fopen('samples_100MHz_fs3200kHz.raw');
s = fread(f, 2*N, 'uint8');
fclose(f);

s = s-127;

% IQ --> complex
wideband_signal = s(1:2:end) + sqrt(-1)*s(2:2:end); clear s;

% Narysuj widmo gęstości mocy sygnału z tunera
figure;
psd(spectrum.welch('Hamming',1024), wideband_signal(1:min(length(wideband_signal),2^18)),'Fs',fs);
title('Widmo gęstości mocy sygnału z tunera');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% Extract carrier of selected service, then shift in frequency the selected service to the baseband
wideband_signal_shifted = wideband_signal .* exp(-sqrt(-1)*2*pi*fc/fs*[0:N-1]');

% Narysuj widmo po przesunięciu
figure;
psd(spectrum.welch('Hamming',1024), wideband_signal_shifted(1:min(length(wideband_signal_shifted),2^18)),'Fs',fs);
title('Widmo gęstości mocy sygnału po przesunięciu do 0 Hz');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% Filter out the service from the wide-band signal
% Implementacja filtru Butterworth LP 4 rzędu o częstotliwości granicznej 80 kHz
[b, a] = butter(4, 2*80e3/fs);
wideband_signal_filtered = filter(b, a, wideband_signal_shifted);

% Narysuj widmo po filtracji LP
figure;
psd(spectrum.welch('Hamming',1024), wideband_signal_filtered(1:min(length(wideband_signal_filtered),2^18)),'Fs',fs);
title('Widmo gęstości mocy sygnału po filtracji LP');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% Down-sample to service bandwidth - bwSERV = new sampling rate
x = wideband_signal_filtered(1:fs/bwSERV:end);

% Narysuj widmo po zmianie częstotliwości próbkowania
figure;
psd(spectrum.welch('Hamming',1024), x(1:min(length(x),2^18)),'Fs',bwSERV);
title('Widmo gęstości mocy sygnału po decymacji do 80kHz');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% FM demodulation
dx = x(2:end).*conj(x(1:end-1));
y = atan2(imag(dx), real(dx));

% Narysuj widmo po demodulacji FM
figure;
psd(spectrum.welch('Hamming',1024), y(1:min(length(y),2^18)),'Fs',bwSERV);
title('Widmo gęstości mocy sygnału po demodulacji FM');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% Decimate to audio signal bandwidth bwAUDIO
% Implementacja filtru antyaliasingowego przed decymacją
[b_lp, a_lp] = butter(4, 2*16e3/bwSERV);
y = filter(b_lp, a_lp, y);

% Narysuj widmo po filtracji antyaliasingowej
figure;
psd(spectrum.welch('Hamming',1024), y(1:min(length(y),2^18)),'Fs',bwSERV);
title('Widmo gęstości mocy sygnału po filtrze antyaliasingowym');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

ym = y(1:bwSERV/bwAUDIO:end);  % decimate (1/5)

% Narysuj widmo po decymacji do 16kHz
figure;
psd(spectrum.welch('Hamming',1024), ym(1:min(length(ym),2^18)),'Fs',bwAUDIO);
title('Widmo gęstości mocy sygnału po decymacji do 16kHz');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% De-emfaza - filtr LP 1 rzedu z granicą 2.1 kHz
[b_de, a_de] = butter(1, 2*2.1e3/bwAUDIO);
ym_de = filter(b_de, a_de, ym);

% Narysuj widmo po de-emfazie
figure;
psd(spectrum.welch('Hamming',1024), ym_de(1:min(length(ym_de),2^18)),'Fs',bwAUDIO);
title('Widmo gęstości mocy sygnału po de-emfazie');
xlabel('Częstotliwość (Hz)');
ylabel('Gęstość mocy (dB/Hz)');

% Listen to the final result
ym_de = ym_de - mean(ym_de);
ym_de = ym_de / (1.001*max(abs(ym_de)));
soundsc(ym_de, bwAUDIO);

% Zadanie opcjonalne - wizualizacja spektrogramu sygnału
figure;
spectrogram(ym_de, hamming(512), 256, 1024, bwAUDIO, 'yaxis');
title('Spektrogram sygnału audio po dekodowaniu');

% Zadanie opcjonalne - filtr de-emfazy z charakterystyką płaską do 2.1 kHz i spadkiem 20dB/dekadę
fc_de = 2.1e3;
[b_de_opt, a_de_opt] = butter(1, 2*fc_de/bwAUDIO);

% Wizualizacja charakterystyki filtru de-emfazy
figure;
freqz(b_de_opt, a_de_opt, 1024, bwAUDIO);
title('Charakterystyka filtru de-emfazy');

% Pre-emfaza - odwrotny do de-emfazy (filtr dla nadajnika)
% Zamieniamy biegun na zero i na odwrót
b_pre = a_de_opt;
a_pre = b_de_opt;

% Wizualizacja charakterystyki filtru pre-emfazy
figure;
freqz(b_pre, a_pre, 1024, bwAUDIO);
title('Charakterystyka filtru pre-emfazy');

% Test filtracji pre-emfazy i de-emfazy w kaskadzie
ym_pre = filter(b_pre, a_pre, ym);
ym_cascade = filter(b_de_opt, a_de_opt, ym_pre);

% Porównanie sygnału oryginalnego i po kaskadzie filtrów
figure;
subplot(2,1,1);
plot(ym(1:1000));
title('Sygnał oryginalny');
subplot(2,1,2);
plot(ym_cascade(1:1000));
title('Sygnał po pre-emfazie i de-emfazie');
