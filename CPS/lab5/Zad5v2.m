clearvars; close all; clc;

% Częstotliwość środkowa (nośna) [Hz]
fc = 96e6;

% 1. FILTR TESTOWY (96 MHz ±1 MHz)
disp('Projektowanie filtru testowego (96 MHz ±1 MHz)');

% Parametry filtru testowego
f_pass1_test = fc - 1e6; % Dolna granica pasma przepustowego [Hz]
f_pass2_test = fc + 1e6; % Górna granica pasma przepustowego [Hz]
f_stop1_test = fc - 2e6; % Dolna granica pasma zaporowego [Hz]
f_stop2_test = fc + 2e6; % Górna granica pasma zaporowego [Hz]
Rp_test = 3;             % Maksymalne zafalowanie w paśmie przepustowym [dB]
Rs_test = 40;            % Minimalne tłumienie w paśmie zaporowym [dB]

% Normalizacja częstotliwości względem częstotliwości Nyquista
f_nyquist_test = fc + 5e6; % Przyjęta częstotliwość Nyquista [Hz]
Wp_test = [f_pass1_test f_pass2_test] / f_nyquist_test;
Ws_test = [f_stop1_test f_stop2_test] / f_nyquist_test;

% Projektowanie filtru Butterwortha
[n_test, Wn_test] = buttord(Wp_test, Ws_test, Rp_test, Rs_test);
[b_test, a_test] = butter(n_test, Wn_test, 'bandpass');

% Wyświetlenie parametrów filtru testowego
disp(['Rząd filtru testowego: ' num2str(n_test)]);
disp(['Znormalizowane częstotliwości graniczne: ' num2str(Wn_test)]);

% Charakterystyka częstotliwościowa filtru testowego
f_test = linspace(fc-5e6, fc+5e6, 10000);
w_test = 2*pi*f_test;
[h_test, ~] = freqs(b_test, a_test, w_test);
mag_test = 20*log10(abs(h_test));

% Rysowanie charakterystyki filtru testowego
figure('Name', 'Filtr testowy (96 MHz ±1 MHz)', 'Position', [100, 100, 900, 600]);
plot(f_test/1e6, mag_test, 'LineWidth', 1.5);
grid on;
xlabel('Częstotliwość [MHz]');
ylabel('Amplituda [dB]');
title(['Charakterystyka częstotliwościowa filtru testowego (rząd ' num2str(n_test) ')']);

% Zaznaczenie charakterystycznych punktów
hold on;
yline(-Rp_test, 'r--', 'Pasmo przepustowe (-3 dB)', 'LineWidth', 1.5);
yline(-Rs_test, 'g--', 'Pasmo zaporowe (-40 dB)', 'LineWidth', 1.5);
xline((fc-1e6)/1e6, 'b--', 'f_{c} - 1 MHz', 'LineWidth', 1.5);
xline((fc+1e6)/1e6, 'b--', 'f_{c} + 1 MHz', 'LineWidth', 1.5);
xline((fc-2e6)/1e6, 'm--', 'f_{c} - 2 MHz', 'LineWidth', 1.5);
xline((fc+2e6)/1e6, 'm--', 'f_{c} + 2 MHz', 'LineWidth', 1.5);
hold off;

% 2. FILTR DOCELOWY (96 MHz ±100 kHz)
disp('Projektowanie filtru docelowego (96 MHz ±100 kHz)');

% Parametry filtru docelowego
f_pass1 = fc - 100e3; % Dolna granica pasma przepustowego [Hz]
f_pass2 = fc + 100e3; % Górna granica pasma przepustowego [Hz]
f_stop1 = fc - 200e3; % Dolna granica pasma zaporowego [Hz]
f_stop2 = fc + 200e3; % Górna granica pasma zaporowego [Hz]
Rp = 3;               % Maksymalne zafalowanie w paśmie przepustowym [dB]
Rs = 40;              % Minimalne tłumienie w paśmie zaporowym [dB]

% Normalizacja częstotliwości względem częstotliwości Nyquista
f_nyquist = fc + 500e3; % Przyjęta częstotliwość Nyquista [Hz]
Wp = [f_pass1 f_pass2] / f_nyquist;
Ws = [f_stop1 f_stop2] / f_nyquist;

% Projektowanie filtru Butterwortha
[n, Wn] = buttord(Wp, Ws, Rp, Rs);
[b, a] = butter(n, Wn, 'bandpass');

% Wyświetlenie parametrów filtru docelowego
disp(['Rząd filtru docelowego: ' num2str(n)]);
disp(['Znormalizowane częstotliwości graniczne: ' num2str(Wn)]);

% Charakterystyka częstotliwościowa filtru docelowego
f = linspace(fc-500e3, fc+500e3, 10000);
w = 2*pi*f;
[h, ~] = freqs(b, a, w);
mag = 20*log10(abs(h));

% Rysowanie charakterystyki filtru docelowego
figure('Name', 'Filtr docelowy (96 MHz ±100 kHz)', 'Position', [100, 100, 900, 600]);
plot(f/1e6, mag, 'LineWidth', 1.5);
grid on;
xlabel('Częstotliwość [MHz]');
ylabel
::contentReference[oaicite:9]{index=9}
 