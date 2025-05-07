%% Zadanie 1 - Filtr Hilberta, demodulacja AM
clear all;
close all;
clc;

%% Wczytanie sygnału z pliku
load('lab08_am.mat');
% Załóżmy, że przedostatnia cyfra legitymacji to 3
signal_idx = 3; % Zmień na swoją przedostatnią cyfrę legitymacji
x = eval(['s', num2str(signal_idx)]);

% Parametry sygnału
fs = 1000; % Częstotliwość próbkowania [Hz]
t = 0:1/fs:(length(x)-1)/fs; % Wektor czasu
fc = 200; % Częstotliwość nośna [Hz]

%% Projektowanie filtru Hilberta (przesuwnik fazowy o -π/2)
N = 63; % Długość filtru (musi być nieparzysta)
n = -(N-1)/2:(N-1)/2; % Indeksy próbek
h_ideal = zeros(1, length(n));

% Idealna odpowiedź impulsowa filtru Hilberta
for i = 1:length(n)
    if n(i) ~= 0
        h_ideal(i) = 2/(pi*n(i)) * sin(pi*n(i)/2)^2;
    else
        h_ideal(i) = 0; % W n=0 odpowiedź impulsowa jest zerowa
    end
end

% Zastosowanie okna Hamminga do poprawy charakterystyki filtru
window = hamming(N)';
h = h_ideal .* window;

% Sprawdzenie charakterystyki częstotliwościowej filtru
[H, w] = freqz(h, 1, 1024);
figure;
subplot(2,1,1);
plot(w/pi, abs(H));
title('Charakterystyka amplitudowa filtru Hilberta');
xlabel('Częstotliwość znormalizowana (\times\pi rad/próbkę)');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(w/pi, unwrap(angle(H)));
title('Charakterystyka fazowa filtru Hilberta');
xlabel('Częstotliwość znormalizowana (\times\pi rad/próbkę)');
ylabel('Faza [rad]');
grid on;

%% Filtracja - obliczenie transformaty Hilberta sygnału
x_hilbert = filter(h, 1, x);

% Kompensacja opóźnienia filtru
delay = (N-1)/2;
x_hilbert = x_hilbert(delay+1:end);
x_sync = x(1:end-delay);

% Obliczenie obwiedni - demodulacja AM
envelope = sqrt(x_sync.^2 + x_hilbert.^2);

%% Analiza widmowa sygnału obwiedni - wyznaczenie parametrów m(t)
t_sync = t(1:length(envelope));
window_env = hann(length(envelope))';
envelope_windowed = envelope .* window_env;

Y = fft(envelope_windowed);
P2 = abs(Y/length(envelope));
P1 = P2(1:floor(length(envelope)/2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(length(envelope)/2))/length(envelope);

figure;
plot(f, P1);
title('Widmo amplitudowe obwiedni (sygnału modulującego)');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda');
grid on;

% Znalezienie największych składowych widma poza składową stałą
[pks, locs] = findpeaks(P1(2:end), 'SortStr', 'descend', 'NPeaks', 3);
f_peaks = f(locs+1);  % +1 ponieważ pominęliśmy składową stałą

% Wyznaczenie parametrów sygnału modulującego
A1 = pks(1) * 2;  % Amplituda 1 (x2 bo widmo jednostronne)
A2 = pks(2) * 2;  % Amplituda 2
A3 = pks(3) * 2;  % Amplituda 3
f1 = f_peaks(1);  % Częstotliwość 1
f2 = f_peaks(2);  % Częstotliwość 2
f3 = f_peaks(3);  % Częstotliwość 3
DC = P1(1);       % Składowa stała

fprintf('Parametry sygnału modulującego m(t):\n');
fprintf('Składowa stała: %.4f\n', DC);
fprintf('A1 = %.4f, f1 = %.2f Hz\n', A1, f1);
fprintf('A2 = %.4f, f2 = %.2f Hz\n', A2, f2);
fprintf('A3 = %.4f, f3 = %.2f Hz\n', A3, f3);

%% Wizualizacja wyników demodulacji
figure;
subplot(3,1,1);
plot(t, x);
title('Sygnał zmodulowany AM');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(3,1,2);
plot(t(1:length(x_hilbert)), x_sync, 'b');
hold on;
plot(t(1:length(x_hilbert)), x_hilbert, 'g');
title('Sygnał wejściowy i jego transformata Hilberta');
legend('Sygnał wejściowy', 'Transformata Hilberta');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(3,1,3);
plot(t_sync, envelope);
title('Obwiednia - sygnał zdemodulowany');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

%% CZĘŚĆ OPCJONALNA: Odtworzenie sygnału x
% Rekonstrukcja sygnału modulującego
t_mod = 0:1/fs:1-1/fs;
m_t = DC + A1*cos(2*pi*f1*t_mod) + A2*cos(2*pi*f2*t_mod) + A3*cos(2*pi*f3*t_mod);

% Rekonstrukcja sygnału zmodulowanego AM
x_reconstructed = m_t .* cos(2*pi*fc*t_mod);

% Porównanie zrekonstruowanego sygnału z oryginalnym
figure;
subplot(2,1,1);
plot(t, x, 'b');
hold on;
plot(t_mod, x_reconstructed, 'r--');
title('Porównanie oryginalnego sygnału AM z rekonstrukcją');
legend('Oryginalny', 'Zrekonstruowany');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(2,1,2);
plot(t, x - x_reconstructed', 'g');
title('Różnica między sygnałami');
xlabel('Czas [s]');
ylabel('Różnica');
grid on;

% Obliczenie błędu średniokwadratowego (MSE)
mse_error = mean((x - x_reconstructed').^2);
fprintf('Błąd średniokwadratowy (MSE): %.6f\n', mse_error);
