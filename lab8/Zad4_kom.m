%% Zadanie 4 - Filtr różniczkujący, dekodowanie FM
% Ten program implementuje trzy różne metody dekodowania sygnału FM
% według opisu z zadania 4 z laboratorium 8.
clear all;
close all;
clc;

%% Wczytanie sygnału FM
load('lab08_fm.mat'); % Wczytanie pliku z sygnałem modulowanym FM
% Plik powinien zawierać sygnał 'x' zawierający próbki sygnału FM

% Parametry sygnału
fs = 2e6;        % Częstotliwość próbkowania sygnału radiowego [Hz]
fc = 200e3;      % Częstotliwość nośna [Hz]
fs_audio = 8000; % Częstotliwość próbkowania sygnału audio [Hz]

% Czas trwania sygnału
t = (0:length(x)-1)/fs; % Wektor czasu dla sygnału FM

%% SPOSÓB 1: Demodulacja przez obliczenie składowych I(n) oraz Q(n)
% Metoda ta używa demodulacji kwadraturowej, obliczając składowe I i Q sygnału,
% a następnie fazę chwilową, której różnica odpowiada sygnałowi zdemodulowanemu.

% Demodulacja kwadraturowa - przemnożenie przez sygnały nośne przesuniete w fazie o 90 stopni
I = lowpass(x .* cos(2*pi*fc*t'), fc+5e3, fs); % Składowa I (w fazie)
Q = lowpass(-x .* sin(2*pi*fc*t'), fc+5e3, fs); % Składowa Q (w kwadraturze)

% Obliczenie sygnału zespolonego z składowych I i Q
y = I + 1j*Q; % Utworzenie sygnału analitycznego

% Obliczenie fazy chwilowej zgodnie ze wzorem (4) i (5) z treści zadania
y_shifted = [y(2:end); 0]; % y(n-1) - przesunięcie o jedną próbkę w prawo
phase_diff = angle(y .* conj(y_shifted)); % Różnica faz między próbkami
phase_diff(end) = phase_diff(end-1); % Uzupełnienie ostatniej wartości

% Skalowanie i normalizacja sygnału zdemodulowanego
x_demod_method1 = phase_diff;
x_demod_method1 = x_demod_method1 - mean(x_demod_method1); % Usunięcie składowej stałej
x_demod_method1 = x_demod_method1 / max(abs(x_demod_method1)); % Normalizacja

% Decymacja do częstotliwości próbkowania audio
decim_factor = round(fs / fs_audio); % Współczynnik decymacji
x_demod_method1_downsampled = x_demod_method1(1:decim_factor:end); % Decymacja

%% SPOSÓB 2: Demodulacja przez różniczkowanie i obwiednię
% Metoda ta wykorzystuje fakt, że po zróżniczkowaniu sygnału FM, 
% obwiednia otrzymanego sygnału zawiera sygnał modulujący.

% Projektowanie filtru różniczkującego w całym paśmie
Ndiff = 101; % Długość odpowiedzi impulsowej filtru różniczkującego
n_diff = -(Ndiff-1)/2:(Ndiff-1)/2; % Środkowanie filtru wokół zera

% Idealna odpowiedź impulsowa filtru różniczkującego
h_diff_ideal = zeros(size(n_diff));
for i = 1:length(n_diff)
    if n_diff(i) == 0
        h_diff_ideal(i) = 0; % Wartość w zerze
    else
        % Wzór na idealną odpowiedź impulsową filtru różniczkującego
        h_diff_ideal(i) = (sin(pi*n_diff(i)) - pi*n_diff(i)*cos(pi*n_diff(i))) / (pi*n_diff(i)^2);
    end
end

% Zastosowanie okna Blackmana do poprawy charakterystyki filtru różniczkującego
window_diff = blackman(Ndiff)'; % Okno Blackmana
h_diff = h_diff_ideal .* window_diff; % Zastosowanie okna

% Projektowanie filtru pasmowo-przepustowego
Nbp = 201; % Długość odpowiedzi impulsowej filtru pasmowo-przepustowego
f_low = (fc - 5e3) / (fs/2);  % Znormalizowana dolna częstotliwość odcięcia
f_high = (fc + 5e3) / (fs/2); % Znormalizowana górna częstotliwość odcięcia
h_bp = fir1(Nbp-1, [f_low, f_high], 'bandpass', blackman(Nbp)); % Filtr BP z oknem Blackmana

% Obliczenie odpowiedzi impulsowej kaskady filtrów różniczkującego i pasmowo-przepustowego
h_diff_bp = conv(h_diff, h_bp); % Splot odpowiedzi impulsowych obu filtrów

% Sprawdzenie charakterystyki częstotliwościowej filtru
[H_diff, w_diff] = freqz(h_diff, 1, 1024); % Charakterystyka filtru różniczkującego
[H_bp, w_bp] = freqz(h_bp, 1, 1024); % Charakterystyka filtru pasmowo-przepustowego
[H_diff_bp, w_diff_bp] = freqz(h_diff_bp, 1, 1024); % Charakterystyka kaskady filtrów

% Wizualizacja charakterystyk filtrów
figure;
subplot(3,1,1);
plot(w_diff/pi*fs/2/1000, abs(H_diff));
title('Charakterystyka amplitudowa filtru różniczkującego');
xlabel('Częstotliwość [kHz]');
ylabel('Amplituda');
grid on;

subplot(3,1,2);
plot(w_bp/pi*fs/2/1000, abs(H_bp));
title('Charakterystyka amplitudowa filtru pasmowo-przepustowego');
xlabel('Częstotliwość [kHz]');
ylabel('Amplituda');
grid on;
xlim([180, 220]); % Ograniczenie zakresu dla lepszej widoczności

subplot(3,1,3);
plot(w_diff_bp/pi*fs/2/1000, abs(H_diff_bp));
title('Charakterystyka amplitudowa kaskady filtrów');
xlabel('Częstotliwość [kHz]');
ylabel('Amplituda');
grid on;
xlim([180, 220]); % Ograniczenie zakresu dla lepszej widoczności

% Filtracja sygnału FM za pomocą kaskady filtrów
x_filtered = filter(h_diff_bp, 1, x); % Różniczkowanie i filtracja pasmowa w jednym kroku

% Odzyskanie obwiedni - kwadraturowanie (podniesienie do kwadratu)
x_squared = x_filtered.^2;

% Filtracja dolnoprzepustowa w celu wyodrębnienia obwiedni
cutoff_freq = 10e3 / (fs/2); % Częstotliwość odcięcia dla filtru LP (10 kHz)
h_lp = fir1(101, cutoff_freq, 'low', blackman(102)); % Filtr dolnoprzepustowy
x_env = filter(h_lp, 1, x_squared); % Filtracja dolnoprzepustowa

% Pierwiastkowanie dla odzyskania obwiedni
x_demod_method2 = sqrt(abs(x_env));

% Normalizacja sygnału
x_demod_method2 = x_demod_method2 - mean(x_demod_method2); % Usunięcie składowej stałej
x_demod_method2 = x_demod_method2 / max(abs(x_demod_method2)); % Normalizacja

% Decymacja do częstotliwości próbkowania audio
x_demod_method2_downsampled = x_demod_method2(1:decim_factor:end);

%% SPOSÓB 3: Demodulacja przez kaskadę oddzielnych filtrów BP i DIFF
% W tej metodzie zamiast stosować jeden filtr różniczkujący pasmowo-przepustowy,
% stosujemy dwa oddzielne filtry: pasmowo-przepustowy i różniczkujący.

% Najpierw filtracja pasmowo-przepustowa
x_bp_filtered = filter(h_bp, 1, x); % Wyodrębnienie sygnału w paśmie częstotliwości nośnej

% Filtr różniczkujący - prosty filtr IIR (b=[-1,1], a=1)
% Implementuje operację różnicy y(n) = x(n) - x(n-1)
b_diff_simple = [-1, 1]; % Współczynniki filtru różniczkującego
a_diff_simple = 1;

% Różniczkowanie przefiltrowanego sygnału
x_diff_filtered = filter(b_diff_simple, a_diff_simple, x_bp_filtered);

% Kwadraturowanie do odzyskania obwiedni
x_squared_method3 = x_diff_filtered.^2;

% Filtracja dolnoprzepustowa
x_env_method3 = filter(h_lp, 1, x_squared_method3);

% Pierwiastkowanie dla odzyskania obwiedni
x_demod_method3 = sqrt(abs(x_env_method3));

% Normalizacja sygnału
x_demod_method3 = x_demod_method3 - mean(x_demod_method3); % Usunięcie składowej stałej
x_demod_method3 = x_demod_method3 / max(abs(x_demod_method3)); % Normalizacja

% Decymacja do częstotliwości próbkowania audio
x_demod_method3_downsampled = x_demod_method3(1:decim_factor:end);

%% Porównanie wyników demodulacji
figure;
subplot(3,1,1);
plot((0:length(x_demod_method1_downsampled)-1)/fs_audio, x_demod_method1_downsampled);
title('Metoda 1: Demodulacja przez obliczenie składowych I(n) oraz Q(n)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(3,1,2);
plot((0:length(x_demod_method2_downsampled)-1)/fs_audio, x_demod_method2_downsampled);
title('Metoda 2: Demodulacja przez różniczkowanie i obwiednię');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(3,1,3);
plot((0:length(x_demod_method3_downsampled)-1)/fs_audio, x_demod_method3_downsampled);
title('Metoda 3: Demodulacja przez kaskadę filtrów BP i DIFF');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

%% Odsłuch zdemodulowanych sygnałów
% Każda z trzech metod jest odtwarzana po kolei
soundsc(x_demod_method1_downsampled, fs_audio);
pause(length(x_demod_method1_downsampled)/fs_audio + 0.5);
soundsc(x_demod_method2_downsampled, fs_audio);
pause(length(x_demod_method2_downsampled)/fs_audio + 0.5);
soundsc(x_demod_method3_downsampled, fs_audio);

%% Bonus: Implementacja filtru różniczkującego pasmowo-przepustowego za pomocą firls
% Metoda ta umożliwia zaprojektowanie filtru o złożonej charakterystyce bezpośrednio,
% bez konieczności kaskadowego łączenia dwóch oddzielnych filtrów.

% Projektowanie filtru różniczkującego pasmowo-przepustowego
Nbp_diff = 201; % Długość odpowiedzi impulsowej
f = [0, f_low-0.01, f_low, f_high, f_high+0.01, 1]; % Punkty częstotliwości
a = [0, 0, 1j, -1j, 0, 0]; % Odpowiedź różniczkująca (mnożenie przez jω)

% Projektowanie filtru za pomocą firls - metoda najmniejszych kwadratów dla filtru FIR
h_bp_diff = firls(Nbp_diff-1, f, abs(a));

% Konwersja odpowiedzi do rzeczywistej (usunięcie części urojonej)
h_bp_diff = real(h_bp_diff);

% Sprawdzenie charakterystyki filtru
[H_bp_diff, w_bp_diff] = freqz(h_bp_diff, 1, 1024);

figure;
plot(w_bp_diff/pi*fs/2/1000, abs(H_bp_diff));
title('Charakterystyka amplitudowa filtru różniczkującego pasmowo-przepustowego (firls)');
xlabel('Częstotliwość [kHz]');
ylabel('Amplituda');
grid on;
xlim([180, 220]);

% Filtracja sygnału FM przy użyciu filtru firls
x_filtered_firls = filter(h_bp_diff, 1, x);

% Kwadraturowanie
x_squared_firls = x_filtered_firls.^2;

% Filtracja dolnoprzepustowa
x_env_firls = filter(h_lp, 1, x_squared_firls);

% Pierwiastkowanie obwiedni
x_demod_firls = sqrt(abs(x_env_firls));

% Normalizacja
x_demod_firls = x_demod_firls - mean(x_demod_firls); % Usunięcie składowej stałej
x_demod_firls = x_demod_firls / max(abs(x_demod_firls)); % Normalizacja

% Decymacja do częstotliwości audio
x_demod_firls_downsampled = x_demod_firls(1:decim_factor:end);

% Porównanie wszystkich czterech metod
figure;
subplot(4,1,1);
plot((0:length(x_demod_method1_downsampled)-1)/fs_audio, x_demod_method1_downsampled);
title('Metoda 1: Demodulacja przez obliczenie składowych I(n) oraz Q(n)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(4,1,2);
plot((0:length(x_demod_method2_downsampled)-1)/fs_audio, x_demod_method2_downsampled);
title('Metoda 2: Demodulacja przez różniczkowanie i obwiednię');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(4,1,3);
plot((0:length(x_demod_method3_downsampled)-1)/fs_audio, x_demod_method3_downsampled);
title('Metoda 3: Demodulacja przez kaskadę filtrów BP i DIFF');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(4,1,4);
plot((0:length(x_demod_firls_downsampled)-1)/fs_audio, x_demod_firls_downsampled);
title('Bonus: Demodulacja z filtrem zaprojektowanym za pomocą firls');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;