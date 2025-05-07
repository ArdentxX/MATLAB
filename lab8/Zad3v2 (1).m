%% Zadanie 3 - Filtr interpolatora i decymatora cyfrowego, mikser audio
clear all;
close all;
clc;

%% Definicja parametrów sygnałów
% x1: f1=1001.2 Hz, fs1= 8000 Hz, t=1 s
% x2: f2= 303.1 Hz, fs2=32000 Hz, t=1 s
% x3: f3=2110.4 Hz, fs3=48000 Hz, t=1 s

% Parametry sygnału wynikowego
fs_target = 48000; % Docelowa częstotliwość próbkowania
t = 1; % Czas trwania sygnału w sekundach

% Parametry sygnałów składowych
f1 = 1001.2; % Hz
fs1 = 8000; % Hz
f2 = 303.1; % Hz
fs2 = 32000; % Hz
f3 = 2110.4; % Hz
fs3 = 48000; % Hz

%% Generacja sygnałów sinusoidalnych
t1 = (0:1/fs1:t-1/fs1)';
t2 = (0:1/fs2:t-1/fs2)';
t3 = (0:1/fs3:t-1/fs3)';
t_target = (0:1/fs_target:t-1/fs_target)';

x1 = sin(2*pi*f1*t1);
x2 = sin(2*pi*f2*t2);
x3 = sin(2*pi*f3*t3);

% Generacja oczekiwanego sygnału wynikowego (analitycznie)
x4_expected = sin(2*pi*f1*t_target) + sin(2*pi*f2*t_target) + sin(2*pi*f3*t_target);

%% Repróbkowanie sygnału x1 z fs1=8000 Hz do fs_target=48000 Hz
% Współczynnik nadpróbkowania (upsampling)
up_factor1 = fs_target / fs1; % = 6

% Nadpróbkowanie przez wstawienie zer (upsampling)
x1_upsampled = zeros(length(x1) * up_factor1, 1);
x1_upsampled(1:up_factor1:end) = x1;

% Projektowanie filtru interpolującego dla x1
% Częstotliwość odcięcia = fs1/2 / fs_target = 8000/2 / 48000 = 1/12
cutoff1 = (fs1/2) / fs_target;
filter_order1 = 48; % Rząd filtru interpolującego 
b_interp1 = fir1(filter_order1, cutoff1, 'low') * up_factor1;

% Filtracja interpolująca
x1_interpolated = filter(b_interp1, 1, x1_upsampled);

% Kompensacja opóźnienia filtru
delay1 = filter_order1 / 2;
x1_interpolated = x1_interpolated(delay1+1:end);
x1_interpolated = [x1_interpolated; zeros(delay1, 1)];

%% Repróbkowanie sygnału x2 z fs2=32000 Hz do fs_target=48000 Hz
% Współczynnik nadpróbkowania (upsampling)
up_factor2 = fs_target / fs2; % = 3/2 = 1.5

% Implementacja ułamkowego repróbkowania (upsampling o czynnik 3, potem decymacja o czynnik 2)
% Najpierw nadpróbkowanie o czynnik 3
x2_upsampled_3 = zeros(length(x2) * 3, 1);
x2_upsampled_3(1:3:end) = x2;

% Projektowanie filtru interpolującego dla x2
cutoff2_up = (fs2/2) / (fs2*3);
filter_order2_up = 36;
b_interp2_up = fir1(filter_order2_up, cutoff2_up, 'low') * 3;

% Filtracja interpolująca
x2_interpolated_3 = filter(b_interp2_up, 1, x2_upsampled_3);

% Kompensacja opóźnienia filtru
delay2_up = filter_order2_up / 2;
x2_interpolated_3 = x2_interpolated_3(delay2_up+1:end);
x2_interpolated_3 = [x2_interpolated_3; zeros(delay2_up, 1)];

% Następnie decymacja o czynnik 2
% Projektowanie filtru decymującego (antyaliasingowego)
cutoff2_down = (fs2*3/2) / (fs2*3);
filter_order2_down = 24;
b_decim2 = fir1(filter_order2_down, cutoff2_down, 'low');

% Filtracja antyaliasingowa przed decymacją
x2_filtered = filter(b_decim2, 1, x2_interpolated_3);

% Kompensacja opóźnienia filtru
delay2_down = filter_order2_down / 2;
x2_filtered = x2_filtered(delay2_down+1:end);
x2_filtered = [x2_filtered; zeros(delay2_down, 1)];

% Decymacja
x2_interpolated = x2_filtered(1:2:end);

%% Sygnał x3 już ma docelową częstotliwość próbkowania fs3=48000 Hz
x3_interpolated = x3;

%% Łączenie sygnałów (przez dodawanie)
% Dopasowanie długości sygnałów
min_length = min([length(x1_interpolated), length(x2_interpolated), length(x3_interpolated)]);
x1_interpolated = x1_interpolated(1:min_length);
x2_interpolated = x2_interpolated(1:min_length);
x3_interpolated = x3_interpolated(1:min_length);

% Dodawanie sygnałów
x4 = x1_interpolated + x2_interpolated + x3_interpolated;

%% Analiza wyników
% Porównanie z oczekiwanym sygnałem
x4_expected = x4_expected(1:min_length);

% Obliczenie błędu średniokwadratowego (MSE)
mse = mean((x4 - x4_expected).^2);
fprintf('Błąd średniokwadratowy między rzeczywistym a oczekiwanym sygnałem: %e\n', mse);

% Wizualizacja wyników
t_plot = t_target(1:min_length);

figure;
subplot(5,1,1);
plot(t_plot(1:500), x1_interpolated(1:500));
title('Sygnał x1 po repróbkowaniu do 48 kHz');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(5,1,2);
plot(t_plot(1:500), x2_interpolated(1:500));
title('Sygnał x2 po repróbkowaniu do 48 kHz');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(5,1,3);
plot(t_plot(1:500), x3_interpolated(1:500));
title('Sygnał x3 (już w 48 kHz)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(5,1,4);
plot(t_plot(1:500), x4(1:500));
title('Sygnał wynikowy x4 = x1 + x2 + x3');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(5,1,5);
plot(t_plot(1:500), x4_expected(1:500));
title('Oczekiwany sygnał x4 (wygenerowany analitycznie)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

% Analiza widmowa sygnału wynikowego
NFFT = 2^nextpow2(length(x4));
X4 = fft(x4, NFFT) / length(x4);
f = fs_target/2 * linspace(0, 1, NFFT/2+1);

figure;
plot(f, 2*abs(X4(1:NFFT/2+1)));
title('Widmo sygnału wynikowego x4');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda');
grid on;
xlim([0, 3000]); % Ograniczenie zakresu częstotliwości do 3 kHz dla lepszej widoczności

%% Odsłuchanie
sound(x4/max(abs(x4)), fs_target);

%% Miksowanie rzeczywistych sygnałów audio
% Wczytanie plików audio
[x1_real, fs1_real] = audioread('x1.wav');
[x2_real, fs2_real] = audioread('x2.wav');

% Repróbkowanie do 48 kHz

% Repróbkowanie x1_real
% Podobnie jak wcześniej, ale dostosowujemy współczynniki do rzeczywistych częstotliwości
up_factor1_real = fs_target / fs1_real;

% Implementacja repróbkowania dla x1_real
if up_factor1_real == round(up_factor1_real) % Jeśli współczynnik jest liczbą całkowitą
    x1_real_upsampled = zeros(length(x1_real) * up_factor1_real, 1);
    x1_real_upsampled(1:up_factor1_real:end) = x1_real;
    
    cutoff1_real = (fs1_real/2) / fs_target;
    filter_order1_real = 48;
    b_interp1_real = fir1(filter_order1_real, cutoff1_real, 'low') * up_factor1_real;
    
    x1_real_interpolated = filter(b_interp1_real, 1, x1_real_upsampled);
    
    delay1_real = filter_order1_real / 2;
    x1_real_interpolated = x1_real_interpolated(delay1_real+1:end);
    x1_real_interpolated = [x1_real_interpolated; zeros(delay1_real, 1)];
else
    % Dla ułamkowego współczynnika repróbkowania używamy resample
    x1_real_interpolated = resample(x1_real, fs_target, fs1_real);
end

% Repróbkowanie x2_real
up_factor2_real = fs_target / fs2_real;

if up_factor2_real == round(up_factor2_real) % Jeśli współczynnik jest liczbą całkowitą
    x2_real_upsampled = zeros(length(x2_real) * up_factor2_real, 1);
    x2_real_upsampled(1:up_factor2_real:end) = x2_real;
    
    cutoff2_real = (fs2_real/2) / fs_target;
    filter_order2_real = 48;
    b_interp2_real = fir1(filter_order2_real, cutoff2_real, 'low') * up_factor2_real;
    
    x2_real_interpolated = filter(b_interp2_real, 1, x2_real_upsampled);
    
    delay2_real = filter_order2_real / 2;
    x2_real_interpolated = x2_real_interpolated(delay2_real+1:end);
    x2_real_interpolated = [x2_real_interpolated; zeros(delay2_real, 1)];
else
    % Dla ułamkowego współczynnika repróbkowania używamy resample
    x2_real_interpolated = resample(x2_real, fs_target, fs2_real);
end

% Upewnij się, że są wektorami kolumnowymi
x1_real_interpolated = x1_real_interpolated(:);
x2_real_interpolated = x2_real_interpolated(:);

% Miksowanie rzeczywistych sygnałów audio
min_length_real = min(length(x1_real_interpolated), length(x2_real_interpolated));
x1_real_interpolated = x1_real_interpolated(1:min_length_real);
x2_real_interpolated = x2_real_interpolated(1:min_length_real);


% Dodawanie
mixed_real = x1_real_interpolated + x2_real_interpolated;


% Normalizacja i zapis do pliku
mixed_real = mixed_real / max(abs(mixed_real));
audiowrite('mixed_real_48kHz.wav', mixed_real, fs_target);

%% CZĘŚĆ OPCJONALNA: Miksowanie do częstotliwości CD-Audio (44.1 kHz)
fs_cd = 44100; % Częstotliwość CD-Audio

% Repróbkowanie zmiksowanych sygnałów do 44.1 kHz
mixed_real_cd = resample(mixed_real, fs_cd, fs_target);

% Normalizacja i zapis do pliku
mixed_real_cd = mixed_real_cd / max(abs(mixed_real_cd));
sound(mixed_real_cd, fs_cd);

%% CZĘŚĆ OPCJONALNA: Miksowanie sygnałów syntetycznych różnymi metodami

% Metoda 1: Interpolacja liniowa
x1_linear = interp1(t1, x1, t_target, 'linear');
x2_linear = interp1(t2, x2, t_target, 'linear');
x3_linear = x3; % już ma właściwą częstotliwość próbkowania

% Dopasowanie długości
min_length_linear = min([length(x1_linear), length(x2_linear), length(x3_linear)]);
x1_linear = x1_linear(1:min_length_linear);
x2_linear = x2_linear(1:min_length_linear);
x3_linear = x3_linear(1:min_length_linear);

% Miksowanie
x4_linear = x1_linear + x2_linear + x3_linear;

% Metoda 2: Rekonstrukcja sygnału metodą splotu z sin(x)/x
% Funkcja sinc
sinc_func = @(x) sin(pi*x)./(pi*x + eps);

% Repróbkowanie x1 metodą sin(x)/x
x1_sinc = zeros(size(t_target));
for i = 1:length(t1)
    % Dla każdej próbki oryginalnego sygnału
    t_diff = t_target - t1(i);
    x1_sinc = x1_sinc + x1(i) * sinc_func(fs1 * t_diff);
end

% Repróbkowanie x2 metodą sin(x)/x
x2_sinc = zeros(size(t_target));
for i = 1:length(t2)
    % Dla każdej próbki oryginalnego sygnału
    t_diff = t_target - t2(i);
    x2_sinc = x2_sinc + x2(i) * sinc_func(fs2 * t_diff);
end

x3_sinc = x3; % już ma właściwą częstotliwość próbkowania

% Dopasowanie długości
min_length_sinc = min([length(x1_sinc), length(x2_sinc), length(x3_sinc)]);
x1_sinc = x1_sinc(1:min_length_sinc);
x2_sinc = x2_sinc(1:min_length_sinc);
x3_sinc = x3_sinc(1:min_length_sinc);

% Miksowanie
x4_sinc = x1_sinc + x2_sinc + x3_sinc;

% Porównanie metod
figure;
subplot(4,1,1);
plot(t_plot(1:500), x4_expected(1:500));
title('Oczekiwany sygnał x4 (wygenerowany analitycznie)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(4,1,2);
plot(t_plot(1:500), x4(1:500));
title('Sygnał x4 - filtr FIR');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(4,1,3);
plot(t_plot(1:500), x4_linear(1:500));
title('Sygnał x4 - interpolacja liniowa');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

subplot(4,1,4);
plot(t_plot(1:500), x4_sinc(1:500));
title('Sygnał x4 - rekonstrukcja sinc');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

% Obliczenie błędów dla poszczególnych metod
mse_fir = mean((x4 - x4_expected(1:length(x4))).^2);
mse_linear = mean((x4_linear - x4_expected(1:length(x4_linear))).^2);
mse_sinc = mean((x4_sinc - x4_expected(1:length(x4_sinc))).^2);

fprintf('Błąd MSE dla metody FIR: %e\n', mse_fir);
fprintf('Błąd MSE dla interpolacji liniowej: %e\n', mse_linear);
fprintf('Błąd MSE dla rekonstrukcji sinc: %e\n', mse_sinc);

% Odsłuch
sound(x4/max(abs(x4)), fs_target);
sound(x4_linear/max(abs(x4_linear)), fs_target);
sound(x4_sinc/max(abs(x4_sinc)), fs_target);
