clear all;
close all;

load("lab_03.mat");

numer_indeksu = 422852;
nazwa_wektora = sprintf('x_%d', mod(numer_indeksu, 16) + 1);
disp(mod(numer_indeksu, 16) + 1);

K = 8;  % Liczba ramek
N = 512; % Długość ramki (bez prefiksu)
M = 32;  % Długość prefiksu
ramki = zeros(N, K); % Macierz do przechowywania ramek

% Ekstrakcja ramek z sygnału (usuwając prefiks)
for m = 0:K-1
    start_idx = m * (N + M) + M + 1;
    ramki(:, m+1) = x_5(start_idx:start_idx + N - 1);
end

% Obliczenie DFT (FFT) dla każdej ramki
widmo = fft(ramki);

% Obliczenie wartości amplitudowych widma
amplitudy = abs(widmo);

% Znalezienie harmonicznych o najwyższej amplitudzie
[~, harmoniczne] = max(amplitudy, [], 1);

% Wyświetlenie wyników
disp('Najsilniejsze harmoniczne w każdej ramce:');
disp(harmoniczne);
