clearvars; close all;
% Parametry początkowe
fs = 8000;  % Częstotliwość próbkowania (Hz)
T = 1;      % Czas trwania (s)
t = 0:1/fs:T-1/fs;  % Wektor czasu

% 1. Wygeneruj sygnał złożony z dwóch harmonicznych
A1 = -0.5; f1 = 34.2;  % Parametry pierwszej sinusoidy
A2 = 1;    f2 = 115.5; % Parametry drugiej sinusoidy
dref1 = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t);     

% 2. Sygnał SFM (Sinusoidalna modulacja częstotliwości)
fc = 1000;      % Częstotliwość nośna (Hz)
deltaf = 500;   % Dewiacja częstotliwości (Hz)
fm = 0.25;      % Częstotliwość modulująca (Hz)
dref2 = sin(2*pi*fc*t + (deltaf/fm) * sin(2*pi*fm*t));

% 3. Dźwięk silnika
[dref3, fs3] = audioread('silnik.wav'); 
dref3 = dref3(1:fs3);
dref3 = resample(dref3, fs, fs3);
dref3 = dref3.';

% Przechowywanie wszystkich sygnałów referencyjnych
drefs = {dref1, dref2, dref3};
titles = {'sinusoidalnych', 'SFM', 'silnika'}; 

% Wartości SNR do testowania
snr_values = [10, 20, 40];

% 2. Dobierz tak parametry filtru adaptacyjnego aby jednym zestawem parametrów
% odszumić jak najlepiej wszystkie 3 wersje sygnału
M = 32;     % Długość filtru - zwiększona dla lepszej adaptacji
mi = 0.008; % Współczynnik adaptacji - zmniejszony dla stabilności

% Przetwarzanie każdego sygnału referencyjnego
for i = 1:length(drefs)
    dref = drefs{i};
    
    % Tworzenie wykresu dla bieżącego sygnału
    figure;
    sgtitle(['Porównanie sygnałów ', titles{i}])
    
    % Przetwarzanie dla każdej wartości SNR
    for j = 1:length(snr_values)
        snr = snr_values(j);
        
        % Dodaj szum do sygnału referencyjnego
        d = awgn(dref, snr, 'measured');      % Zaszumiony sygnał
        x = [d(1), d(1:end-1)];               % Opóźniony sygnał wejściowy
        
        % Inicjalizacja zmiennych filtru adaptacyjnego
        y = zeros(size(d));    % Wyjście filtru
        e = zeros(size(d));    % Sygnał błędu
        bx = zeros(M, 1);      % Bufor wejściowy
        h = zeros(M, 1);       % Współczynniki filtru
        
        % Implementacja algorytmu LMS
        for n = 1:length(x)
            bx = [x(n); bx(1:M-1)];   % pobierz nową próbkę x[n] do bufora
            y(n) = h' * bx;           % oblicz y[n] = sum(h .* bx) – filtr FIR
            e(n) = d(n) - y(n);       % oblicz e[n]
            h = h + mi * e(n) * bx;   % LMS
            % h = h + mi * e(n) * bx / (bx'*bx + 1e-10); % NLMS
        end
        
        % Obliczenie SNR zgodnie z definicją w zadaniu
        SNR_db = 10 * log10(sum(dref.^2) / sum((dref - y).^2));
        fprintf('Sygnał: %s, SNR wejściowy: %d dB, SNR po odszumianiu: %.2f dB\n', ...
            titles{i}, snr, SNR_db);
        
        % Wyświetlenie wyników
        subplot(3, 1, j);
        plot(t, dref, t, d, t, y);
        legend('referencyjny', 'odniesienia', 'rezultat')
        xlabel('Czas [s]');
        ylabel('Amplituda');
        title(sprintf('SNR = %d dB, Wynik SNR = %.2f dB', snr, SNR_db))
    end
    pause;
end

% Optional: Experiment with NLMS (Normalized LMS) algorithm
% Odkomentuj poniższy kod, aby uruchomić wersję NLMS i porównać wyniki

% % Przykład tylko dla sygnału sinusoidalnego przy 20dB
% dref = drefs{1};
% snr = 20;
% d = awgn(dref, snr, 'measured');
% x = [d(1), d(1:end-1)];
% 
% % Inicjalizacja zmiennych
% y_nlms = zeros(size(d));
% e_nlms = zeros(size(d));
% bx = zeros(M, 1);
% h = zeros(M, 1);
% 
% % Algorytm NLMS
% for n = 1:length(x)
%     bx = [x(n); bx(1:M-1)];
%     y_nlms(n) = h' * bx;
%     e_nlms(n) = d(n) - y_nlms(n);
%     
%     % Aktualizacja NLMS z normalizacją i małym epsilon, aby zapobiec dzieleniu przez zero
%     norm_factor = bx' * bx + 1e-10;
%     h = h + mi * (e_nlms(n) / norm_factor) * bx;
% end