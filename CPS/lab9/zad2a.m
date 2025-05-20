%% Zadanie 2a: Identyfikacja obiektu

clear;
close all;

%% Wczytanie próbki mowy
[x, fs] = audioread('mowa8000.wav');

%% Tworzenie rzeczywistej odpowiedzi impulsowej h_rzecz
N = 256;
h_rzecz = zeros(N, 1);
h_rzecz(31) = 0.1;    % i=30
h_rzecz(121) = -0.5;  % i=120
h_rzecz(256) = 0.8;   % i=255

%% Filtrowanie sygnału mowy zadanym filtrem
d = filter(h_rzecz, 1, x);

%% Implementacja filtru adaptacyjnego - identyfikacja obiektu
M = 256;             % Długość filtru (równa rzeczywistej odpowiedzi impulsowej)
mi = 0.01;           % Współczynnik szybkości adaptacji
y = zeros(size(x));  % Wyjście filtru
e = zeros(size(x));  % Sygnał błędu
bx = zeros(M, 1);    % Bufor na próbki wejściowe
h = zeros(M, 1);     % Początkowe wagi filtru

% Główna pętla adaptacji
for n = 1:length(x)
    bx = [x(n); bx(1:M-1)];  % Pobierz nową próbkę x[n] do bufora
    y(n) = h' * bx;          % Oblicz y[n] = sum(h .* bx) – filtr FIR
    e(n) = d(n) - y(n);      % Oblicz błąd e[n]
    h = h + mi * e(n) * bx;  % Aktualizacja wag (LMS)
end

%% Analiza wyników
figure('Name', 'Identyfikacja obiektu - sygnał mowy', 'Position', [100, 100, 1000, 800]);

% Wykres odpowiedzi impulsowych
subplot(2, 1, 1);
stem(0:N-1, h_rzecz, 'b', 'LineWidth', 1.5);
hold on;
stem(0:M-1, h, 'r:', 'LineWidth', 1);
hold off;
grid on;
title('Porównanie odpowiedzi impulsowych');
xlabel('Próbka');
ylabel('Amplituda');
legend('Rzeczywista odpowiedź impulsowa', 'Estymowana odpowiedź impulsowa');

% Wykres błędu adaptacji
subplot(2, 1, 2);
plot(10*log10(e.^2));
grid on;
title('Błąd kwadratowy (dB)');
xlabel('Próbka');
ylabel('Błąd [dB]');

%% Test z szumem białym
disp('Test z szumem białym:');

% Generowanie szumu białego
len = length(x);
white_noise = randn(len, 1);

% Filtrowanie szumu białego
d_noise = filter(h_rzecz, 1, white_noise);

% Adaptacja z szumem białym
y_noise = zeros(size(white_noise));
e_noise = zeros(size(white_noise));
bx_noise = zeros(M, 1);
h_noise = zeros(M, 1);

% Główna pętla adaptacji dla szumu białego
for n = 1:length(white_noise)
    bx_noise = [white_noise(n); bx_noise(1:M-1)];
    y_noise(n) = h_noise' * bx_noise;
    e_noise(n) = d_noise(n) - y_noise(n);
    h_noise = h_noise + mi * e_noise(n) * bx_noise;
end

%% Analiza wyników dla szumu białego
figure('Name', 'Identyfikacja obiektu - szum biały', 'Position', [100, 100, 1000, 800]);

% Wykres odpowiedzi impulsowych
subplot(2, 1, 1);
stem(0:N-1, h_rzecz, 'b', 'LineWidth', 1.5);
hold on;
stem(0:M-1, h_noise, 'r:', 'LineWidth', 1);
hold off;
grid on;
title('Porównanie odpowiedzi impulsowych (szum biały)');
xlabel('Próbka');
ylabel('Amplituda');
legend('Rzeczywista odpowiedź impulsowa', 'Estymowana odpowiedź impulsowa (szum biały)');

% Wykres błędu adaptacji
subplot(2, 1, 2);
plot(10*log10(e_noise.^2));
grid on;
title('Błąd kwadratowy (dB) - szum biały');
xlabel('Próbka');
ylabel('Błąd [dB]');

% Porównanie dokładności estymacji
mse_speech = mean((h - h_rzecz).^2);
mse_noise = mean((h_noise - h_rzecz).^2);

disp(['MSE dla estymacji z sygnału mowy: ', num2str(mse_speech)]);
disp(['MSE dla estymacji z szumu białego: ', num2str(mse_noise)]);

% Porównanie estymacji dla różnych współczynników adaptacji (opcjonalnie)
mi_values = [0.001, 0.01, 0.05, 0.1];
figure('Name', 'Wpływ współczynnika adaptacji', 'Position', [100, 100, 1200, 800]);

for i = 1:length(mi_values)
    mi_test = mi_values(i);
    h_test = zeros(M, 1);
    bx_test = zeros(M, 1);
    
    for n = 1:length(white_noise)
        bx_test = [white_noise(n); bx_test(1:M-1)];
        y_test = h_test' * bx_test;
        e_test = d_noise(n) - y_test;
        h_test = h_test + mi_test * e_test * bx_test;
    end
    
    subplot(2, 2, i);
    stem(0:N-1, h_rzecz, 'b', 'LineWidth', 1.5);
    hold on;
    stem(0:M-1, h_test, 'r:', 'LineWidth', 1);
    hold off;
    grid on;
    title(['Współczynnik adaptacji mi = ', num2str(mi_test)]);
    xlabel('Próbka');
    ylabel('Amplituda');
    legend('Rzeczywista', 'Estymowana');
end