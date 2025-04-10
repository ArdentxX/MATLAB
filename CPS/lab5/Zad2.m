clearvars; clc; close all;

%% Parametry filtru
f_odc = 100;                   % częstotliwość odcięcia [Hz]
omega_odc = 2*pi*f_odc;        % częstotliwość odcięcia [rad/s]
rzedy = [2, 4, 6, 8];         % rzędy filtru do zaprojektowania

%% Projekt filtru Butterwortha
% Inicjalizacja komórkowej tablicy funkcji przenoszenia
filtry = cell(length(rzedy), 1);

for i = 1:length(rzedy)
    N = rzedy(i);
    
    % Wektoryzacja obliczenia biegunów:
    % Obliczamy kąty wg: π/2 + (2k-1)*π/(2N)
    idx = 1:N;
    bieguny = omega_odc * exp(1j * (pi/2 + (2*idx - 1)*pi/(2*N))).';
    
    % Transmitancja filtru Butterwortha:
    % Licznik to omega_odc^N (zera są w nieskończoności)
    licznik = omega_odc^N;
    % Mianownik jako współczynniki wielomianu o biegunach – wymuszamy rzeczywistość współczynników
    mianownik = real(poly(bieguny));
    
    % Funkcja przenoszenia
    filtry{i} = tf(licznik, mianownik);
end

%% Definicja zakresów częstotliwości
% Dla obliczeń opartych na freqs (podobnie jak w podanym przykładzie)
freq_range = 0:10:1000;   % [Hz]

%% Wykresy amplitudowe – dwa wykresy na jednej figurze
figure;
% Wykres amplitudowy – skala liniowa
subplot(2,1,1);
for i = 1:length(rzedy)
    % Używamy freqs na podstawie współczynników funkcji przenoszenia
    [H_odp, ~] = freqs(filtry{i}.Numerator{1}, filtry{i}.Denominator{1}, 2*pi*freq_range);
    mag = abs(H_odp);
    mag_dB = 20*log10(mag);
    plot(freq_range, mag_dB, 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(rzedy(i))]);
    hold on;
end
grid on;
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
title('Charakterystyka amplitudowa - skala liniowa');
legend show;

% Wykres amplitudowy – skala logarytmiczna
subplot(2,1,2);
for i = 1:length(rzedy)
    [H_odp, ~] = freqs(filtry{i}.Numerator{1}, filtry{i}.Denominator{1}, 2*pi*freq_range);
    mag = abs(H_odp);
    mag_dB = 20*log10(mag);
    semilogx(freq_range, mag_dB, 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(rzedy(i))]);
    hold on;
end
grid on;
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
title('Charakterystyka amplitudowa - skala logarytmiczna');
legend show;

%% Wykres charakterystyki fazowej – zgodnie z podanym przykładem
figure;
for i = 1:length(rzedy)
    [H_odp, ~] = freqs(filtry{i}.Numerator{1}, filtry{i}.Denominator{1}, 2*pi*freq_range);
    faza_deg = angle(H_odp) * (180/pi);  % faza w stopniach
    plot(freq_range, unwrap(faza_deg), 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(rzedy(i))]);
    hold on;
end
grid on;
xlabel('Częstotliwość [Hz]');
ylabel('Faza [°]');
title('Charakterystyka fazowa filtru Butterwortha');
legend show;

%% Odpowiedzi czasowe dla filtru 4. rzędu (N = 4)
figure;
% Impulsowa odpowiedź czasowa
subplot(2,1,1);
impulse(filtry{2});
title('Odpowiedź impulsowa (N = 4)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;

% Skokowa odpowiedź czasowa
subplot(2,1,2);
step(filtry{2});
title('Odpowiedź skokowa (N = 4)');
xlabel('Czas [s]');
ylabel('Amplituda');
grid on;
