clc; close all; clear;

%% Zera i bieguny
z = [1j*5, -1j*5, 1j*15, -1j*15]; % zera: ±j5, ±j15
p = [-0.5+1j*9.5, -0.5-1j*9.5, ...
     -1+1j*10, -1-1j*10, ...
     -0.5+1j*10.5, -0.5-1j*10.5]; % bieguny

% Licznik i mianownik transmitancji
b = poly(z);
a = poly(p);

% Zakres częstotliwości
omega = linspace(0, 30, 1000); % od 0.1 do 100 rad/s
s = 1j * omega;

% Obliczenie transmitancji H(jω)
H = polyval(b, s) ./ polyval(a, s);
amplitude = abs(H);
phase = angle(H);
amplitude_dB = 20 * log10(amplitude);

% Normalizacja: wzmocnienie max = 1
[max_gain, idx_max] = max(amplitude);
b = b / max_gain;
H = polyval(b, s) ./ polyval(a, s);
amplitude = abs(H);
amplitude_dB = 20 * log10(amplitude);
phase = angle(H) * 180 / pi;

%% 1. Zera i bieguny
figure;
plot(real(z), imag(z), 'ro', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(real(p), imag(p), 'b*', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Część rzeczywista');
ylabel('Część urojona');
title('Zera (o) i Bieguny (*) na płaszczyźnie zespolonej');
grid on;
axis equal;
legend('Zera', 'Bieguny');

%% 2. Charakterystyka amplitudowa - linia
figure;
plot(omega, amplitude, 'k');
xlabel('\omega [rad/s]');
ylabel('|H(j\omega)|');
title('Charakterystyka amplitudowa (skala liniowa)');
grid on;

%% 3. Charakterystyka amplitudowa - decybelowa
figure;
semilogx(omega, amplitude_dB, 'm');
xlabel('\omega [rad/s]');
ylabel('20log_{10}|H(j\omega)| [dB]');
title('Charakterystyka amplitudowa (skala dB)');
grid on;

%% 4. Charakterystyka fazowa
figure;
semilogx(omega, phase, 'b');
xlabel('\omega [rad/s]');
ylabel('Faza [°]');
title('Charakterystyka fazowo-częstotliwościowa');
grid on;

%% 5. Analiza tłumienia
[min_amp, idx_min] = min(amplitude);
fprintf('Maksymalne wzmocnienie (po normalizacji): %.3f przy ω = %.2f rad/s\n', ...
        max(amplitude), omega(idx_max));
fprintf('Minimalne tłumienie: %.2f dB przy ω = %.2f rad/s\n', ...
        amplitude_dB(idx_min), omega(idx_min));