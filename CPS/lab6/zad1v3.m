v% ZADANIE 1 - Filtr IIR z poprawionym pre-warpingiem i zgodnością z poleceniem

% Wczytanie danych filtru analogowego
load('butter.mat');  % zawiera: z, p, k

% Parametry projektu
fs = 16000;           % częstotliwość próbkowania
T = 1/fs;
f_low_analog = 1189;  % dolna częstotliwość graniczna (analogowa)
f_high_analog = 1229; % górna częstotliwość graniczna (analogowa)

% Konwersja z-p-k do transmitancji H(s)
[b_a, a_a] = zp2tf(z, p, k);

% Konwersja H(s) -> H(z) z użyciem transformacji biliniowej
[b_d, a_d] = bilinear(b_a, a_a, fs);

% Rysowanie charakterystyk filtru analogowego i cyfrowego
figure;
[h_analog, f_analog] = freqs(b_a, a_a, 1024);
semilogx(f_analog/(2*pi), 20*log10(abs(h_analog)), 'LineWidth', 2, 'Color', [0 0.447 0.741]); hold on;
[h_digital, f_digital] = freqz(b_d, a_d, 1024, fs);
plot(f_digital, 20*log10(abs(h_digital)), 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);

xline([f_low_analog, f_high_analog], '--k', 'LineWidth', 1.2);
xlim([1 fs/2]); % Zakres częstotliwości od 1 Hz do Nyquista

legend({'$H(s)$', '$H(z)$', 'Granice analogowe'}, 'Interpreter', 'latex', 'Location', 'Southwest');
title('Charakterystyki filtru analogowego i cyfrowego', 'FontSize', 14);
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
grid on; grid minor;

% Zadanie opcjonalne - Pre-warping
% Korekcja pre-warping dla częstotliwości granicznych
w_low_analog = 2*pi*f_low_analog;
w_high_analog = 2*pi*f_high_analog;

% Przeliczenie częstotliwości analogowych z uwzględnieniem pre-warpingu
w_low_prewarped = 2/T * tan(w_low_analog * T / 2);
w_high_prewarped = 2/T * tan(w_high_analog * T / 2);

% Projekt nowego filtru analogowego z poprawionymi częstotliwościami
[n, Wn] = buttord([w_low_prewarped w_high_prewarped], [w_low_prewarped*0.9 w_high_prewarped*1.1], 3, 40, 's');
[z_new, p_new, k_new] = butter(n, Wn, 'bandpass', 's');

% Konwersja do transmitancji H(s)
[b_a_new, a_a_new] = zp2tf(z_new, p_new, k_new);

% Konwersja H(s) -> H(z) z użyciem transformacji biliniowej
[b_d_new, a_d_new] = bilinear(b_a_new, a_a_new, fs);

% Rysowanie charakterystyk z zaznaczonymi częstotliwościami granicznymi
figure;
semilogx(f_analog/(2*pi), 20*log10(abs(h_analog)), 'LineWidth', 1.5, 'Color', [0 0.447 0.741]); hold on;
plot(f_digital, 20*log10(abs(h_digital)), 'LineWidth', 1.5, 'Color', [0.85 0.325 0.098]);

semilogx(f_analog_new/(2*pi), 20*log10(abs(h_analog_new)), 'LineWidth', 1.5, 'Color', [0.466 0.674 0.188]);
plot(f_digital_new, 20*log10(abs(h_digital_new)), 'LineWidth', 1.5, 'Color', [0.494 0.184 0.556]);

xline([f_low_analog, f_high_analog], '--k', 'LineWidth', 1.2);
xline(f_digital_bounds, '--r', 'LineWidth', 1.2);
xlim([1 fs/2]); % Zakres częstotliwości od 1 Hz do Nyquista

legend({'$H(s)$', '$H(z)$', '$H_{w}(s)$', '$H_{w}(z)$', 'Analogowe granice', 'Cyfrowe granice'}, ...
       'Interpreter', 'latex', 'Location', 'Southwest');
title('Charakterystyki z poprawionym pre-warpingiem', 'FontSize', 14);
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
grid on; grid minor;

% Poprawiono błąd w nawiasie
% Obliczenie częstotliwości cyfrowych z pre-warpingiem
f_digital_bounds = 2*fs/pi * atan(pi*[f_low_analog, f_high_analog]/fs);
xline(f_digital_bounds, '--r', 'Cyfrowe granice');
legend('H(s)', 'H(z)', 'H_w(s)', 'H_w(z)', 'Location', 'Southwest');
title('Charakterystyki z poprawionym pre-warpingiem');
xlabel('Częstotliwość (Hz)');
ylabel('Amplituda (dB)');
grid on;

% Generowanie sygnału testowego
t = 0:T:1; % 1 sekunda
x = sin(2*pi*1209*t) + sin(2*pi*1272*t); % składowe 1209 Hz i 1272 Hz

% Implementacja ręczna filtru IIR (znormalizowane współczynniki)
y_manual = zeros(size(x));
a_d_norm = a_d / a_d(1); % Normalizacja (a_d(1) = 1)
b_d_norm = b_d / a_d(1);

% Poprawiona implementacja ręczna filtru
for n = 1:length(x)
    % Inicjalizacja wartości wyjściowej dla bieżącego n
    y_n = 0;
    
    % Sumowanie członów b (numerator) - część wejściowa
    for k = 1:length(b_d_norm)
        if n-k+1 >= 1
            y_n = y_n + b_d_norm(k) * x(n-k+1);
        end
    end
    
    % Sumowanie członów a (denominator) - sprzężenie zwrotne
    for k = 2:length(a_d_norm)  % Zaczynamy od 2, ponieważ a_d_norm(1) = 1
        if n-k+1 >= 1
            y_n = y_n - a_d_norm(k) * y_manual(n-k+1);
        end
    end
    
    % Zapisanie wyniku
    y_manual(n) = y_n;
end

% Filtracja z użyciem filter()
y_filter = filter(b_d_norm, a_d_norm, x);

% Porównanie wyników w dziedzinie czasu
figure;
subplot(2,1,1);
plot(t, y_manual, 'LineWidth', 1.5, 'Color', [0.850 0.325 0.098]);
title('Filtracja ręczna', 'FontSize', 12);
xlabel('Czas [s]'); ylabel('Amplituda');
xlim([1 fs/2]); % Zakres częstotliwości od 1 Hz do Nyquista

grid on; grid minor;

subplot(2,1,2);
plot(t, y_filter, 'LineWidth', 1.5, 'Color', [0 0.447 0.741]);
title('Filtracja z użyciem filter()', 'FontSize', 12);
xlabel('Czas [s]'); ylabel('Amplituda');
xlim([1 fs/2]); % Zakres częstotliwości od 1 Hz do Nyquista

grid on; grid minor;


% Porównanie widma
figure;
plot(f, 20*log10(X_manual), 'r', 'LineWidth', 1.5); hold on;
plot(f, 20*log10(X_filter), 'b--', 'LineWidth', 1.5);
xline([f_low_analog, f_high_analog], '--k', 'LineWidth', 1.2);
legend({'Ręczna', 'filter()'}, 'Location', 'Southwest');
title('Widmo sygnału wyjściowego', 'FontSize', 14);
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda [dB]');
xlim([1 fs/2]); % Zakres częstotliwości od 1 Hz do Nyquista
ylim([-100 max(20*log10(X_manual))+5]);
grid on; grid minor;
