%% Zadanie 3: Pętla PLL

clear;
close all;

%% Parametry symulacji
fs = 192000;         % Częstotliwość próbkowania [Hz]
fpilot = 19000;      % Częstotliwość pilota [Hz]
duration = 10;       % Czas trwania sygnału [s]
t = 0:1/fs:(duration-1/fs);  % Wektor czasu
N = length(t);       % Liczba próbek

%% 1. Generacja sygnału pilota o stałym przesunięciu fazowym
phase_shift = pi/4;  % Przesunięcie fazowe pilota [rad]
pilot_const = sin(2*pi*fpilot*t + phase_shift);

%% Pętla PLL dla sygnału o stałej częstotliwości
% Inicjalizacja parametrów PLL
freq_const = 2*pi*fpilot/fs;  % Znormalizowana częstotliwość początkowa
theta_const = zeros(1, N+1);  % Wektor fazy
alpha = 1e-2;                 % Parametr alpha pętli PLL
beta = alpha^2/4;             % Parametr beta pętli PLL

% Główna pętla PLL
for n = 1:N
    perr = -pilot_const(n)*sin(theta_const(n));  % Błąd fazowy
    theta_const(n+1) = theta_const(n) + freq_const + alpha*perr;  % Aktualizacja fazy
    freq_const = freq_const + beta*perr;  % Aktualizacja częstotliwości
end

% Generacja sygnałów wyjściowych
c19_const = sin(theta_const(1:end-1));  % Odtworzony pilot 19 kHz
c38_const = sin(2*theta_const(1:end-1));  % Odtworzony sygnał 38 kHz
c57_const = cos(3*theta_const(1:end-1));  % Odtworzony sygnał 57 kHz

% Analiza wyników
figure('Name', 'PLL - Pilot o stałej częstotliwości', 'Position', [100, 100, 1200, 800]);

% Porównanie sygnałów
subplot(3, 1, 1);
plot(t(1:min(5000, N)), pilot_const(1:min(5000, N)), 'b');
hold on;
plot(t(1:min(5000, N)), c19_const(1:min(5000, N)), 'r--');
hold off;
grid on;
title('Porównanie sygnału pilota i odtworzonego sygnału (pierwsze 5000 próbek)');
xlabel('Czas [s]');
ylabel('Amplituda');
legend('Sygnał odniesienia', 'Sygnał odtworzony przez PLL');

% Przebiegi harmonicznych
subplot(3, 1, 2);
plot(t(1:min(1000, N)), c19_const(1:min(1000, N)), 'b');
hold on;
plot(t(1:min(1000, N)), c38_const(1:min(1000, N)), 'r');
plot(t(1:min(1000, N)), c57_const(1:min(1000, N)), 'g');
hold off;
grid on;
title('Przebiegi harmonicznych (pierwsze 1000 próbek)');
xlabel('Czas [s]');
ylabel('Amplituda');
legend('c19 (19 kHz)', 'c38 (38 kHz)', 'c57 (57 kHz)');

% Błąd fazy
subplot(3, 1, 3);
phase_diff = unwrap(theta_const(1:end-1) - (2*pi*fpilot*t + phase_shift));
plot(t, phase_diff);
grid on;
title('Różnica fazowa między sygnałem referencyjnym a odtworzonym');
xlabel('Czas [s]');
ylabel('Różnica fazy [rad]');

%% 2. Generator sygnału pilota o zmiennej częstotliwości
df = 10;  % Odchylenie częstotliwości [Hz]
fm = 0.1; % Częstotliwość modulacji [Hz]

% Generacja sygnału o zmiennej częstotliwości
freq_var = fpilot + df*sin(2*pi*fm*t);  % Częstotliwość chwilowa
phase_var = 2*pi*cumsum(freq_var)/fs;   % Faza chwilowa
pilot_var = sin(phase_var);             % Sygnał pilota

% Pętla PLL dla sygnału o zmiennej częstotliwości
freq_var_pll = 2*pi*fpilot/fs;  % Znormalizowana częstotliwość początkowa
theta_var = zeros(1, N+1);      % Wektor fazy

% Główna pętla PLL
for n = 1:N
    perr = -pilot_var(n)*sin(theta_var(n));  % Błąd fazowy
    theta_var(n+1) = theta_var(n) + freq_var_pll + alpha*perr;  % Aktualizacja fazy
    freq_var_pll = freq_var_pll + beta*perr;  % Aktualizacja częstotliwości
end

% Generacja sygnałów wyjściowych
c19_var = sin(theta_var(1:end-1));  % Odtworzony pilot 19 kHz
c38_var = sin(2*theta_var(1:end-1));  % Odtworzony sygnał 38 kHz
c57_var = cos(3*theta_var(1:end-1));  % Odtworzony sygnał 57 kHz

% Analiza wyników
figure('Name', 'PLL - Pilot o zmiennej częstotliwości', 'Position', [100, 100, 1200, 800]);

% Porównanie sygnałów
subplot(4, 1, 1);
plot(t(1:min(5000, N)), pilot_var(1:min(5000, N)), 'b');
hold on;
plot(t(1:min(5000, N)), c19_var(1:min(5000, N)), 'r--');
hold off;
grid on;
title('Porównanie sygnału pilota i odtworzonego sygnału (pierwsze 5000 próbek)');
xlabel('Czas [s]');
ylabel('Amplituda');
legend('Sygnał odniesienia', 'Sygnał odtworzony przez PLL');

% Zmiana częstotliwości pilota
subplot(4, 1, 2);
plot(t, freq_var);
grid on;
title('Zmiana częstotliwości pilota w czasie');
xlabel('Czas [s]');
ylabel('Częstotliwość [Hz]');

% Estymacja częstotliwości przez PLL
subplot(4, 1, 3);
freq_est = diff(theta_var)/(2*pi)*fs;
plot(t(1:end-1), freq_est);
grid on;
title('Estymacja częstotliwości przez PLL');
xlabel('Czas [s]');
ylabel('Estymowana częstotliwość [Hz]');

% Błąd fazy
subplot(4, 1, 4);
phase_diff_var = unwrap(theta_var(1:end-1) - phase_var);
plot(t, phase_diff_var);
grid on;
title('Różnica fazowa między sygnałem referencyjnym a odtworzonym');
xlabel('Czas [s]');
ylabel('Różnica fazy [rad]');

%% 3. Badanie szybkości zbieżności pętli PLL
% Zgodnie z poleceniem, używamy mocy szumu: 0, 5, 10 i 20 dB
noise_power_values = [0, 5, 10, 20];  % Wartości mocy szumu [dB]
convergence_times = zeros(size(noise_power_values));  % Czasy zbieżności [próbki]

figure('Name', 'PLL - Szybkość zbieżności', 'Position', [100, 100, 1200, 800]);

for i = 1:length(noise_power_values)
    noise_power_dB = noise_power_values(i);
    
    % Generacja zaszumionego sygnału pilota
    if noise_power_dB == 0
        pilot_noisy = pilot_const;  % Bez szumu
    else
        pilot_power = 10*log10(mean(pilot_const.^2));  % Moc sygnału pilota [dB]
        noise_amp = 10^((pilot_power - noise_power_dB)/20);  % Amplituda szumu
        noise = noise_amp * randn(size(pilot_const));  % Szum AWGN
        pilot_noisy = pilot_const + noise;  % Zaszumiony sygnał pilota
    end
    
    % Pętla PLL dla zaszumionego sygnału
    freq_noisy = 2*pi*fpilot/fs;  % Znormalizowana częstotliwość początkowa
    theta_noisy = zeros(1, N+1);  % Wektor fazy
    phase_error = zeros(1, N);  % Wektor błędu fazowego
    
    % Główna pętla PLL
    for n = 1:N
        perr = -pilot_noisy(n)*sin(theta_noisy(n));  % Błąd fazowy
        phase_error(n) = perr;  % Zapisanie błędu fazowego
        theta_noisy(n+1) = theta_noisy(n) + freq_noisy + alpha*perr;  % Aktualizacja fazy
        freq_noisy = freq_noisy + beta*perr;  % Aktualizacja częstotliwości
    end
    
    % Generacja sygnału odtworzonego
    c19_noisy = sin(theta_noisy(1:end-1));  % Odtworzony pilot 19 kHz
    
    % Obliczenie czasu zbieżności (gdy błąd fazowy spada poniżej progu)
    phase_error_smooth = movmean(abs(phase_error), fs/10);  % Wygładzenie błędu fazowego
    threshold = 0.01;  % Próg zbieżności
    conv_indices = find(phase_error_smooth < threshold);
    
    if ~isempty(conv_indices)
        convergence_times(i) = conv_indices(1);
    else
        convergence_times(i) = Inf;
    end
    
    % Wykresy dla każdego poziomu szumu
    subplot(length(noise_power_values), 1, i);
    plot(t(1:min(fs, N)), pilot_noisy(1:min(fs, N)), 'b');
    hold on;
    plot(t(1:min(fs, N)), c19_noisy(1:min(fs, N)), 'r--');
    hold off;
    grid on;
    if noise_power_dB == 0
        title(['Moc szumu = 0 dB (bez szumu), Czas zbieżności = ', num2str(convergence_times(i)), ' próbek']);
    else
        title(['Moc szumu = ', num2str(noise_power_dB), ' dB, Czas zbieżności = ', num2str(convergence_times(i)), ' próbek']);
    end
    xlabel('Czas [s]');
    ylabel('Amplituda');
    legend('Sygnał odniesienia', 'Sygnał odtworzony przez PLL');
end

% Tabela czasów zbieżności
figure('Name', 'PLL - Czasy zbieżności', 'Position', [100, 100, 600, 400]);
bar(convergence_times);
set(gca, 'XTickLabel', {'0', '5', '10', '20'});
grid on;
title('Czasy zbieżności PLL dla różnych wartości mocy szumu');
xlabel('Moc szumu [dB]');
ylabel('Czas zbieżności [próbki]');

% Wnioski
disp('Wyniki analizy szybkości zbieżności PLL:');
for i = 1:length(noise_power_values)
    if noise_power_values(i) == 0
        disp(['Moc szumu = 0 dB (bez szumu): ', num2str(convergence_times(i)), ' próbek (', num2str(convergence_times(i)/fs*1000), ' ms)']);
    else
        disp(['Moc szumu = ', num2str(noise_power_values(i)), ' dB: ', num2str(convergence_times(i)), ' próbek (', num2str(convergence_times(i)/fs*1000), ' ms)']);
    end
end