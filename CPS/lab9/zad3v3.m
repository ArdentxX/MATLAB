%% Zadanie 3: Pętla PLL
% Implementacja cyfrowej pętli PLL do odtworzenia sygnałów o stałej i zmiennej częstotliwości
% oraz badanie szybkości zbieżności w zależności od poziomu szumu

clear;
close all;

%% Parametry symulacji
fs = 192000;         % Częstotliwość próbkowania [Hz]
fpilot = 19000;      % Częstotliwość pilota [Hz]
duration = 10;       % Czas trwania sygnału [s]
t = 0:1/fs:(duration-1/fs);  % Wektor czasu
N = length(t);       % Liczba próbek

%% 1. Generacja sygnału pilota o stałym przesunięciu fazowym
% Wygeneruj sygnał udający pilota 19 kHz o stałym przesunięciu fazowym 
% i sprawdź, czy adaptacyjny oscylator się do niego dostroi.
phase_shift = pi/4;  % Przesunięcie fazowe pilota [rad]
pilot_const = sin(2*pi*fpilot*t + phase_shift);

% Pętla PLL dla sygnału o stałej częstotliwości
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
c38_const = sin(2*theta_const(1:end-1));  % Odtworzony sygnał 38 kHz (dla stereo)
c57_const = cos(3*theta_const(1:end-1));  % Odtworzony sygnał 57 kHz (dla RDS)

% Analiza wyników
figure('Name', 'PLL - Pilot o stałej częstotliwości');
% Porównanie sygnałów
subplot(2, 1, 1);
plot(t(1:min(5000, N)), pilot_const(1:min(5000, N)), 'b');
hold on;
plot(t(1:min(5000, N)), c19_const(1:min(5000, N)), 'r--');
hold off;
grid on;
title('Porównanie sygnału pilota i odtworzonego sygnału (pierwsze 5000 próbek)');
xlabel('Czas [s]');
ylabel('Amplituda');
legend('Sygnał odniesienia', 'Sygnał odtworzony przez PLL');

% Błąd fazy
subplot(2, 1, 2);
phase_diff = unwrap(theta_const(1:end-1) - (2*pi*fpilot*t + phase_shift));
plot(t, phase_diff);
grid on;
title('Różnica fazowa między sygnałem referencyjnym a odtworzonym');
xlabel('Czas [s]');
ylabel('Różnica fazy [rad]');

%% 2. Generator sygnału pilota o zmiennej częstotliwości
% Wygeneruj sygnał j.w. tylko niech częstotliwość pilota dodatkowo
% wolno się zmienia sinusoidalnie: ±10 Hz (df=10 Hz) jeden raz na 10 sekund (fm=0.1 Hz)
df = 10;  % Odchylenie częstotliwości [Hz]
fm = 0.1; % Częstotliwość modulacji [Hz]

% Generacja sygnału o zmiennej częstotliwości
freq_var = fpilot + df*sin(2*pi*fm*t);  % Częstotliwość chwilowa
phase_var = 2*pi*cumsum(freq_var)/fs;   % Faza chwilowa
pilot_var = sin(phase_var);             % Sygnał pilota

% Pętla PLL dla sygnału o zmiennej częstotliwości
% Sprawdzenie czy adaptacyjny oscylator dostroi się do pilota o zmiennej częstotliwości
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

% Analiza wyników
figure('Name', 'PLL - Pilot o zmiennej częstotliwości');
% Porównanie sygnałów
subplot(3, 1, 1);
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
subplot(3, 1, 2);
plot(t, freq_var);
grid on;
title('Zmiana częstotliwości pilota w czasie');
xlabel('Czas [s]');
ylabel('Częstotliwość [Hz]');

% Błąd fazy
subplot(3, 1, 3);
phase_diff_var = unwrap(theta_var(1:end-1) - phase_var);
plot(t, phase_diff_var);
grid on;
title('Różnica fazowa między sygnałem referencyjnym a odtworzonym');
xlabel('Czas [s]');
ylabel('Różnica fazy [rad]');

%% 3. Badanie szybkości zbieżności pętli PLL
% Sprawdź szybkość zbieżności pętli PLL. W tym celu do sygnału 
% z pkt. 1 dodaj szum AWGN o mocy: 0, 5, 10 i 20 dB. Znając sygnał oczekiwany (wzorcowy), 
% określ po ilu próbkach oscylator dostroił się do sygnału.
SNR_values = [inf, 20, 10, 5, 0];  % Wartości SNR [dB]
convergence_times = zeros(size(SNR_values));  % Czasy zbieżności [próbki]

figure('Name', 'PLL - Szybkość zbieżności');

for i = 1:length(SNR_values)
    SNR = SNR_values(i);
    
    % Generacja zaszumionego sygnału pilota
    if isinf(SNR)
        pilot_noisy = pilot_const;  % Bez szumu
    else
        pilot_power = mean(pilot_const.^2);  % Moc sygnału pilota
        noise_power = pilot_power / (10^(SNR/10));  % Wymagana moc szumu
        noise = sqrt(noise_power) * randn(size(pilot_const));  % Szum AWGN
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
    
    % Wykresy dla każdego SNR
    subplot(length(SNR_values), 1, i);
    plot(t(1:min(fs, N)), pilot_noisy(1:min(fs, N)), 'b');
    hold on;
    plot(t(1:min(fs, N)), c19_noisy(1:min(fs, N)), 'r--');
    hold off;
    grid on;
    if isinf(SNR)
        title(['SNR = \infty dB, Czas zbieżności = ', num2str(convergence_times(i)), ' próbek']);
    else
        title(['SNR = ', num2str(SNR), ' dB, Czas zbieżności = ', num2str(convergence_times(i)), ' próbek']);
    end
    xlabel('Czas [s]');
    ylabel('Amplituda');
    legend('Sygnał zaszumiony', 'Sygnał odtworzony przez PLL');
end

% Tabela czasów zbieżności
figure('Name', 'PLL - Czasy zbieżności');
bar(convergence_times);
set(gca, 'XTickLabel', {'∞', '20', '10', '5'});
grid on;
title('Czasy zbieżności PLL dla różnych wartości SNR');
xlabel('SNR [dB]');
ylabel('Czas zbieżności [próbki]');