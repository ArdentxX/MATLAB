clearvars; close all;

% Parametry filtru
fs = 256e3; % Częstotliwość próbkowania [Hz]
f3dB = 64e3; % Pasmo przenoszenia [Hz]
fn = fs / 2; % Połowa częstotliwości próbkowania (Nyquist) [Hz]
An = 40; % Minimalne tłumienie w paśmie zaporowym [dB]

% Zamiana częstotliwości na pulsacje [rad/s]
omega3dB = 2 * pi * f3dB;
omega_nyquist = 2 * pi * fn;

% Projektowanie filtrów
filter_types = {'Butterworth', 'Czebyszew 1', 'Czebyszew 2', 'Eliptyczny'};
figure;

for idx = 1:length(filter_types)
    filter_type = filter_types{idx};
    switch filter_type
        case 'Butterworth'
            [N, Wn] = buttord(omega3dB, omega_nyquist, 3, An, 's'); % Minimalny rząd
            [b, a] = butter(N, Wn, 'low', 's'); % Filtr LP Butterworth
            
        case 'Czebyszew 1'
            [N, Wn] = cheb1ord(omega3dB, omega_nyquist, 3, An, 's');
            [b, a] = cheby1(N, 3, Wn, 'low', 's'); % Filtr LP Czebyszew 1
            
        case 'Czebyszew 2'
            [N, Wn] = cheb2ord(omega3dB, omega_nyquist, 3, An, 's');
            [b, a] = cheby2(N, An, Wn, 'low', 's'); % Filtr LP Czebyszew 2
            
        case 'Eliptyczny'
            [N, Wn] = ellipord(omega3dB, omega_nyquist, 3, An, 's');
            [b, a] = ellip(N, 3, An, Wn, 'low', 's'); % Filtr LP eliptyczny
    end
    
    % Rysowanie biegunów
    subplot(2, length(filter_types), idx);
    zplane([], roots(a)); % Bieguny
    title([filter_type ' - Bieguny']);
    xlabel('Re'); ylabel('Im'); grid on;
    
    % Charakterystyka częstotliwościowa
    f = linspace(0, fs / 2, 1000); % Zakres częstotliwości [Hz]
    omega = 2 * pi * f; % Zamiana na pulsacje [rad/s]
    H = freqs(b, a, omega); % Obliczenie transmitancji
    
    subplot(2, length(filter_types), idx + length(filter_types));
    semilogx(f, 20 * log10(abs(H)), 'LineWidth', 2); % Wykres w skali log
    hold on;
    plot([f3dB f3dB], [-100 3], '--r'); % Linia f3dB
    plot([fn fn], [-100 3], '--g'); % Linia fs/2
    title([filter_type ' - Charakterystyka |H(j\omega)|']);
    xlabel('Częstotliwość [Hz]');
    ylabel('20log10(|H(j\omega)|) [dB]');
    grid on;
end

%{
Filtr eliptyczny jest najkorzystniejszy, 
ponieważ minimalizuje zarówno złożoność układu, jak i spełnia założone wymagania projektowe. 
Jednakże wybór filtru może również zależeć od specyficznych wymagań projektowych, 
takich jak tolerancja na zafalowania w paśmie przepustowym czy możliwości realizacji w układzie analogowym.
%}