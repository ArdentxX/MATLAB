clearvars; close all; clc;

% Wymagania dla projektu filtra
N = 8; % liczba biegunów filtra prototypowego
f0 = 100; % dla filtrów LowPass i HighPass [Hz]
f1 = 10; f2 = 100; % dla filtrów BandPass i BandStop [Hz]
Rp = 3; % oscylacje w paśmie przepustowym [dB]
Rs = 100; % tłumienie w paśmie zaporowym [dB]

% Typy prototypów filtrów
prototype_types = {'Butterworth', 'Czebyszew-I', 'Czebyszew-II', 'Eliptyczny'};
filter_transformations = {'LowPass', 'HighPass', 'BandPass', 'BandStop'};

for proto_idx = 1:length(prototype_types)
    prototype_type = prototype_types{proto_idx};

    % **Projekt analogowego filtra prototypowego**
    switch prototype_type
        case 'Butterworth'
            [z, p, wzm] = buttap(N); % Butterworth
        case 'Czebyszew-I'
            [z, p, wzm] = cheb1ap(N, Rp); % Czebyszew typu I
        case 'Czebyszew-II'
            [z, p, wzm] = cheb2ap(N, Rs); % Czebyszew typu II
        case 'Eliptyczny'
            [z, p, wzm] = ellipap(N, Rp, Rs); % Eliptyczny
    end

    b = wzm * poly(z); % Licznik transmitancji
    a = poly(p); % Mianownik transmitancji

    % **Charakterystyka przed transformacją (prototyp)**
    f = 0:0.1:1000; % Częstotliwości [Hz]
    w = 2 * pi * f; % Pulsacje [rad/s]
    H = freqs(b, a, w); % Charakterystyka częstotliwościowa filtra

    figure;
    semilogx(f, 20*log10(abs(H))); grid;
    xlabel('Częstotliwość [Hz]');
    ylabel('|H(f)| [dB]');
    title(['Charakterystyka Prototypu: ' prototype_type]);
    pause;
    
    figure;
    % **Transformacje filtra prototypowego na różne typy**
    for trans_idx = 1:length(filter_transformations)
        filter_type = filter_transformations{trans_idx};
        switch filter_type
            case 'LowPass'
                [b_t, a_t] = lp2lp(b, a, 2*pi*f0); % Transformacja LowPass na LowPass
            case 'HighPass'
                [b_t, a_t] = lp2hp(b, a, 2*pi*f0); % Transformacja LowPass na HighPass
            case 'BandPass'
                [b_t, a_t] = lp2bp(b, a, 2*pi*sqrt(f1*f2), 2*pi*(f2-f1)); % LowPass na BandPass
            case 'BandStop'
                [b_t, a_t] = lp2bs(b, a, 2*pi*sqrt(f1*f2), 2*pi*(f2-f1)); % LowPass na BandStop
        end

        % Charakterystyka po transformacji
        H_t = freqs(b_t, a_t, w);
        subplot(2, 2, trans_idx);
        semilogx(f, 20*log10(abs(H_t)), 'LineWidth', 1.5); grid;
        xlabel('Częstotliwość [Hz]');
        ylabel('|H(f)| [dB]');
        title(['Charakterystyka: ' filter_type ' (Prototyp: ' prototype_type ')']);
    end
    pause;
end