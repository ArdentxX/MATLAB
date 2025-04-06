clearvars; close all;

f0 = 100; % Częstotliwość graniczna w Hz
omega_3dB = 2 * pi * f0; % Częstotliwość graniczna w rad/s
frequencies = 0:10:10000; % Zakres częstotliwości w Hz

% Pętla dla każdego rzędu filtra
for N = 2:2:8
    % Oblicz bieguny filtra Butterwortha
    poles = zeros(1, N);
    for k = 1:N
        theta_k = pi / 2 + (k - 1) * pi / N + (1 / (2 * N)) * pi;
        poles(k) = omega_3dB * exp(1j * theta_k);
    end
    
    % Oblicz współczynniki filtra
    [b, a] = zp2tf([], poles, omega_3dB^N);
    
    % Wyznacz funkcję transmitancji
    H = tf(b, a);

    % Obliczenie charakterystyki amplitudowej i fazowej
    [H_w, w] = freqs(b, a, 2 * pi * frequencies);
    amplitude = 20 * log10(abs(H_w)); % Charakterystyka amplitudowa
    phase = angle(H_w) * (180 / pi); % Charakterystyka fazowa w stopniach

    % Rysowanie charakterystyk amplitudowych
    figure(1);
    subplot(2, 1, 1);
    plot(frequencies, amplitude, 'DisplayName', ['N = ' num2str(N)]); hold on;
    xlabel('Częstotliwość [Hz]'); ylabel('Amplituda [dB]');
    title('Charakterystyka amplitudowa (oś liniowa)');
    legend show; grid on;
    
    subplot(2, 1, 2);
    semilogx(frequencies, amplitude, 'DisplayName', ['N = ' num2str(N)]); hold on;
    xlabel('Częstotliwość [Hz]'); ylabel('Amplituda [dB]');
    title('Charakterystyka amplitudowa (oś logarytmiczna)');
    legend show; grid on;

    % Rysowanie charakterystyk fazowych
    figure(2);
    plot(frequencies, phase, 'DisplayName', ['N = ' num2str(N)]); hold on;
    xlabel('Częstotliwość [Hz]'); ylabel('Faza [stopnie]');
    title('Charakterystyka fazowa');
    legend show; grid on;
    
    % Dla filtru N=4: Wyznacz odpowiedź impulsową i skokową
    if N == 4
        figure(3);
        % Odpowiedź impulsowa
        subplot(2, 1, 1);
        impulse(H);
        title('Odpowiedź impulsowa filtru N=4');
        xlabel('Czas [s]'); ylabel('Amplituda');
        grid on;

        % Odpowiedź na skok jednostkowy
        subplot(2, 1, 2);
        step(H);
        title('Odpowiedź na skok jednostkowy filtru N=4');
        xlabel('Czas [s]'); ylabel('Amplituda');
        grid on;
    end
    pause;
end