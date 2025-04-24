% ZADANIE 2 - Dekodowanie DTMF

% Wczytanie sygnału wzorcowego
[s, fs] = audioread('s.wav');

% Rysowanie spektrogramu
figure;
spectrogram(s, 4096, 4096-512, 0:5:2000, fs, 'yaxis');
title('Spektrogram sygnalu DTMF');

% Filtrowanie sygnału przy uzyciu filtru z zadania 1 (b_d, a_d)
load('butter.mat');
[b_a, a_a] = zp2tf(z, p, k);
[b_d, a_d] = bilinear(b_a, a_a, fs);

s_filtered = filter(b_d, a_d, s);

% Porownanie spektrogramow przed i po filtracji
figure;
subplot(2,1,1);
spectrogram(s, 4096, 4096-512, 0:5:2000, fs, 'yaxis');
title('Przed filtracja');

subplot(2,1,2);
spectrogram(s_filtered, 4096, 4096-512, 0:5:2000, fs, 'yaxis');
title('Po filtracji');

% OPCJA: Analiza Goertzla
frequencies = [697 770 852 941 1209 1336 1477];
N = 1024;
window = s(1:N);  % fragment sygnalu
result = zeros(size(frequencies));

for i = 1:length(frequencies)
    k = round(frequencies(i)*N/fs);
    omega = 2*pi*k/N;
    coeff = 2*cos(omega);
    q1 = 0; q2 = 0;

    for n = 1:N
        q0 = window(n) + coeff*q1 - q2;
        q2 = q1;
        q1 = q0;
    end

    result(i) = q1^2 + q2^2 - coeff*q1*q2;
end

figure;
bar(frequencies, result);
xlabel('Czestotliwosc [Hz]');
ylabel('Energia');
title('Wynik algorytmu Goertzla');
