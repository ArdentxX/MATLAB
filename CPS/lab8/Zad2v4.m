%% === ZADANIE 2: Modulacja i demodulacja AM (DSB-C, DSB-SC, SSB-SC) ===
clear; clc; close all;

%% === [1] Wczytanie i przygotowanie sygnałów audio ===
[x1, fsx] = audioread('mowa8000.wav');
x2 = flipud(x1);

fs = 400000;
fc1 = 100000;
fc2 = 110000;
dA = 0.25;

x1_up = resample(x1, fs, fsx);
x2_up = resample(x2, fs, fsx);
t = (0:length(x1_up)-1)' / fs;

%% === [2] Filtr Hilberta FIR ===
N = 101;
n = -(N-1)/2 : (N-1)/2;
h = zeros(size(n));
h(n ~= 0) = (2 ./ (pi * n(n ~= 0))) .* sin(pi * n(n ~= 0) / 2);
h = h .* hamming(N)';

x1h = conv(x1_up, h, 'same');
x2h = conv(x2_up, h, 'same');

%% === [3] Modulacje ===
x1_dsb_c = (1 + dA * x1_up) .* cos(2*pi*fc1*t);
x2_dsb_c = (1 + dA * x2_up) .* cos(2*pi*fc2*t);
yDSB_C = x1_dsb_c + x2_dsb_c;

x1_dsb_sc = dA * x1_up .* cos(2*pi*fc1*t);
x2_dsb_sc = dA * x2_up .* cos(2*pi*fc2*t);
yDSB_SC = x1_dsb_sc + x2_dsb_sc;

x1_ssb = 0.5 * x1_up .* cos(2*pi*fc1*t) + 0.5 * x1h .* sin(2*pi*fc1*t);
x2_ssb = 0.5 * x2_up .* cos(2*pi*fc2*t) - 0.5 * x2h .* sin(2*pi*fc2*t);
ySSB_SC = x1_ssb + x2_ssb;

%% === [4] Filtr dolnoprzepustowy do demodulacji ===
lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', 4000, ...
    'StopbandFrequency', 6000, ...
    'SampleRate', fs, ...
    'DesignMethod', 'equiripple');

%% === [5] DEMODULACJA ===
% DSB-C
s1 = filter(lpFilt, yDSB_C .* cos(2*pi*fc1*t));
s1h = conv(s1, h, 'same');
dec1_dsb_c = (sqrt(s1.^2 + s1h.^2) - 1) / dA;
dec1_dsb_c = dec1_dsb_c / max(abs(dec1_dsb_c));
dec1_dsb_c_down = resample(dec1_dsb_c, fsx, fs);

s2 = filter(lpFilt, yDSB_C .* cos(2*pi*fc2*t));
s2h = conv(s2, h, 'same');
dec2_dsb_c = (sqrt(s2.^2 + s2h.^2) - 1) / dA;
dec2_dsb_c = dec2_dsb_c / max(abs(dec2_dsb_c));
dec2_dsb_c_down = resample(dec2_dsb_c, fsx, fs);

% DSB-SC
s1 = filter(lpFilt, yDSB_SC .* cos(2*pi*fc1*t));
s1h = conv(s1, h, 'same');
dec1_dsb_sc = (sqrt(s1.^2 + s1h.^2) - 1) / dA;
dec1_dsb_sc = dec1_dsb_sc / max(abs(dec1_dsb_sc));
dec1_dsb_sc_down = resample(dec1_dsb_sc, fsx, fs);

s2 = filter(lpFilt, yDSB_SC .* cos(2*pi*fc2*t));
s2h = conv(s2, h, 'same');
dec2_dsb_sc = (sqrt(s2.^2 + s2h.^2) - 1) / dA;
dec2_dsb_sc = dec2_dsb_sc / max(abs(dec2_dsb_sc));
dec2_dsb_sc_down = resample(dec2_dsb_sc, fsx, fs);

% SSB-SC
s1 = filter(lpFilt, ySSB_SC .* cos(2*pi*fc1*t));
s1h = conv(s1, h, 'same');
dec1_ssb = s1 - s1h;
dec1_ssb = dec1_ssb / max(abs(dec1_ssb));
dec1_ssb_down = resample(dec1_ssb, fsx, fs);

s2 = filter(lpFilt, ySSB_SC .* cos(2*pi*fc2*t));
s2h = conv(s2, h, 'same');
dec2_ssb = s2 + s2h;
dec2_ssb = dec2_ssb / max(abs(dec2_ssb));
dec2_ssb_down = resample(dec2_ssb, fsx, fs);

%% === [6] ODSŁUCH ===
disp('Odtwarzanie: DSB-C Stacja 1');
sound(dec1_dsb_c_down, fsx); pause(length(x1)/fsx + 1);

disp('Odtwarzanie: DSB-C Stacja 2');
sound(flipud(dec2_dsb_c_down), fsx); pause(length(x1)/fsx + 1);

disp('Odtwarzanie: DSB-SC Stacja 1');
sound(dec1_dsb_sc_down, fsx); pause(length(x1)/fsx + 1);

disp('Odtwarzanie: DSB-SC Stacja 2');
sound(flipud(dec2_dsb_sc_down), fsx); pause(length(x1)/fsx + 1);

disp('Odtwarzanie: SSB-SC Stacja 1');
sound(dec1_ssb_down, fsx); pause(length(x1)/fsx + 1);

disp('Odtwarzanie: SSB-SC Stacja 2');
sound(flipud(dec2_ssb_down), fsx); pause(length(x1)/fsx + 1);

%% === [7] Opcja dodatkowa: dwie stacje na jednej nośnej ===
ySSB_sameCarrier = ...
    0.5 * x1_up .* cos(2*pi*fc1*t) + 0.5 * x1h .* sin(2*pi*fc1*t) + ...
    0.5 * x2_up .* cos(2*pi*fc1*t) - 0.5 * x2h .* sin(2*pi*fc1*t);

s1 = filter(lpFilt, ySSB_sameCarrier .* cos(2*pi*fc1*t));
s1h = conv(s1, h, 'same');
dec_usb = s1 - s1h;
dec_usb = dec_usb / max(abs(dec_usb));
dec_usb_down = resample(dec_usb, fsx, fs);

s2 = filter(lpFilt, ySSB_sameCarrier .* cos(2*pi*fc1*t));
s2h = conv(s2, h, 'same');
dec_lsb = s2 + s2h;
dec_lsb = dec_lsb / max(abs(dec_lsb));
dec_lsb_down = resample(dec_lsb, fsx, fs);

disp('Odtwarzanie: Opcja dodatkowa - Stacja 1 (USB)');
sound(dec_usb_down, fsx); pause(length(x1)/fsx + 1);

disp('Odtwarzanie: Opcja dodatkowa - Stacja 2 (LSB)');
sound(flipud(dec_lsb_down), fsx); pause(length(x1)/fsx + 1);

%% === [8] PORÓWNANIE OBWIEDNI — WYKRESY ===
figure('Name', 'Obwiednie DSB-C'); 
subplot(2,1,1);
plot(dec1_dsb_c(1:4000)); title('DSB-C Stacja 1'); xlabel('Próbki'); ylabel('Amplituda');
subplot(2,1,2);
plot(dec2_dsb_c(1:4000)); title('DSB-C Stacja 2'); xlabel('Próbki'); ylabel('Amplituda');

figure('Name', 'Obwiednie DSB-SC'); 
subplot(2,1,1);
plot(dec1_dsb_sc(1:4000)); title('DSB-SC Stacja 1'); xlabel('Próbki'); ylabel('Amplituda');
subplot(2,1,2);
plot(dec2_dsb_sc(1:4000)); title('DSB-SC Stacja 2'); xlabel('Próbki'); ylabel('Amplituda');

figure('Name', 'Obwiednie SSB-SC'); 
subplot(2,1,1);
plot(dec1_ssb(1:4000)); title('SSB-SC Stacja 1 (USB)'); xlabel('Próbki'); ylabel('Amplituda');
subplot(2,1,2);
plot(dec2_ssb(1:4000)); title('SSB-SC Stacja 2 (LSB)'); xlabel('Próbki'); ylabel('Amplituda');

figure('Name', 'SSB-SC dwie stacje na jednej nośnej');
subplot(2,1,1);
plot(dec_usb(1:4000)); title('Stacja 1 (USB) - jedna nośna'); xlabel('Próbki'); ylabel('Amplituda');
subplot(2,1,2);
plot(dec_lsb(1:4000)); title('Stacja 2 (LSB) - jedna nośna'); xlabel('Próbki'); ylabel('Amplituda');
