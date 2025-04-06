% cps_07_analog_tranform.m
clear all; close all;
% Wymagania
N = 8; % liczba biegunow transmitanci prototypu analogowego
f0 = 100; % dla filtrow LowPass oraz HighPass
f1 = 10; f2=100; % dla filtrow BandPass oraz BandStop
Rp = 3; % dozwolony poziom oscylacji w pasmie przepustowym (w dB) - ripples in pass
Rs = 100; % dozwolony poziom oscylacji w pasmie zaporowym (w dB) - ripples in stop

% Projekt analogowego filtra prototypowego - dolnoprzepustowy z w0 = 1
[z,p,wzm] = buttap(N); % analogowy prototyp Butterwotha
% [z,p,wzm] = cheb1ap(N,Rp); % analogowy prototyp Czebyszewa typu-I
% [z,p,wzm] = cheb2ap(N,Rs); % analogowy prototyp Czebyszewa typu-II
% [z,p,wzm] = ellipap(N,Rp,Rs); % analogowy prototyp Cauera (eliptyczny)
b = wzm*poly(z); % [z,wzm] --> b
a = poly(p); % p --> a
f = 0 : 0.01 : 1000; w=2*pi*f; % zakres czestotliwosci
H = freqs(b,a,w); % odpowiedz czestotliwosciowa filtra (funkcja Matlaba)
fi=0:pi/1000:2*pi; c=cos(fi); s=sin(fi);
figure; semilogx(w,20*log10(abs(H))); grid; xlabel('w [rad*Hz]'); title('Analog Proto |H(f)'); pause
figure; plot(real(z),imag(z),'ro',real(p),imag(p),'b*',c,s,'k-'); grid; title('Analog Proto ZP'); pause

% Transformacja czestotliwosci: znormalizowany (w0=1) LowPass --> xxxPass, xxxStop
% Funkcje xx2yy z Signal Processing Toolbox
[b,a] = lp2lp(b,a,2*pi*f0); % LowPass na LowPass
% [b,a] = lp2hp(b,a,2*pi*f0); % LowPass na HighPass
% [b,a] = lp2bp(b,a,2*pi*sqrt(f1*f2),2*pi*(f2-f1)); % LowPass na BandPass
% [b,a] = lp2bs(b,a,2*pi*sqrt(f1*f2),2*pi*(f2-f1)); % LowPass na BandStop
z=roots(b); p=roots(a);
% ... kontynuuj cps_07_analog_intro.m