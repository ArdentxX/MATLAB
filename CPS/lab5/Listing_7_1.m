% cps_07_analog_intro.m
clear all; close all;

% Zaprojektuj/dobierz wspolczynniki transmitancji filtra analogowego
if(0) % dobor wartosci wspolczynnikow wielomianow zmiennej 's' transmitancji
    b = [3,2]; % [ b1, b0 ]
    a = [4,3,2,1]; % [ a3, a2, a1, a0=1]
    z = roots(b); p = roots(a); % [b,a] --> [z,p]
else % dobor wartosci pierwiastkow wielomianow zmiennej 's' transmitancji
    wzm = 0.001;
    z = j*2*pi*[ 600,800 ]; z = [z conj(z)];
    p = [-1,-1] + j*2*pi*[100,200]; p = [p conj(p)];
    b = wzm*poly(z); a = poly(p); % [z,p] --> [b,a]
end
figure; plot(real(z),imag(z),'bo', real(p),imag(p),'r*'); grid;

title('Zera (o) i Bieguny (*)'); xlabel('Real()'); ylabel('Imag()'); pause

% Weryfikacja odpowiedzi (charakterystyki) filtra: amplitudowej fazowej, impulsowej, skokowej
f = 0 : 0.1 : 1000; % czestotliwosc w hercach
w = 2*pi*f; % czestotliwosc katowa, pulsacja
s = j*w; % zmienna transformacji Laplace'a
H = polyval(b,s) ./ polyval(a,s); % h(s) --> H(f) dla s=j*w: iloraz dwoch wielomianow
figure; plot(f,20*log10(abs(H))); xlabel('f [Hz]'); title('|H(f)| [dB]'); grid; pause
figure; plot(f,unwrap(angle(H))); xlabel('f [Hz]'); title('angle(H(f)) [rad]'); grid; pause
figure; impulse(b,a); pause % odpowiedz filtra na pobudzenie impulsowe (z Control Toolbox)
figure; step(b,a); pause % odpowiedz filtra na pobudzenie skokowe (z Control Toolbox)
