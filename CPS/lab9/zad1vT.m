clearvars; close all;

fs = 8000; 
T = 1;
t = 0:1/fs:T-1/fs;
                
A1 = -0.5; f1 = 34.2;
A2 = 1; f2 = 115.5;
dref1 = A1*sin(2*pi*f1*t) + A2*sin(2*pi*f2*t);     

fc = 1000; 
deltaf = 500; 
fm = 0.25; 
dref2 = sin(2*pi*fc*t + deltaf/fm * sin(2*pi*fm*t));

[dref3, fs3] = audioread('engine.wav'); 
dref3 = dref3(1:fs3);
dref3 = resample(dref3, fs, fs3);
dref3 = dref3.';

drefs = {dref1, dref2, dref3};
titles = {'sinusoidalnych', 'SFM', 'silnika'}; 

for i = 1:3
    dref = drefs{i};

    figure;
    sgtitle(['Porównanie sygnałów ', titles{i}])
    for snr = [10, 20, 40]
        d = awgn( dref, snr, 'measured' );      % WE: sygnał odniesienia dla sygnału x
        x = [ d(1) d(1:end-1) ];                % WE: sygnał filtrowany, teraz opóźniony d
        
        M = 15;             % długość filtru
        mi = 0.01;        % współczynnik szybkości adaptacji
    
        y = []; e = [];     % sygnały wyjściowe z filtra
        bx = zeros(M,1);    % bufor na próbki wejściowe x
        h = zeros(M,1);     % początkowe (puste) wagi filtru
        
        for n = 1 : length(x)
            bx = [ x(n); bx(1:M-1) ];   % pobierz nową próbkę x[n] do bufora
            y(n) = h' * bx;             % oblicz y[n] = sum( x .* bx) – filtr FIR
            e(n) = d(n) - y(n);         % oblicz e[n]
            h = h + mi * e(n) * bx;     % LMS
            % h = h + mi * e(n) * bx /(bx'*bx); % NLMS
        end
        
        y = y(:)';
    
        SNRdb = 10 * log10( sum(dref.^2) / sum((dref - y).^2) ),
    
        subplot(3,1,find([10, 20, 40] == snr))
        plot(t, dref, t, d, t, y);
        legend('referencyjny', 'odniesienia', 'rezultat')
        xlabel('Czas [s]');
        ylabel('Amplituda');
        title(sprintf('SNR = %d dB', snr))
    end
    pause;
end