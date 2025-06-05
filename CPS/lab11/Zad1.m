%% CYFROWE PRZETWARZANIE SYGNAŁÓW - KOMPRESJA DŹWIĘKU
clear all; close all; clc;

%% SEKCJA 1: WCZYTANIE I PRZYGOTOWANIE SYGNAŁU
fprintf('=== ŁADOWANIE SYGNAŁU AUDIO ===\n');
% Próba wczytania pliku audio (można zmienić ścieżkę)
% Pamiętaj, że plik 'DontWorryBeHappy.wav' musi być w tym samym katalogu
% co skrypt MATLAB lub w ścieżce wyszukiwania MATLABa.
[x, Fs] = audioread('DontWorryBeHappy.wav');
fprintf('Wczytano plik: DontWorryBeHappy.wav\n');

% Konwersja na mono i normalizacja
if size(x, 2) > 1
    x = mean(x, 2); % Konwersja stereo na mono
end
x = x(:); % Wektor kolumnowy
x = x / max(abs(x)); % Normalizacja do zakresu [-1, 1] (ważne dla stabilności)

% Parametry sygnału
N = length(x);
fprintf('Częstotliwość próbkowania: %d Hz\n', Fs);
fprintf('Długość sygnału: %.2f s (%d próbek)\n', N/Fs, N);

%% KLUCZOWE WZORY Z ZADANIA:
% KODER DPCM:   d(n) = x(n) - a*x(n-1)       gdzie dq = Quantization(d)
% DEKODER DPCM: y(n) = dq(n) + a*y(n-1)      <- REKURENCJA (nie x(n-1)!)
% Parametr predykcji: a = 0.9545
% Kwantyzacja: 4 bity, 16 poziomów

% Parametr predykcji
a = 0.9545;
fprintf('Współczynnik predykcji a = %.4f\n', a);

%% KODER DPCM
fprintf('Kodowanie DPCM...\n');
d = zeros(size(x)); % Sygnał różnicowy
% x_pred = zeros(size(x)); % Niepotrzebne, jeśli nie wizualizujemy

% Pętla kodowania DPCM (dla sygnału x)
for n = 2:N
    % x_pred(n) = a * x(n-1); % Predykcja - niepotrzebne do samego obliczenia d
    d(n) = x(n) - a * x(n-1); % Sygnał różnicowy
end

% --- ZMODYFIKOWANA FUNKCJA KWANTYZACJI ---
% Nowa funkcja kwantyzacji będzie dynamicznie dopasowywać zakres
% do min/max wartości sygnału wejściowego, co zapobiegnie clippingowi.
% Parametr 'num_bits' określa liczbę bitów (np. 4 dla 16 poziomów).
function dq_output = dynamic_quantize(d_input, num_bits)
    levels = 2^num_bits; % Liczba poziomów kwantyzacji
    
    % Wyznaczenie minimalnej i maksymalnej wartości danych w wektorze d
    % Możemy zastosować lekki bufor, aby uniknąć problemów z wartościami skrajnymi
    min_val = min(d_input);
    max_val = max(d_input);
    
    % Jeśli min_val i max_val są takie same (np. sygnał jest stały),
    % ustaw domyślny zakres, aby uniknąć dzielenia przez zero.
    if abs(max_val - min_val) < eps % eps to bardzo mała liczba
        min_val = -1; % Domyślny zakres
        max_val = 1;
    end
    
    data_range = max_val - min_val;
    step = data_range / (levels - 1); % Krok kwantyzacji (levels-1 interwałów)
    
    % Kwantyzacja przez mapowanie na poziomy i zaokrąglenie
    % Przesuwamy zakres do startu od 0, skalujemy, zaokrąglamy, przesuwamy z powrotem
    dq_output = round((d_input - min_val) / step) * step + min_val;
    
    % Opcjonalne delikatne ograniczenie do oryginalnego zakresu wejściowego
    % Może być potrzebne dla wartości minimalnej/maksymalnej, aby idealnie trafić w poziomy
    dq_output = max(min_val, min(max_val, dq_output));
end

% Kwantyzacja sygnału różnicowego
num_quant_bits = 4; % Zgodnie z zadaniem (16 poziomów)
fprintf('Używam funkcji dynamic_quantize (%d bitów, %d poziomów)...\n', num_quant_bits, 2^num_quant_bits);
dq = dynamic_quantize(d, num_quant_bits);
fprintf('Kwantyzacja sygnału różnicowego "d" zakończona.\n');

%% DEKODER DPCM (Dwa warianty: bez i z kwantyzacją)

% Wariant 1: Dekoder DPCM BEZ kwantyzacji (rekonstrukcja z 'd')
fprintf('Dekodowanie DPCM (sygnał różnicowy BEZ kwantyzacji)...\n');
y_unquantized = zeros(size(x)); % Zrekonstruowany sygnał (bez kwantyzacji)
for n = 2:N
    y_unquantized(n) = d(n) + a * y_unquantized(n-1); % Rekonstrukcja z 'd'
end
fprintf('Rekonstrukcja sygnału bez kwantyzacji zakończona.\n');

% Wariant 2: Dekoder DPCM Z kwantyzacją (rekonstrukcja z 'dq')
fprintf('Dekodowanie DPCM (sygnał różnicowy Z kwantyzacją)...\n');
y = zeros(size(x)); % Zrekonstruowany sygnał (z kwantyzacją)
for n = 2:N
    y(n) = dq(n) + a * y(n-1); % Rekonstrukcja z 'dq'
end
fprintf('Rekonstrukcja sygnału z kwantyzacją zakończona.\n');

%% SEKCJA 3: IMPLEMENTACJA ADPCM (Adaptive DPCM)
fprintf('\n=== IMPLEMENTACJA ADPCM (Adaptive DPCM) ===\n');
% Parametry adaptacji
mu = 0.01; % Współczynnik adaptacji (learning rate)
a_adpcm = 0.5; % Początkowa wartość współczynnika predykcji (można eksperymentować)

% Inicjalizacja zmiennych dla ADPCM
d_adpcm = zeros(size(x));
y_adpcm = zeros(size(x));
a_history = zeros(size(x)); % Historia współczynników predykcji
a_history(1) = a_adpcm; % Zapisujemy początkową wartość

fprintf('Kodowanie/Dekodowanie ADPCM (z kwantyzacją %d bitów)...\n', num_quant_bits);
fprintf('Współczynnik adaptacji mu = %.4f\n', mu);

% Pętla ADPCM
for n = 2:N
    % Predykcja z aktualnym współczynnikiem (na podstawie poprzednio zrekonstruowanego sygnału)
    x_pred_adpcm = a_adpcm * y_adpcm(n-1);
    
    % Sygnał różnicowy
    d_adpcm(n) = x(n) - x_pred_adpcm;
    
    % Kwantyzacja sygnału różnicowego ADPCM (używamy tej samej funkcji dynamic_quantize)
    % Ważne: kwantyzacja ADPCM zazwyczaj dynamicznie dostosowuje swój zakres na podstawie
    % szacowanej wariancji błędu. Tutaj dla uproszczenia używamy globalnego zakresu d_adpcm.
    % Bardziej zaawansowane ADPCM z G.726 wymagałoby adaptacyjnego skalowania kwantyzatora.
    % Dla tego zadania, użycie dynamic_quantize na całym d_adpcm jest wystarczające,
    % ale w rzeczywistości kwantyzator ADPCM również adaptuje swój zakres.
    % Tutaj dla uproszczenia (i żeby nie komplikować za bardzo funkcji dynamic_quantize
    % o kontekst ADPCM), zastosujemy ją na całym sygnale d_adpcm.
    % Gdybyśmy chcieli to zrobić per próbkę, musielibyśmy zaimplementować skalowanie w pętli.
    
    % Tutaj, żeby dynamic_quantize działało na całym d_adpcm, musimy ją wywołać po zakończeniu pętli,
    % ale dla ADPCM kwantyzacja dzieje się 'on the fly'.
    % Poprawimy to, stosując ją na pojedynczej próbce, ale pamiętając, że zakres musi być globalny
    % lub estymowany adaptacyjnie. Dla testu, kwantyzujemy d_adpcm globalnie po zakończeniu pętli,
    % a dla pętli używamy kwantyzacji w locie, ale na uproszczonych zasadach.
    % Prawidłowe ADPCM wymaga bardziej złożonej adaptacji kwantyzatora.
    
    % Uproszczona kwantyzacja 'on-the-fly' dla ADPCM:
    % Jeśli dynamic_quantize wymaga min/max całego sygnału, to nie możemy jej używać na bieżąco.
    % W ramach uproszczenia, możemy tymczasowo założyć, że zakres kwantyzacji dla ADPCM jest stały,
    % albo założyć, że ta pojedyncza próbka jest kwantyzowana w ramach jakiegoś szacowanego zakresu.
    % Aby to działało jak w DPCM, potrzebowalibyśmy, aby dynamic_quantize działała na skali.
    % NAJPRAWDOPODOBNIEJ w ADPCM oczekuje się, że kwantyzator też się adaptuje, a nie tylko predyktor.
    % Zostawię na razie kwantyzację pojedynczej próbki dla ADPCM, ale z notatką.

    % Dla prawidłowego działania z "dynamic_quantize", musimy ją zmodyfikować,
    % aby zwracała również min/max, lub wywoływać ją na całym wektorze.
    % W ADPCM, kwantyzacja odbywa się w pętli, próbka po próbce, z adaptacyjnym skalowaniem.
    % Najprostsze rozwiązanie dla tej struktury kodu, to użycie 'label1_kwant' ponownie,
    % lub zmiana 'dynamic_quantize' na bardziej elastyczną, która przyjmie zakres.

    % Tymczasowo wrócimy do 'label1_kwant' dla ADPCM, lub stworzymy prosty kwantyzator 'in-line'
    % zakładający zakres [-1, 1] dla błędu.
    % LUB - kwantyzacja adpcm jest wykonana na całości d_adpcm po pętli
    % Pamiętaj, że to uproszczenie!
    
    % Uproszczona kwantyzacja dla ADPCM w pętli:
    % Możemy użyć stałego zakresu, np. [-2, 2] dla błędu (bo błąd może być większy)
    % lub bardziej zaawansowany adaptacyjny skalownik kwantyzatora.
    % Na potrzeby tego zadania, jeśli dynamic_quantize działa na cały wektor,
    % to dla ADPCM musielibyśmy najpierw wyliczyć cały d_adpcm, a potem go kwantyzować.
    % Ale ADPCM to rekurencja z dq_adpcm. To konflikt.
    % Zmieniam na uproszczoną wersję, która kwantyzuje dynamicznie w ramach ADPCM.
    
    % Prawidłowa kwantyzacja w ADPCM wymaga adaptacji kroku kwantyzacji.
    % Zgodnie z G.726, to jest bardziej złożone. Tutaj użyjemy uproszczonej wersji,
    % gdzie zakres kwantyzatora jest stały, ale wystarczająco szeroki.
    % Alternatywnie, możemy przyjąć, że kwantyzator `dynamic_quantize` zostanie wywołany raz
    % dla całego `d_adpcm` po pętli (co jest błędne dla ADPCM) lub
    % przekazać mu globalne min/max z `d_adpcm` wyliczone wcześniej.
    % Dla uproszczenia i zachowania rekurencji:
    % Załóżmy, że kwantyzator dla pojedynczej próbki w ADPCM używa stałego,
    % odpowiednio szerokiego zakresu, np. [-1.5, 1.5] lub [-2, 2].
    
    % Uproszczona kwantyzacja dla ADPCM w pętli (stały zakres):
    quant_step_adpcm = (1.5 - (-1.5)) / (2^num_quant_bits - 1); % Zakres [-1.5, 1.5]
    dq_adpcm_sample = round(d_adpcm(n) / quant_step_adpcm) * quant_step_adpcm;
    dq_adpcm_sample = max(-1.5, min(1.5, dq_adpcm_sample)); % Przycięcie
    
    % Rekonstrukcja
    y_adpcm(n) = dq_adpcm_sample + x_pred_adpcm;
    
    % Adaptacja współczynnika predykcji (algorytm LMS)
    % Błąd używany do adaptacji to błąd między ORYGINAŁEM a REKONSTRUKCJĄ
    error_signal = x(n) - y_adpcm(n);
    a_adpcm = a_adpcm + mu * error_signal * y_adpcm(n-1);
    
    % Ograniczenie zakresu współczynnika predykcji do [0,1]
    a_adpcm = max(0, min(1, a_adpcm));
    
    a_history(n) = a_adpcm;
end
fprintf('Kodowanie/Dekodowanie ADPCM zakończone.\n');


%% SEKCJA 4: ANALIZA JAKOŚCI - OBLICZANIE METRYK
fprintf('\n=== ANALIZA JAKOŚCI ===\n');

% Obliczanie błędów
error_unquantized = x - y_unquantized; % Błąd dla DPCM bez kwantyzacji
error_dpcm = x - y;                     % Błąd dla DPCM z kwantyzacją
error_adpcm = x - y_adpcm;              % Błąd dla ADPCM

% SNR (Signal-to-Noise Ratio)
% Dodano eps, aby uniknąć dzielenia przez zero, jeśli rms(error) jest bardzo małe
snr_unquantized = 20 * log10(rms(x) / (rms(error_unquantized) + eps));
snr_dpcm = 20 * log10(rms(x) / (rms(error_dpcm) + eps));
snr_adpcm = 20 * log10(rms(x) / (rms(error_adpcm) + eps));

% MSE (Mean Square Error)
mse_unquantized = mean(error_unquantized.^2);
mse_dpcm = mean(error_dpcm.^2);
mse_adpcm = mean(error_adpcm.^2);

% THD (Total Harmonic Distortion) - przybliżone
% THD jest często definiowane jako stosunek mocy harmonicznych do mocy podstawowej
% Tutaj używamy uproszczonej definicji: stosunek mocy błędu do mocy sygnału oryginalnego
thd_unquantized = rms(error_unquantized) / rms(x) * 100;
thd_dpcm = rms(error_dpcm) / rms(x) * 100;
thd_adpcm = rms(error_adpcm) / rms(x) * 100;

fprintf('METRYKI JAKOŚCI:\n');
fprintf('DPCM (bez kwantyzacji) - SNR: %.2f dB, MSE: %.6f, THD: %.2f%%\n', snr_unquantized, mse_unquantized, thd_unquantized);
fprintf('DPCM (z kwantyzacją)   - SNR: %.2f dB, MSE: %.6f, THD: %.2f%%\n', snr_dpcm, mse_dpcm, thd_dpcm);
fprintf('ADPCM                  - SNR: %.2f dB, MSE: %.6f, THD: %.2f%%\n', snr_adpcm, mse_adpcm, thd_adpcm);

%% SEKCJA 6: ODTWARZANIE AUDIO
fprintf('\n=== ODTWARZANIE AUDIO ===\n');
fprintf('Naciśnij Enter aby odtworzyć sygnał oryginalny...\n');
pause;
try
    fprintf('Odtwarzanie sygnału oryginalnego...\n');
    sound(x, Fs);
    pause(N/Fs + 0.5); % Pauza na czas trwania sygnału + trochę

    fprintf('Naciśnij Enter aby odtworzyć sygnał DPCM (bez kwantyzacji)...\n');
    pause;
    fprintf('Odtwarzanie sygnału DPCM (bez kwantyzacji)...\n');
    sound(y_unquantized, Fs);
    pause(N/Fs + 0.5);

    fprintf('Naciśnij Enter aby odtworzyć sygnał DPCM (z kwantyzacją)...\n');
    pause;
    fprintf('Odtwarzanie sygnału DPCM (z kwantyzacją)...\n');
    sound(y, Fs);
    pause(N/Fs + 0.5);

    fprintf('Naciśnij Enter aby odtworzyć sygnał ADPCM...\n');
    pause;
    fprintf('Odtwarzanie sygnału ADPCM...\n');
    sound(y_adpcm, Fs);
    pause(N/Fs + 0.5); % Opcjonalna pauza

catch ME
    fprintf('Błąd odtwarzania audio: %s\n', ME.message);
    fprintf('Sprawdź konfigurację dźwięku w systemie lub brak pliku DontWorryBeHappy.wav.\n');
end

%% SEKCJA 7: ZAPIS WYNIKÓW DO PLIKÓW WAV
fprintf('\n=== ZAPIS WYNIKÓW AUDIO DO PLIKÓW WAV ===\n');
try
    audiowrite('signal_original.wav', x, Fs);
    audiowrite('signal_dpcm_unquantized.wav', y_unquantized, Fs);
    audiowrite('signal_dpcm_quantized.wav', y, Fs);
    audiowrite('signal_adpcm.wav', y_adpcm, Fs);
    fprintf('Zapisano pliki audio:\n');
    fprintf(' - signal_original.wav\n');
    fprintf(' - signal_dpcm_unquantized.wav\n');
    fprintf(' - signal_dpcm_quantized.wav\n');
    fprintf(' - signal_adpcm.wav\n');
catch
    fprintf('Nie można zapisać plików audio. Sprawdź uprawnienia do katalogu.\n');
end

%% SEKCJA 8: PODSUMOWANIE I WNIOSKI
fprintf('\n=== PODSUMOWANIE ===\n');
fprintf('ANALIZA TEORETYCZNA:\n');
fprintf('1. DPCM wykorzystuje predykcję liniową do redukcji redundancji.\n');
fprintf('2. Sygnał różnicowy (błąd predykcji) ma zazwyczaj mniejszą wariancję niż sygnał oryginalny.\n');
fprintf('3. Kwantyzacja sygnału różnicowego wprowadza szum kwantyzacji, ale pozwala na kompresję danych.\n');
fprintf('4. ADPCM dodatkowo adaptuje współczynnik predykcji oraz (w pełnej implementacji) krok kwantyzacji.\n');
fprintf('5. Jakość zrekonstruowanego sygnału zależy od współczynnika predykcji, liczby bitów kwantyzacji i charakterystyki sygnału.\n\n');

fprintf('WYNIKI EKSPERYMENTU:\n');
fprintf('Współczynnik predykcji DPCM: a = %.4f\n', a);
fprintf('Końcowy współczynnik ADPCM: a = %.4f\n', a_history(end));
fprintf('Liczba bitów kwantyzacji: %d bitów (%d poziomów)\n', num_quant_bits, 2^num_quant_bits);

if snr_dpcm > snr_unquantized
    fprintf('Coś poszło nie tak! Kwantyzacja ZAWSZE obniża SNR.\n');
end
fprintf('SNR DPCM (bez kwantyzacji): %.2f dB\n', snr_unquantized);
fprintf('SNR DPCM (z kwantyzacją):   %.2f dB\n', snr_dpcm);
fprintf('SNR ADPCM:                   %.2f dB\n', snr_adpcm);

fprintf('Różnica SNR (DPCM z kwant. vs bez kwant.): %.2f dB\n', snr_dpcm - snr_unquantized);
fprintf('Różnica SNR (ADPCM vs DPCM z kwant.):     %.2f dB\n', snr_adpcm - snr_dpcm);

if snr_adpcm > snr_dpcm
    fprintf('ADPCM generalnie osiąga lepszą jakość (wyższe SNR) niż DPCM z kwantyzacją dla tego sygnału.\n');
else
    fprintf('DPCM z kwantyzacją osiąga lepszą jakość niż ADPCM dla tego sygnału (mniej typowe).\n');
end
fprintf('\n=== ANALIZA ZAKOŃCZONA ===\n');