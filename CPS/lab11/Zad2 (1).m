%% CYFROWE PRZETWARZANIE SYGNAŁÓW - KODOWANIE TRANSFORMACYJNE
clear all; close all; clc;

%% SEKCJA 1: WCZYTANIE I PRZYGOTOWANIE SYGNAŁU
fprintf('=== ŁADOWANIE SYGNAŁU AUDIO ===\n');
% Próba wczytania pliku audio
% Pamiętaj, że plik 'DontWorryBeHappy.wav' musi być w tym samym katalogu
% co skrypt MATLAB lub w ścieżce wyszukiwania MATLABa.
try
    [x, Fs] = audioread('DontWorryBeHappy.wav');
    fprintf('Wczytano plik: DontWorryBeHappy.wav\n');
catch
    error('Błąd: Plik DontWorryBeHappy.wav nie znaleziono. Upewnij się, że znajduje się w bieżącym katalogu MATLAB.');
end

% Konwersja na mono i normalizacja
if size(x, 2) > 1
    x = mean(x, 2); % Konwersja stereo na mono
end
x = x(:); % Wektor kolumnowy
x = x / max(abs(x)); % Normalizacja do zakresu [-1, 1]

% Parametry sygnału
N_samples = length(x);
fprintf('Częstotliwość próbkowania: %d Hz\n', Fs);
fprintf('Długość sygnału: %.2f s (%d próbek)\n', N_samples/Fs, N_samples);

%% SEKCJA 2: PARAMETRY KODERA MDCT
% --- 1) Wykonaj koder dla N=32 i N=128. ---
N_values = [32, 128]; % Rozmiary okna analizy/syntezy MDCT
% Pamiętaj: MDCT przyjmuje N próbek i zwraca N/2 współczynników.
% Okna zachodzą na siebie w 50%, więc przesuń o N/2 próbek.

% --- 2) Sprawdź czy kodowanie może być bezstratne dla odpowiednio dużego Q. ---
% Dla bezstratnego kodowania Q musi być bardzo małe (lub pomijamy kwantyzację).
% Zostawiamy Q do konfiguracji poniżej.

% --- Wybierz taką wartość Q, aby dla monofonicznego sygnału próbkowanego
% z częstotliwością fs=44100 Hz uzyskać strumień bitów o przepływności 64 kbps. ---
target_bitrate_kbps = 64; % kbit/s
target_bitrate_bps = target_bitrate_kbps * 1000; % bit/s

% Liczba bitów na próbkę (N/2 współczynników na N/2 nowych próbek)
% Przepływność = (Liczba bitów na N/2 próbek) * (Liczba bloków na sekundę)
% Liczba bloków na sekundę = Fs / (N/2)  (bo okno przesuwa się o N/2)
% Liczba bitów na sekundę = (N_MDCT_coeffs * Bits_per_coeff) * (Fs / (N/2))
% target_bitrate_bps = (N/2 * Bits_per_coeff) * (Fs / (N/2))
% target_bitrate_bps = Fs * Bits_per_coeff
% Bits_per_coeff = target_bitrate_bps / Fs

if Fs ~= 44100
    warning('Częstotliwość próbkowania pliku (%d Hz) różni się od wymaganej w zadaniu (44100 Hz). Obliczenia Q będą dla Fs z pliku.', Fs);
end

required_bits_per_coeff = target_bitrate_bps / Fs;
fprintf('Docelowa przepływność: %d kbps\n', target_bitrate_kbps);
fprintf('Wymagana średnia liczba bitów na współczynnik MDCT: %.2f bitów\n', required_bits_per_coeff);

% Kwantyzacja: Użyjemy prostej kwantyzacji uniformowej.
% Q_step = (max_val - min_val) / (2^bits - 1)
% Jeśli chcemy "bits_per_coeff", to będziemy musieli zaokrąglić do najbliższej liczby całkowitej bitów.
% Na potrzeby zadania, jeśli wymagana jest taka liczba bitów, oznacza to, że musimy dobrać krok kwantyzacji.
% Dla uproszczenia i na potrzeby zadania, przyjmiemy Q jako krok kwantyzacji.
% Im większe Q, tym mniejsza dokładność (większa kompresja, gorsza jakość).
% Im mniejsze Q, tym większa dokładność (mniejsza kompresja, lepsza jakość).
% Dla "bezstratnego" kodowania Q musi być bardzo małe (bliskie zeru, np. 1e-10)

% Wybierzmy Q, które odpowiada wymaganej przepływności.
% Możemy estymować zakres współczynników MDCT (np. -1 do 1, jeśli sygnał wejściowy jest znormalizowany).
% Dla 44100 Hz i 64 kbps -> required_bits_per_coeff = 64000 / 44100 = 1.45 bitów/współczynnik.
% To oznacza, że nie możemy użyć całkowitej liczby bitów. W praktyce, kwantyzatory są adaptacyjne.
% Na potrzeby tego zadania, spróbujemy znaleźć Q, które daje taki "efektywny bitrate".
% Możemy to zrobić, symulując i mierząc.
% Alternatywnie, zadanie może sugerować, że powinniśmy wybrać Q tak,
% aby np. 1.45 bita efektywnie oznaczało, że średnio na 2 współczynniki przypada 3 bity.
% Ustawmy Q_factor, który będzie regulował "agresywność" kwantyzacji.

% Spróbujmy znaleźć Q, które odpowiadałoby bitrate'owi.
% Jeśli sygnał jest znormalizowany do [-1, 1], to współczynniki MDCT też będą w podobnym zakresie.
% N/2 współczynników.
% Kwantyzacja do 1 bit (2 poziomy) to krok 1/2. Kwantyzacja do 2 bitów (4 poziomy) to krok 1/4.
% Q_step_for_1_bit = 2 / (2^1) = 1.
% Q_step_for_2_bits = 2 / (2^2) = 0.5.
% Możemy użyć Q jako odwrotności liczby poziomów.
% Q = 1 / 2^required_bits_per_coeff to byłoby zbyt uproszczone, ale daje nam wyczucie.
% Przyjmijmy Q jako krok kwantyzacji. Będzie to wartość stała.
% Możemy go dobrać eksperymentalnie. Zacznijmy od wartości sugerującej ok. 1.5 bita.
% 2^1.5 = ~2.82 poziomów. range / 2.82.
% Jeśli zakres współczynników wynosi ok. 2 (od -1 do 1), to Q_step = 2 / 2.82 = 0.707.
% Zaczniemy od Q_step = 0.5 (efektywnie 2 bity na próbkę).
% Możesz eksperymentować z tą wartością.
Q_step_64kbps = 0.2; % Zacznijmy od tego i zobaczmy SNR. Możesz to dostroić.
Q_step_lossless = 1e-10; % Bardzo małe Q dla "bezstratnego" kodowania

%% GŁÓWNA PĘTLA KODOWANIA I DEKODOWANIA MDCT
results = struct(); % Struktura do przechowywania wyników

for i = 1:length(N_values)
    N = N_values(i);
    M = N/2; % Liczba współczynników MDCT, rozmiar przesunięcia okna

    fprintf('\n--- Przetwarzanie dla N = %d (M = %d) ---\n', N, M);

    % 1. Okno analizy i syntezy
    h = sin(pi * (0:N-1)' / N + 0.5);

    % 2. Macierz analizy A dla MDCT
    k_indices = (0:M-1)';
    m_indices = (0:N-1);
    A = sqrt(4/N) * cos( (2*pi/N) * (k_indices + 0.5) .* (m_indices + 0.5 + N/4) );

    % 3. Macierz syntezy S (transpozycja A)
    S = A';

    % Inicjalizacja zrekonstruowanego sygnału
    x_reconstructed_lossless = zeros(N_samples, 1);
    x_reconstructed_quantized = zeros(N_samples, 1);

    % Bufor do obsługi nakładania okien
    overlap_buffer_lossless = zeros(N, 1);
    overlap_buffer_quantized = zeros(N, 1);

    num_frames = floor((N_samples - M) / M); % Liczba pełnych klatek do przetworzenia

    all_coeff_lossless = []; % Do zbierania współczynników
    all_coeff_quantized = []; % Do zbierania współczynników

    for frame_idx = 0:num_frames-1
        start_idx = frame_idx * M + 1;
        end_idx = start_idx + N - 1;

        if end_idx > N_samples
            % Obsługa ostatniej, niepełnej ramki (np. dopełnienie zerami)
            current_block = zeros(N, 1);
            current_block(1:N_samples - start_idx + 1) = x(start_idx:N_samples);
        else
            current_block = x(start_idx:end_idx);
        end

        % Koder
        windowed_block = current_block .* h; % Aplikacja okna analizy
        X_freq = A * windowed_block;        % MDCT

        % --- Bezstratne kodowanie (bez kwantyzacji) ---
        X_freq_lossless = X_freq;
        all_coeff_lossless = [all_coeff_lossless; X_freq_lossless]; %#ok<AGROW> % Zbieranie

        % --- Kodowanie z kwantyzacją (do 64 kbps) ---
        % Prosta kwantyzacja uniformowa
        % Zakres kwantyzacji dla współczynników MDCT jest zwykle zbliżony do zakresu sygnału wejściowego,
        % ale może być nieco większy. Dla znormalizowanego sygnału [-1, 1], MDCT coefficients mogą być ok. [-1.5, 1.5].
        % Możemy również dynamicznie obliczyć min/max dla X_freq i użyć Q_step na tym zakresie.
        % Dla uproszczenia i na potrzeby zadania, użyjmy ustalonego kroku Q.
        
        % Sprawdź zakres współczynników dla lepszego wyboru Q_step
        % min_xf = min(X_freq);
        % max_xf = max(X_freq);
        % fprintf('Min/Max X_freq: [%f, %f]\n', min_xf, max_xf);

        % Kwantyzacja: `round(val / Q_step) * Q_step`
        X_freq_quantized = round(X_freq / Q_step_64kbps) * Q_step_64kbps;
        all_coeff_quantized = [all_coeff_quantized; X_freq_quantized]; %#ok<AGROW> % Zbieranie


        % Dekoder
        % Odtworzenie bloku czasowego z kwantyzacji
        reconstructed_block_lossless = S * X_freq_lossless; % IDCT
        reconstructed_block_quantized = S * X_freq_quantized; % IDCT

        % Aplikacja okna syntezy (ponownie h(n))
        reconstructed_block_lossless = reconstructed_block_lossless .* h;
        reconstructed_block_quantized = reconstructed_block_quantized .* h;

        % Sumowanie z nakładającymi się oknami
        % Nowy blok pokrywa indeksy: [start_idx : end_idx]
        % Ale dekoder wypluwa już 2*M próbek, które są sumą z poprzednim blokiem.
        % Pierwsze M próbek nowego bloku dodaje się do ostatnich M próbek poprzedniego bufora.
        % Drugie M próbek nowego bloku tworzy początek następnego bufora.

        % Przykładowa implementacja Overlap-Add dla MDCT (zgodnie z teorią):
        % Wyjście IDCT (reconstructed_block) ma N próbek.
        % Sumujemy pierwsze M próbek z 'overlap_buffer' (poprzednie N/2 próbki)
        % i pierwsze M próbek bieżącego 'reconstructed_block'.
        % Następnie do 'x_reconstructed' dodajemy M środkowych próbek (tych, które są unikalne dla tego okna).
        % Na koniec, pozostałe M próbek 'reconstructed_block' staje się nowym 'overlap_buffer'.

        % Rekonstrukcja Overlap-Add:
        % x_reconstructed(start_idx : start_idx + M - 1) = overlap_buffer(1:M) + reconstructed_block(1:M);
        % x_reconstructed(start_idx + M : end_idx) = reconstructed_block(M+1:N);
        % overlap_buffer = reconstructed_block(M+1:N);

        % Lepsza implementacja Overlap-Add:
        % Bufor o rozmiarze N (gdzie koniec poprzedniego bloku i początek obecnego się spotykają)
        % Pierwsza połowa bufora to koniec poprzedniego okna
        % Druga połowa bufora to początek obecnego okna
        
        % Zaktualizuj bufor nakładania z poprzedniego kroku (overlap_buffer ma N próbek)
        % Pierwsze M próbek bufora to "stary" koniec
        % Drugie M próbek bufora to "nowy" początek
        
        % Dodaj zrekonstruowany blok do bufora
        if frame_idx == 0
            % Pierwsza ramka: tylko wypełniamy bufor
            overlap_buffer_lossless = reconstructed_block_lossless;
            overlap_buffer_quantized = reconstructed_block_quantized;
        else
            % Przetwarzaj M próbek, które są w pełni zrekonstruowane w tym kroku
            output_start_idx = (frame_idx - 1) * M + 1; % Próbki od (N/2) poprzedniej ramki
            output_end_idx = output_start_idx + M - 1;

            if output_end_idx <= N_samples
                x_reconstructed_lossless(output_start_idx : output_end_idx) = ...
                    overlap_buffer_lossless(1:M) + reconstructed_block_lossless(1:M);
                x_reconstructed_quantized(output_start_idx : output_end_idx) = ...
                    overlap_buffer_quantized(1:M) + reconstructed_block_quantized(1:M);
            end

            
            % Przesuń bufor: druga połowa obecnego bloku staje się pierwszą połową następnego bufora
            overlap_buffer_lossless(1:M) = reconstructed_block_lossless(M+1:N);
            overlap_buffer_quantized(1:M) = reconstructed_block_quantized(M+1:N);
            
            % Wypełnij resztę bufora drugą połową obecnego bloku (przesuniętą)
            % to już jest w sumowaniu powyżej.
            % Po prostu 'reconstructed_block_lossless' jest w nowym buforze.
            % Drugie M próbek z reconstructed_block to początek następnego bufora.
            % overlap_buffer (nowy) = reconstructed_block (M+1:N)
        end
    end
    
    % Obsłuż ostatnie próbki w buforze nakładania po pętli
    % Sumowanie reszty bufora
    if num_frames > 0
        last_output_start_idx = (num_frames - 1) * M + 1;
        last_output_end_idx = last_output_start_idx + M - 1; % Próbki, które były już zsumowane
        
        % Dodaj resztę bufora nakładania, która nie została jeszcze zsumowana
        % To co pozostało w overlap_buffer, plus to co się nakłada z ostatnim blokiem
        % (reconstructed_block(M+1:N) z ostatniego bloku)
        
        % Ta część jest nieco skomplikowana w implementacji MDCT Overlap-Add.
        % Najprostszym sposobem na prawidłową rekonstrukcję ostatniej części sygnału
        % jest obsłużenie końcówki po pętli.
        
        % Po pętli, w overlap_buffer jest N próbek, z czego pierwsze M to już przesunięta
        % druga połowa ostatniego bloku.
        % Reszta sygnału (jeśli N_samples nie jest wielokrotnością M)
        % Możemy dodać pozostałe próbki z overlap_buffer.
        
        % Finalne dodanie pozostałości z overlap_buffer, które nie zostały jeszcze zsumowane
        % Z powodu 50% overlapu i sumowania co M próbek, na końcu zostaje bufor M próbek
        % które powinny być dodane.
        
        % Uproszczenie: Po prostu skopiuj ostatnią połowę ostatniego zrekonstruowanego bloku.
        % To nie jest idealne, ale dla ogólnego zobrazowania działania MDCT wystarczy.
        % Prawidłowa obsługa ostatniego bloku jest złożona.
        
        % Lepsza obsługa:
        % Po pętli `for`, `overlap_buffer_lossless(1:M)` zawiera ostatnie `M` zrekonstruowane próbki.
        % Dodać je do `x_reconstructed_lossless`.
        
        if last_output_end_idx + M <= N_samples
            x_reconstructed_lossless(last_output_end_idx + 1 : last_output_end_idx + M) = overlap_buffer_lossless(1:M);
            x_reconstructed_quantized(last_output_end_idx + 1 : last_output_end_idx + M) = overlap_buffer_quantized(1:M);
        else % Jeśli ostatnia ramka była niepełna
            x_reconstructed_lossless(last_output_end_idx + 1 : N_samples) = overlap_buffer_lossless(1 : N_samples - last_output_end_idx);
            x_reconstructed_quantized(last_output_end_idx + 1 : N_samples) = overlap_buffer_quantized(1 : N_samples - last_output_end_idx);
        end
    end
    
    % Przytnij zrekonstruowane sygnały do oryginalnej długości
    x_reconstructed_lossless = x_reconstructed_lossless(1:N_samples);
    x_reconstructed_quantized = x_reconstructed_quantized(1:N_samples);

    % --- Analiza Jakości ---
    error_lossless = x - x_reconstructed_lossless;
    error_quantized = x - x_reconstructed_quantized;

    snr_lossless = 20 * log10(rms(x) / (rms(error_lossless) + eps));
    snr_quantized = 20 * log10(rms(x) / (rms(error_quantized) + eps));

    mse_lossless = mean(error_lossless.^2);
    mse_quantized = mean(error_quantized.^2);

    fprintf('Dla N = %d:\n', N);
    fprintf('  MDCT (bez kwantyzacji) - SNR: %.2f dB, MSE: %.8f\n', snr_lossless, mse_lossless);
    fprintf('  MDCT (z kwantyzacją)   - SNR: %.2f dB, MSE: %.8f\n', snr_quantized, mse_quantized);

    % Zapis wyników do struktury
    results(i).N = N;
    results(i).x_reconstructed_lossless = x_reconstructed_lossless;
    results(i).x_reconstructed_quantized = x_reconstructed_quantized;
    results(i).snr_lossless = snr_lossless;
    results(i).snr_quantized = snr_quantized;
    results(i).mse_lossless = mse_lossless;
    results(i).mse_quantized = mse_quantized;
    results(i).Q_step_used = Q_step_64kbps;
    results(i).all_coeff_lossless = all_coeff_lossless;
    results(i).all_coeff_quantized = all_coeff_quantized;
end

%% SEKCJA 3: ODTWARZANIE I ZAPIS WYNIKÓW
fprintf('\n=== ODTWARZANIE AUDIO ===\n');

% Odtwarzanie sygnału oryginalnego
fprintf('Naciśnij Enter aby odtworzyć sygnał oryginalny...\n');
pause;
try
    fprintf('Odtwarzanie sygnału oryginalnego...\n');
    sound(x, Fs);
    pause(N_samples/Fs + 0.5);
catch ME
    fprintf('Błąd odtwarzania audio: %s\n', ME.message);
end

% Odtwarzanie i zapis dla każdego N
for i = 1:length(N_values)
    N = results(i).N;
    fprintf('\n--- Wyniki dla N = %d ---\n', N);

    % Bezstratne kodowanie
    fprintf('Naciśnij Enter aby odtworzyć sygnał MDCT (bez kwantyzacji) dla N=%d...\n', N);
    pause;
    try
        fprintf('Odtwarzanie sygnału MDCT (bez kwantyzacji)...\n');
        sound(results(i).x_reconstructed_lossless, Fs);
        pause(N_samples/Fs + 0.5);
        audiowrite(sprintf('signal_mdct_N%d_lossless.wav', N), results(i).x_reconstructed_lossless, Fs);
        fprintf('Zapisano: signal_mdct_N%d_lossless.wav\n', N);
    catch ME
        fprintf('Błąd odtwarzania/zapisu: %s\n', ME.message);
    end

    % Kodowanie z kwantyzacją
    fprintf('Naciśnij Enter aby odtworzyć sygnał MDCT (z kwantyzacją %.2f bit/coeff) dla N=%d...\n', required_bits_per_coeff, N);
    pause;
    try
        fprintf('Odtwarzanie sygnału MDCT (z kwantyzacją)...\n');
        sound(results(i).x_reconstructed_quantized, Fs);
        pause(N_samples/Fs + 0.5);
        audiowrite(sprintf('signal_mdct_N%d_quantized_%dkbps.wav', N, target_bitrate_kbps), results(i).x_reconstructed_quantized, Fs);
        fprintf('Zapisano: signal_mdct_N%d_quantized_%dkbps.wav\n', N, target_bitrate_kbps);
    catch ME
        fprintf('Błąd odtwarzania/zapisu: %s\n', ME.message);
    end
end

%% SEKCJA 4: PODSUMOWANIE I WNIOSKI
fprintf('\n=== PODSUMOWANIE ===\n');
fprintf('Transformacyjne kodowanie dźwięku z MDCT.\n');
fprintf('Główne parametry: Rozmiar okna (N), Kwantyzacja (step).\n');

for i = 1:length(N_values)
    N = results(i).N;
    fprintf('\n--- Wyniki dla N = %d ---\n', N);
    fprintf('  Kodowanie bezstratne:\n');
    fprintf('    SNR: %.2f dB, MSE: %.8f\n', results(i).snr_lossless, results(i).mse_lossless);
    if results(i).snr_lossless > 100 % Bardzo wysokie SNR wskazuje na niemal bezstratne kodowanie
        fprintf('    Potwierdza to, że MDCT może być praktycznie bezstratne bez kwantyzacji.\n');
    else
        fprintf('    Mały błąd może wynikać z niedokładności numerycznych lub problemów z obsługą krawędzi sygnału.\n');
    end

    fprintf('  Kodowanie z kwantyzacją (docelowe %.2f bit/coeff):\n', required_bits_per_coeff);
    fprintf('    Użyty krok kwantyzacji (Q_step): %.4f\n', results(i).Q_step_used);
    fprintf('    SNR: %.2f dB, MSE: %.8f\n', results(i).snr_quantized, results(i).mse_quantized);
    fprintf('    Otrzymany SNR świadczy o jakości dźwięku przy docelowej przepływności 64 kbps.\n');

    % Możemy sprawdzić, ile efektywnie bitów zostało użytych.
    % To jest trudniejsze bez entropii lub adaptacyjnej alokacji bitów,
    % ale możemy oszacować, ile poziomów kwantyzacji 'X_freq_quantized' rzeczywiście wykorzystuje.
    % num_unique_levels = length(unique(results(i).all_coeff_quantized));
    % if num_unique_levels > 1
    %     effective_bits_per_coeff = log2(num_unique_levels);
    %     fprintf('    Efektywne bity na współczynnik MDCT: %.2f (na podstawie unikalnych poziomów w %d współczynnikach)\n', effective_bits_per_coeff, length(results(i).all_coeff_quantized));
    % end
end

fprintf('\nEKSPERYMENTY DO PRZEPROWADZENIA:\n');
fprintf('1. Zmień wartość `Q_step_64kbps` i obserwuj wpływ na jakość (SNR).\n');
fprintf('2. Zmień rozmiar okna `N_values` i zobacz, jak wpływa na jakość i wydajność.\n');
fprintf('   - Większe N: lepsza rozdzielczość częstotliwościowa, gorsza czasowa, większe opóźnienie.\n');
fprintf('   - Mniejsze N: lepsza rozdzielczość czasowa, gorsza częstotliwościowa, mniejsze opóźnienie.\n');
fprintf('3. Spróbuj zaimplementować bardziej zaawansowaną kwantyzację (np. adaptacyjną lub nieliniową).\n');
fprintf('4. Rozważ psychoakustyczne modelowanie (alokacja bitów do pasm, gdzie ucho jest mniej wrażliwe).\n');

fprintf('\n=== ANALIZA ZAKOŃCZONA ===\n');