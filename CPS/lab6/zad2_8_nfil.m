%% KOMPLETNY DEKODER DTMF
clc; clear; close all;

%% 1. Wczytanie sygnału DTMF
filename = 's.wav'; % Tutaj wstaw właściwy plik .wav
[s, fs] = audioread(filename);

% Wyświetlenie oryginalnego spektrogramu
figure('Name', 'Spektrogram oryginalnego sygnału');
spectrogram(s, 4096, 4096-512, 0:5:2000, fs, 'yaxis');
title(['Spektrogram sygnalu DTMF - ' filename]);

%% 2. Definicja parametrów DTMF
% Częstotliwości DTMF
row_freq = [697, 770, 852, 941];  % częstotliwości wierszy
col_freq = [1209, 1336, 1477];    % częstotliwości kolumn
all_freq = [row_freq, col_freq];  % wszystkie częstotliwości

% Mapowanie częstotliwości na znaki
dtmf_keys = ['1', '2', '3';
             '4', '5', '6';
             '7', '8', '9';
             '*', '0', '#'];

%% 3. FILTRACJA PASMOWA (600-1500 Hz)
disp('=== FILTRACJA PASMOWA ===');

% Projektowanie filtru Butterwortha 4-go rzędu
low_cutoff = 697 / (fs/2);
high_cutoff = 1477 / (fs/2);
[b, a] = butter(4, [low_cutoff, high_cutoff], 'bandpass');

% Filtracja sygnału
s_filtered = filter(b, a, s);

% Obliczenie opóźnienia z zabezpieczeniami
try
    delay = round(mean(grpdelay(b, a)));
    if ~isfinite(delay) || isempty(delay)
        delay = 0;
    end
catch
    delay = 0;
end

% Zabezpieczenie przed nieprawidłowymi wartościami delay
if ~isnumeric(delay) || ~isscalar(delay) || ~isfinite(delay)
    delay = 0;
end

% Kompensacja opóźnienia z pełnym zabezpieczeniem
try
    if delay > 0
        % Ograniczenie delay do maksymalnej możliwej wartości
        delay = min(delay, length(s_filtered)-1);
        s_filtered = [s_filtered(delay+1:end); zeros(delay, 1)];
    elseif delay < 0
        % Obsługa ujemnego opóźnienia (teoretycznie nie powinno wystąpić)
        delay = abs(min(delay, length(s_filtered)-1));
        s_filtered = [zeros(delay, 1); s_filtered(1:end-delay)];
    end
catch ME
    warning('Problem z kompensacją opóźnienia: %s', ME.message);
    % W przypadku błędu używamy filtfilt który nie wymaga kompensacji
    s_filtered = filtfilt(b, a, s);
end

% Wizualizacja efektów filtracji
figure('Name', 'Porównanie przed i po filtracji');
subplot(2,1,1);
spectrogram(s, 4096, 4096-512, 0:5:2000, fs, 'yaxis');
title('Przed filtracją');

subplot(2,1,2);
spectrogram(s_filtered, 4096, 4096-512, 0:5:2000, fs, 'yaxis');
title('Po filtracji');

figure('Name', 'Porównanie sygnałów w dziedzinie czasu');
t = (0:length(s)-1)/fs;
plot(t, s, 'b', t, s_filtered, 'r');
legend('Oryginalny', 'Przefiltrowany');
title('Porównanie sygnałów w dziedzinie czasu');
xlabel('Czas [s]'); ylabel('Amplituda');

%% 4. DEKODOWANIE ALGORYTMEM GOERTZELA
disp('=== ALGORYTM GOERTZELA ===');

% Parametry analizy
window_size = round(0.04 * fs);  % 40 ms okno (optymalne dla DTMF)
window_step = round(0.01 * fs);  % 10 ms krok
num_windows = floor((length(s_filtered) - window_size) / window_step) + 1;

% Prekomputacja współczynników Goertzela
row_coeffs = zeros(1, length(row_freq));
col_coeffs = zeros(1, length(col_freq));

for i = 1:length(row_freq)
    k = round(row_freq(i) * window_size / fs);
    omega = 2 * pi * k / window_size;
    row_coeffs(i) = 2 * cos(omega);
end

for i = 1:length(col_freq)
    k = round(col_freq(i) * window_size / fs);
    omega = 2 * pi * k / window_size;
    col_coeffs(i) = 2 * cos(omega);
end

% Przygotowanie macierzy do przechowywania wyników analizy Goertzela
freq_energies = zeros(num_windows, length(all_freq));
window_times = zeros(1, num_windows);

% Analiza okienkowa z algorytmem Goertzela
for w = 1:num_windows
    start_idx = (w-1)*window_step + 1;
    end_idx = min(start_idx + window_size - 1, length(s_filtered));
    
    % Zapewnienie stałej długości okna
    if end_idx - start_idx + 1 < window_size
        window_data = [s_filtered(start_idx:end_idx); zeros(window_size - (end_idx - start_idx + 1), 1)];
    else
        window_data = s_filtered(start_idx:end_idx);
    end
    
    % Zapisz czas środka okna
    window_times(w) = (start_idx + window_size/2)/fs;
    
    % Obliczenie energii dla częstotliwości wierszy
    for i = 1:length(row_freq)
        q1 = 0; q2 = 0;
        for n = 1:length(window_data)
            q0 = window_data(n) + row_coeffs(i)*q1 - q2;
            q2 = q1;
            q1 = q0;
        end
        freq_energies(w, i) = sqrt(q1^2 + q2^2 - q1*q2*row_coeffs(i));
    end
    
    % Obliczenie energii dla częstotliwości kolumn
    for i = 1:length(col_freq)
        q1 = 0; q2 = 0;
        for n = 1:length(window_data)
            q0 = window_data(n) + col_coeffs(i)*q1 - q2;
            q2 = q1;
            q1 = q0;
        end
        freq_energies(w, i+length(row_freq)) = sqrt(q1^2 + q2^2 - q1*q2*col_coeffs(i));
    end
    
    % Normalizacja energii w każdym oknie
    max_energy = max(freq_energies(w, :));
    if max_energy > 0
        freq_energies(w, :) = freq_energies(w, :) / max_energy;
    end
end

% Wizualizacja wyników algorytmu Goertzela
figure('Name', 'Wyniki algorytmu Goertzela');
imagesc(window_times, 1:length(all_freq), freq_energies');
colorbar;
title('Energie częstotliwości DTMF w czasie');
xlabel('Czas [s]');
ylabel('Indeks częstotliwości');
set(gca, 'YTick', 1:length(all_freq), 'YTickLabel', num2str(all_freq'));

% Wybierz kilka okien czasowych i pokaż spektrum częstotliwości
sample_windows = round(linspace(1, num_windows, min(5, num_windows)));
figure('Name', 'Przykładowe spektra Goertzela');
for i = 1:length(sample_windows)
    subplot(length(sample_windows), 1, i);
    bar(all_freq, freq_energies(sample_windows(i), :));
    title(['Energia częstotliwości w czasie t = ' num2str(window_times(sample_windows(i)), '%.2f') ' s']);
    xlabel('Częstotliwość [Hz]');
    ylabel('Znormalizowana energia');
    ylim([0 1.1]);
end

%% 5. UPROSZCZONY ALGORYTM DECYZYJNY
disp('=== UPROSZCZONY ALGORYTM DECYZYJNY ===');

% Zoptymalizowane parametry
threshold = 0.25;           % Niższy próg detekcji
min_contrast = 1.0;         % Niższy wymagany kontrast
min_duration = 0.02;        % Minimalny czas trwania tonu (30 ms)

% Uproszczona segmentacja
energy_sum = sum(freq_energies.^2, 2);  % Suma energii dla każdego okna
is_active = energy_sum > threshold;     % Maska aktywności

% Znajdź segmenty (grupy aktywnych okien czasowych)
changes = diff([0; is_active; 0]);
rising_edges = find(changes == 1);
falling_edges = find(changes == -1) - 1;
segments = [rising_edges, falling_edges];

% Filtracja krótkich segmentów
valid_segments = [];
for i = 1:size(segments, 1)
    segment_duration = (segments(i, 2) - segments(i, 1) + 1) * window_step/fs;
    if segment_duration >= min_duration
        valid_segments = [valid_segments; segments(i,:)];
    end
end
segments = valid_segments;

% Dekodowanie cyfr
detected_digits = [];
detected_times = [];

for i = 1:size(segments, 1)
    start_w = segments(i, 1);
    end_w = segments(i, 2);
    
    % Znajdź okno z najlepszym stosunkiem sygnału w segmencie
    segment_energies = freq_energies(start_w:end_w, :);
    [~, best_window] = max(sum(segment_energies.^2, 2));
    best_window_idx = start_w + best_window - 1;
    
    % Analizuj energię częstotliwości w najlepszym oknie
    row_energies = freq_energies(best_window_idx, 1:length(row_freq));
    col_energies = freq_energies(best_window_idx, length(row_freq)+1:end);
    
    % Znajdź najsilniejsze częstotliwości
    [max_row_val, row_idx] = max(row_energies);
    [max_col_val, col_idx] = max(col_energies);
    
    % Oblicz kontrast (stosunek najwyższej wartości do drugiej najwyższej)
    row_sorted = sort(row_energies, 'descend');
    col_sorted = sort(col_energies, 'descend');
    
    % Unikaj dzielenia przez zero
    row_contrast = row_sorted(1) / (row_sorted(2) + 0.0001);
    col_contrast = col_sorted(1) / (col_sorted(2) + 0.0001);
    
    % Sprawdź czy wyniki są wiarygodne
    if max_row_val > 0.4 && max_col_val > 0.4 && row_contrast > min_contrast && col_contrast > min_contrast
        detected_digit = dtmf_keys(row_idx, col_idx);
        detected_digits = [detected_digits, detected_digit];
        detected_times = [detected_times, window_times(best_window_idx)];
        
        fprintf('Wykryto cyfrę: %s (czas: %.3fs)\n', detected_digit, window_times(best_window_idx));
    end
end

% Usuń duplikaty cyfr (te same cyfry wykryte blisko siebie)
if length(detected_digits) > 1
    unique_digits = detected_digits(1);
    unique_times = detected_times(1);
    min_time_diff = 0.3; % 300ms
    
    for i = 2:length(detected_digits)
        if detected_times(i) - unique_times(end) > min_time_diff || detected_digits(i) ~= unique_digits(end)
            unique_digits = [unique_digits, detected_digits(i)];
            unique_times = [unique_times, detected_times(i)];
        end
    end
    
    detected_digits = unique_digits;
    detected_times = unique_times;
end

% Wyświetl wyniki
if ~isempty(detected_digits)
    disp(['Zdekodowana sekwencja: ' detected_digits]);
else
    disp('Nie wykryto żadnych cyfr DTMF.');
end

% Wizualizacja wyników
figure('Name', 'Uproszczona segmentacja i wykryte cyfry DTMF');

% Panel 1: Sygnał i segmentacja
subplot(2,1,1);
plot((0:length(s_filtered)-1)/fs, s_filtered, 'Color', [0.7 0.7 0.7]);
hold on;

% Zaznacz segmenty
for i = 1:size(segments, 1)
    start_time = window_times(segments(i, 1));
    end_time = window_times(segments(i, 2));
    y_max = 1.1 * max(abs(s_filtered));
    
    % Zaznacz segment prostokątem
    x_rect = [start_time, end_time, end_time, start_time];
    y_rect = [-y_max, -y_max, y_max, y_max];
    patch(x_rect, y_rect, [0.9, 0.9, 0.3], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

% Dodaj wykryte cyfry
for i = 1:length(detected_times)
    x_pos = detected_times(i);
    y_val = 0.9 * max(abs(s_filtered));
    plot([x_pos x_pos], [-y_val y_val], 'r--', 'LineWidth', 1);
    text(x_pos+0.05, y_val, detected_digits(i), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
end

title('Wykryte segmenty i cyfry DTMF');
xlabel('Czas [s]');
ylabel('Amplituda');
hold off;

% Panel 2: Energie częstotliwości w czasie
subplot(2,1,2);
imagesc(window_times, 1:length(all_freq), freq_energies');
colorbar;
hold on;

% Dodaj oznaczenia wykrytych cyfr
for i = 1:length(detected_times)
    plot([detected_times(i) detected_times(i)], [0.5, length(all_freq)+0.5], 'r--', 'LineWidth', 1);
    text(detected_times(i), length(all_freq)+0.7, detected_digits(i), 'Color', 'r', 'FontWeight', 'bold');
end

title('Energie częstotliwości DTMF w czasie');
xlabel('Czas [s]');
ylabel('Indeks częstotliwości');
set(gca, 'YTick', 1:length(all_freq), 'YTickLabel', num2str(all_freq'));
hold off;
%% 6. WALIDACJA WYNIKÓW I TESTOWANIE SEKWENCJI
disp('=== WALIDACJA I TESTOWANIE ===');

% Walidacja z sekwencją wzorcową dla pliku s.wav
if contains(filename, 's.wav')
    expected_sequence = '123456789*0#';
    disp(['Oczekiwana sekwencja: ' expected_sequence]);
    
    % Porównanie z wykrytą sekwencją
    if strcmp(detected_digits, expected_sequence)
        disp('SUKCES: Wykryta sekwencja zgadza się z oczekiwaną!');
        success_rate = 100;
    else
        % Oblicz dokładność rozpoznawania
        len_detected = length(detected_digits);
        len_expected = length(expected_sequence);
        
        % Oblicz odległość Levenshteina (metryka różnicy między sekwencjami)
        lev_dist = zeros(len_detected+1, len_expected+1);
        
        for i = 0:len_detected
            lev_dist(i+1, 1) = i;
        end
        
        for j = 0:len_expected
            lev_dist(1, j+1) = j;
        end
        
        for i = 1:len_detected
            for j = 1:len_expected
                if detected_digits(i) == expected_sequence(j)
                    cost = 0;
                else
                    cost = 1;
                end
                lev_dist(i+1, j+1) = min([lev_dist(i, j+1) + 1, ...
                                         lev_dist(i+1, j) + 1, ...
                                         lev_dist(i, j) + cost]);
            end
        end
        
        % Oblicz dokładność na podstawie odległości edycji
        error_rate = lev_dist(end, end) / max(len_detected, len_expected);
        success_rate = 100 * (1 - error_rate);
        
        disp(['UWAGA: Wykryta sekwencja nie zgadza się z oczekiwaną.']);
        disp(['Odległość edycji: ' num2str(lev_dist(end, end))]);
        disp(['Dokładność rozpoznawania: ' num2str(success_rate, '%.1f') '%']);
    end
end