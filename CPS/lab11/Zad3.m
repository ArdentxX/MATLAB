%% ZADANIE 3: Podpasmowe kodowanie dźwięku (1 pkt)

function [y, compression_ratio] = subband_encode_decode(x, num_subbands, bits_per_band)
    % Uproszczone podpasmowe kodowanie
    
    % Podział na podpasma (uproszczony - używamy filtrów prostokątnych)
    N = length(x);
    X = fft(x);
    
    % Podział widma na podpasma
    band_size = floor(length(X) / num_subbands);
    subband_coeffs = cell(num_subbands, 1);
    
    for band = 1:num_subbands
        start_idx = (band-1) * band_size + 1;
        if band == num_subbands
            end_idx = length(X);
        else
            end_idx = band * band_size;
        end
        
        subband_coeffs{band} = X(start_idx:end_idx);
    end
    
    % Kwantyzacja każdego podpasma
    quantized_coeffs = cell(num_subbands, 1);
    original_bits = 16 * N; % 16 bitów na próbkę
    used_bits = 0;
    
    for band = 1:num_subbands
        if length(bits_per_band) == 1
            bits = bits_per_band;
        else
            bits = bits_per_band(min(band, length(bits_per_band)));
        end
        
        coeffs = subband_coeffs{band};
        
        % Kwantyzacja
        if bits > 0
            max_val = max(abs(coeffs));
            if max_val > 0
                levels = 2^bits;
                step = 2 * max_val / levels;
                quantized_coeffs{band} = round(coeffs / step) * step;
            else
                quantized_coeffs{band} = coeffs;
            end
            used_bits = used_bits + bits * length(coeffs);
        else
            quantized_coeffs{band} = zeros(size(coeffs));
        end
    end
    
    % Rekonstrukcja
    Y = zeros(size(X));
    for band = 1:num_subbands
        start_idx = (band-1) * band_size + 1;
        if band == num_subbands
            end_idx = length(Y);
        else
            end_idx = band * band_size;
        end
        
        Y(start_idx:end_idx) = quantized_coeffs{band};
    end
    
    y = real(ifft(Y));
    compression_ratio = original_bits / used_bits;
end

% Testowanie różnych wariantów
x_test = x(1:44100); % 1 sekunda

% Wariant 1: 8 podpasm, 6 bitów każde
[y1, cr1] = subband_encode_decode(x_test, 8, 6);

% Wariant 2: 32 podpasma, 6 bitów każde  
[y2, cr2] = subband_encode_decode(x_test, 32, 6);

% Wariant 3: 32 podpasma, zmienna liczba bitów
bits_variable = [8, 8, 7, 6, 4, 4, 4, 4]; % powtarzaj wzorzec
bits_full = repmat(bits_variable, 1, ceil(32/length(bits_variable)));
bits_full = bits_full(1:32);
[y3, cr3] = subband_encode_decode(x_test, 32, bits_full);

fprintf('Współczynniki kompresji:\n');
fprintf('8 podpasm, 6 bitów: %.2f\n', cr1);
fprintf('32 podpasma, 6 bitów: %.2f\n', cr2);  
fprintf('32 podpasma, zmienne bity: %.2f\n', cr3);

% Spektrogramy
figure(2);
subplot(2,2,1); spectrogram(x_test, 256, 128, 256, 44100, 'yaxis'); title('Oryginalny');
subplot(2,2,2); spectrogram(y1, 256, 128, 256, 44100, 'yaxis'); title('8 podpasm');
subplot(2,2,3); spectrogram(y2, 256, 128, 256, 44100, 'yaxis'); title('32 podpasma'); 
subplot(2,2,4); spectrogram(y3, 256, 128, 256, 44100, 'yaxis'); title('32 podpasma zmienne');

% OPCJONALNIE: Adaptacyjna alokacja bitów na podstawie energii
function [y_adaptive, bit_map] = adaptive_subband_coding(x, num_subbands, total_bits_per_frame)
    frame_size = 1024;
    num_frames = floor(length(x) / frame_size);
    y_adaptive = zeros(size(x));
    bit_map = zeros(num_subbands, num_frames);
    
    for frame = 1:num_frames
        start_idx = (frame-1) * frame_size + 1;
        end_idx = frame * frame_size;
        frame_data = x(start_idx:end_idx);
        
        % Analiza energii w podpasmach
        X = fft(frame_data);
        band_size = floor(length(X) / num_subbands);
        energies = zeros(num_subbands, 1);
        
        for band = 1:num_subbands
            band_start = (band-1) * band_size + 1;
            if band == num_subbands
                band_end = length(X);
            else
                band_end = band * band_size;
            end
            energies(band) = sum(abs(X(band_start:band_end)).^2);
        end
        
        % Alokacja bitów proporcjonalnie do energii
        total_energy = sum(energies);
        if total_energy > 0
            bit_map(:, frame) = round(total_bits_per_frame * energies / total_energy);
            bit_map(bit_map(:, frame) < 1, frame) = 1; % minimum 1 bit
        else
            bit_map(:, frame) = ones(num_subbands, 1);
        end
        
        % Kodowanie z adaptacyjnymi bitami
        [frame_reconstructed, ~] = subband_encode_decode(frame_data, num_subbands, bit_map(:, frame));
        y_adaptive(start_idx:end_idx) = frame_reconstructed;
    end
end