%% ZADANIE 4: Szacowanie ilości informacji - Entropia (+0.1 pkt)

function H = calculate_entropy(x)
    % Obliczenie entropii sygnału
    
    % Znajdź unikalne symbole
    unique_symbols = unique(x);
    N = length(x);
    
    % Oblicz prawdopodobieństwa
    probabilities = zeros(size(unique_symbols));
    for i = 1:length(unique_symbols)
        probabilities(i) = sum(x == unique_symbols(i)) / N;
    end
    
    % Oblicz entropię
    H = -sum(probabilities .* log2(probabilities));
    
    fprintf('Symbole: '); disp(unique_symbols');
    fprintf('Prawdopodobieństwa: '); disp(probabilities');
    fprintf('Entropia H = %.4f bitów/symbol\n', H);
end

% Testowe sygnały
x1 = [0, 1, 2, 3, 3, 2, 1, 0];
x2 = [0, 7, 0, 2, 0, 2, 0, 7, 4, 2];
x3 = [0, 0, 0, 0, 0, 0, 0, 15];

fprintf('Entropia dla x1:\n');
H1 = calculate_entropy(x1);

fprintf('\nEntropia dla x2:\n');
H2 = calculate_entropy(x2);

fprintf('\nEntropia dla x3:\n');
H3 = calculate_entropy(x3);
