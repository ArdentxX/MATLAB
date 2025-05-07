# Filtry specjalne i modulacje w cyfrowym przetwarzaniu sygnałów

## Spis treści

1. [Wprowadzenie](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#wprowadzenie)
2. [Transformata Hilberta i jej zastosowania](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#transformata-hilberta-i-jej-zastosowania)
    - [Teoria transformaty Hilberta](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#teoria-transformaty-hilberta)
    - [Implementacja filtru Hilberta](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#implementacja-filtru-hilberta)
3. [Modulacja amplitudy (AM)](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#modulacja-amplitudy-am)
    - [DSB-C: Double Side Band with Carrier](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#dsb-c-double-side-band-with-carrier)
    - [DSB-SC: Double Side Band with Suppressed Carrier](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#dsb-sc-double-side-band-with-suppressed-carrier)
    - [SSB-SC: Single Side Band with Suppressed Carrier](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#ssb-sc-single-side-band-with-suppressed-carrier)
4. [Demodulacja sygnałów AM](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#demodulacja-sygna%C5%82%C3%B3w-am)
    - [Demodulacja obwiedniowa](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#demodulacja-obwiedniowa)
    - [Metody demodulacji dla różnych typów modulacji](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#metody-demodulacji-dla-r%C3%B3%C5%BCnych-typ%C3%B3w-modulacji)
5. [Filtry interpolacyjne i decymacyjne](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#filtry-interpolacyjne-i-decymacyjne)
    - [Nadpróbkowanie i interpolacja](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#nadpr%C3%B3bkowanie-i-interpolacja)
    - [Decymacja i filtry antyaliasingowe](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#decymacja-i-filtry-antyaliasingowe)
6. [Modulacja częstotliwości (FM)](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#modulacja-cz%C4%99stotliwo%C5%9Bci-fm)
    - [Teoria modulacji FM](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#teoria-modulacji-fm)
    - [Metody demodulacji FM](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#metody-demodulacji-fm)
7. [Szczegółowa analiza zadań laboratoryjnych](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#szczeg%C3%B3%C5%82owa-analiza-zada%C5%84-laboratoryjnych)
    - [Zadanie 3: Filtry interpolatora i decymatora, mikser audio](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#zadanie-3-filtry-interpolatora-i-decymatora-mikser-audio)
    - [Zadanie 4: Filtr różniczkujący, dekodowanie FM](https://claude.ai/chat/ff20823e-2d7a-4142-8a18-6dbe61929698#zadanie-4-filtr-r%C3%B3%C5%BCniczkuj%C4%85cy-dekodowanie-fm)

## Wprowadzenie

Cyfrowe przetwarzanie sygnałów (CPS) jest dziedziną inżynierii, która zajmuje się analizą i modyfikacją sygnałów cyfrowych. W ramach CPS szczególne znaczenie mają filtry cyfrowe, które pozwalają na wyodrębnienie pożądanych cech sygnału lub usunięcie zakłóceń. Filtry specjalne, takie jak filtr Hilberta, filtry różniczkujące, interpolacyjne i decymacyjne, mają unikalne właściwości przydatne w konkretnych zastosowaniach, takich jak modulacje sygnałów, obróbka dźwięku czy transmisja danych.

Niniejszy dokument koncentruje się na filtrach specjalnych i ich zastosowaniach w kontekście modulacji i demodulacji sygnałów analogowych, ze szczególnym uwzględnieniem modulacji amplitudy (AM) i częstotliwości (FM).

## Transformata Hilberta i jej zastosowania

### Teoria transformaty Hilberta

Transformata Hilberta to operacja, która przekształca funkcję rzeczywistą w funkcję zespoloną. Dla sygnału rzeczywistego x(t), jego transformata Hilberta ℋ[x(t)] jest zdefiniowana jako:

$$\mathcal{H}[x(t)] = \frac{1}{\pi} \int_{-\infty}^{\infty} \frac{x(\tau)}{t - \tau} d\tau$$

Jest to splot sygnału x(t) z funkcją 1/(πt).

W dziedzinie częstotliwości, transformata Hilberta odpowiada przesunięciu fazowemu o -π/2 dla częstotliwości dodatnich i +π/2 dla częstotliwości ujemnych, bez zmiany amplitudy sygnału:

$$\mathcal{H}[x(t)] \Leftrightarrow -j \cdot \text{sgn}(f) \cdot X(f)$$

gdzie X(f) jest transformatą Fouriera sygnału x(t), a sgn(f) jest funkcją znaku częstotliwości:

$$\text{sgn}(f) = \begin{cases} 1, & \text{dla } f > 0 \ 0, & \text{dla } f = 0 \ -1, & \text{dla } f < 0 \end{cases}$$

### Implementacja filtru Hilberta

W praktyce, idealna transformata Hilberta nie jest realizowalna ze względu na nieskończoną odpowiedź impulsową. Dlatego stosuje się przybliżenie w postaci filtru FIR (o skończonej odpowiedzi impulsowej). Odpowiedź impulsowa filtru Hilberta dla dyskretnych sygnałów jest określona jako:

$$h[n] = \begin{cases} \frac{2}{\pi n} \sin^2(\frac{\pi n}{2}), & \text{dla } n \neq 0 \ 0, & \text{dla } n = 0 \end{cases}$$

Ponieważ odpowiedź impulsowa jest nieskończona, musi być przycięta i pomnożona przez okno (np. okno Hamminga) w celu uzyskania praktycznego filtru:

$$h_{\text{practical}}[n] = h[n] \cdot w[n]$$

gdzie w[n] jest funkcją okna.

Filtr Hilberta ma kilka istotnych zastosowań:

- Wyznaczanie obwiedni sygnału (używane w demodulacji AM)
- Generowanie sygnału analitycznego
- Modulacje SSB (Single Side Band)
- Wykrywanie fazowe w komunikacji cyfrowej

## Modulacja amplitudy (AM)

Modulacja amplitudy (AM) jest jedną z najstarszych i najprostszych technik modulacji, gdzie informacja jest kodowana poprzez zmiany amplitudy sygnału nośnego o wysokiej częstotliwości. W modulacji AM, sygnał informacyjny (modulujący) wpływa na amplitudę sygnału nośnego, zachowując jego częstotliwość.

### DSB-C: Double Side Band with Carrier

DSB-C to standardowa modulacja AM, w której sygnał nośny jest obecny w sygnale modulowanym, a widmo sygnału zawiera dwie wstęgi boczne (po obu stronach częstotliwości nośnej). Matematycznie, sygnał DSB-C można wyrazić jako:

$$y_{DSB-C}(t) = (1 + dA \cdot x(t)) \cdot \cos(2\pi f_c t)$$

gdzie:

- x(t) to sygnał modulujący
- f_c to częstotliwość nośna
- dA to współczynnik głębokości modulacji (0 < dA ≤ 1)

Współczynnik dA określa stopień modulacji i wpływa na efektywność wykorzystania mocy. Gdy dA > 1, występuje przemodelowanie, co prowadzi do zniekształceń sygnału.

### DSB-SC: Double Side Band with Suppressed Carrier

W modulacji DSB-SC, sygnał nośny jest tłumiony (nie jest obecny w sygnale zmodulowanym), co zwiększa efektywność energetyczną transmisji. Formula sygnału DSB-SC:

$$y_{DSB-SC}(t) = x(t) \cdot \cos(2\pi f_c t)$$

DSB-SC jest bardziej wydajny energetycznie niż DSB-C, ponieważ energia nie jest tracona na transmisję fali nośnej, która nie zawiera informacji. Jednakże demodulacja DSB-SC jest bardziej złożona i wymaga synchronizacji fazy.

### SSB-SC: Single Side Band with Suppressed Carrier

SSB-SC to zaawansowana forma modulacji AM, gdzie nośna jest tłumiona i transmitowana jest tylko jedna wstęga boczna (górna - USB lub dolna - LSB). Dzięki temu szerokość pasma sygnału jest zredukowana o połowę w porównaniu do DSB-SC. Matematyczna reprezentacja SSB-SC wykorzystuje transformatę Hilberta:

$$y_{SSB-SC}(t) = 0.5 \cdot x(t) \cdot \cos(2\pi f_c t) \pm 0.5 \cdot \mathcal{H}[x(t)] \cdot \sin(2\pi f_c t)$$

Znak "+" daje dolną wstęgę boczną (LSB), a znak "-" górną wstęgę boczną (USB).

SSB-SC oferuje największą efektywność spektralną spośród modulacji AM, ale wymaga dokładnej filtracji lub implementacji transformaty Hilberta.

## Demodulacja sygnałów AM

### Demodulacja obwiedniowa

Dla sygnałów DSB-C, najprostszą metodą demodulacji jest detekcja obwiedni. Proces ten wymaga wyznaczenia obwiedni sygnału zmodulowanego, która odpowiada sygnałowi modulującemu plus stała składowa.

Jedną z metod wyznaczania obwiedni jest wykorzystanie transformaty Hilberta. Dla sygnału rzeczywistego x(t), obwiednia jest obliczana jako:

$$\text{envelope}(t) = \sqrt{x^2(t) + \mathcal{H}[x(t)]^2}$$

W kontekście sygnału DSB-C, obwiednia zawiera oryginalny sygnał modulujący plus składową stałą, którą można usunąć.

### Metody demodulacji dla różnych typów modulacji

1. **Demodulacja DSB-C**:
    
    - Detekcja obwiedni (przy użyciu transformaty Hilberta)
    - Odjęcie składowej stałej
    - Skalowanie według głębokości modulacji
2. **Demodulacja DSB-SC**:
    
    - Mnożenie przez lokalnie generowaną falę nośną o tej samej częstotliwości i fazie
    - Filtracja dolnoprzepustowa wynikowego sygnału
    - Wymaga dokładnej synchronizacji fazy z nośną
3. **Demodulacja SSB-SC**:
    
    - Metoda I/Q (odbieranie sygnału jako zespolonego)
    - Wykorzystanie właściwości transformaty Hilberta do odzyskania oryginalnego sygnału
    - Dla USB: x(t) = demod - ℋ[demod]
    - Dla LSB: x(t) = demod + ℋ[demod]

gdzie demod to sygnał po przemnożeniu przez lokalną nośną i filtracji dolnoprzepustowej.

## Filtry interpolacyjne i decymacyjne

### Nadpróbkowanie i interpolacja

Nadpróbkowanie (upsampling) to proces zwiększania częstotliwości próbkowania sygnału cyfrowego. Podstawowa operacja nadpróbkowania polega na wstawieniu L-1 zer między każdą parę próbek oryginalnego sygnału, gdzie L jest współczynnikiem nadpróbkowania.

Proces interpolacji wykorzystuje filtry dolnoprzepustowe do wygładzenia sygnału po wstawieniu zer. Filtr interpolacyjny zapełnia wstawione zera wartościami, które przybliżają oryginalny sygnał ciągły. Odpowiedź impulsowa filtru interpolacyjnego powinna mieć kształt funkcji sinc, która jest idealnym filtrem dolnoprzepustowym:

$$h_{ideal}[n] = \text{sinc}(n/L)$$

W praktyce stosuje się przybliżenia filtru idealnego, takie jak filtry FIR zaprojektowane metodami okien lub optymalizacji.

### Decymacja i filtry antyaliasingowe

Decymacja (downsampling) to proces zmniejszania częstotliwości próbkowania poprzez usunięcie próbek. W najprostszej formie, decymacja o współczynniku M polega na zachowaniu tylko co M-tej próbki.

Przed wykonaniem decymacji konieczne jest zastosowanie filtru antyaliasingowego (dolnoprzepustowego), aby zapobiec nakładaniu widma (aliasing). Filtr ten powinien mieć częstotliwość odcięcia nie większą niż π/M, gdzie π odpowiada połowie częstotliwości próbkowania.

Efektywne implementacje często łączą filtrację i decymację, redukując liczbę operacji. Przykładem są filtry polifazowe, które dzielą filtr na M podfiltrów, z których każdy przetwarza tylko te próbki, które zostaną zachowane po decymacji.

## Modulacja częstotliwości (FM)

### Teoria modulacji FM

W modulacji częstotliwości (FM), informacja jest kodowana poprzez zmiany częstotliwości sygnału nośnego, proporcjonalne do sygnału modulującego. Matematyczna reprezentacja sygnału FM:

$$x_{FM}(t) = \cos\left(2\pi f_c t + 2\pi K \int_{0}^{t} x(\tau) d\tau\right) = \cos(2\pi f_c t + \phi(t))$$

gdzie:

- f_c to częstotliwość nośna
- K to współczynnik dewiacji częstotliwości
- x(t) to sygnał modulujący
- φ(t) to faza chwilowa

Dewiacja częstotliwości Δf = K·max(|x(t)|) określa maksymalne odchylenie częstotliwości chwilowej od częstotliwości nośnej.

### Metody demodulacji FM

Istnieje kilka metod demodulacji sygnałów FM:

1. **Detekcja obwiedni po różniczkowaniu**:
    
    - Różniczkowanie sygnału FM: $$\frac{dx_{FM}(t)}{dt} = -(2\pi f_c + 2\pi K \cdot x(t)) \sin(2\pi f_c t + \phi(t))$$
    - Wyznaczenie obwiedni zróżniczkowanego sygnału
    - Odjęcie składowej stałej i skalowanie
    
    Proces można zrealizować za pomocą:
    
    - Filtru różniczkującego o charakterystyce pasmowo-przepustowej
    - Kaskady filtru pasmowo-przepustowego i różniczkującego
    - Detekcji kwadraturowej z obliczeniem obwiedni
2. **Demodulacja I/Q (kwadraturowa)**:
    
    - Wyznaczenie składowych I (in-phase) i Q (quadrature): $$I(n) = \text{LowPass}(x_{FM}(t) \cdot \cos(2\pi f_c t))$$ $$Q(n) = \text{LowPass}(-x_{FM}(t) \cdot \sin(2\pi f_c t))$$
    - Utworzenie sygnału analitycznego: y[n] = I[n] + j·Q[n]
    - Obliczenie różnicy faz między próbkami: x[n] ∝ arg(y[n]·conj(y[n-1]))
3. **Pętla fazowa (PLL)**:
    
    - Wykorzystanie pętli fazowej do śledzenia zmian częstotliwości
    - Sygnał sterujący VCO (voltage-controlled oscillator) zawiera zdemodulowany sygnał

## Szczegółowa analiza zadań laboratoryjnych

### Zadanie 3: Filtry interpolatora i decymatora, mikser audio

Celem zadania 3 jest zrozumienie i implementacja procesu zmiany częstotliwości próbkowania (resampling) sygnałów audio oraz ich miksowanie. Zadanie bazuje na koncepcjach nadpróbkowania, interpolacji i decymacji.

#### Teoretyczne podstawy zadania

Zmiana częstotliwości próbkowania jest kluczowym procesem w systemach cyfrowego przetwarzania dźwięku, szczególnie gdy miksuje się sygnały pochodzące z różnych źródeł o różnych częstotliwościach próbkowania. Główne kroki procesu repróbkowania obejmują:

1. **Nadpróbkowanie (upsampling)** - zwiększenie częstotliwości próbkowania przez wstawienie zer między próbki oryginalne
2. **Interpolacja** - filtracja dolnoprzepustowa dla wygładzenia sygnału
3. **Decymacja (downsampling)** - zmniejszenie częstotliwości próbkowania przez usunięcie próbek (w przypadku gdy współczynnik konwersji nie jest liczbą całkowitą)

#### Analiza implementacji w MATLAB

W zadaniu 3 wykorzystano dwie metody repróbkowania:

##### Metoda 1: Funkcja resample

Pierwsza metoda wykorzystuje gotową funkcję `resample` dostępną w MATLAB, która automatycznie wykonuje wszystkie kroki potrzebne do zmiany częstotliwości próbkowania:

```matlab
[p1, q1] = rat(fs_target/fs1); % Obliczenie współczynników interpolacji/decymacji
[p2, q2] = rat(fs_target/fs2);

x1_resampled = resample(x1, p1, q1); % Repróbkowanie sygnałów
x2_resampled = resample(x2, p2, q2);
```

Funkcja `rat` oblicza przybliżenie ułamkowe współczynnika zmiany próbkowania jako p/q, gdzie:

- p to współczynnik nadpróbkowania (upsampling)
- q to współczynnik decymacji (downsampling)

Dla przykładu, dla konwersji z 8000 Hz do 48000 Hz, współczynnik p=6, q=1, co oznacza sześciokrotne nadpróbkowanie.

##### Metoda 2: Ręczna implementacja

Druga metoda to ręczna implementacja procesu nadpróbkowania i interpolacji, co umożliwia lepsze zrozumienie zachodzących procesów:

1. **Nadpróbkowanie** - wstawienie zer między próbki:

```matlab
up_factor1 = fs_target / fs1;
x1_upsampled = zeros(length(x1) * up_factor1, 1);
x1_upsampled(1:up_factor1:end) = x1;
```

2. **Projektowanie filtru interpolacyjnego**:

```matlab
cutoff1 = (fs1/2) / fs_target;
filter_order1 = 48;
b_interp1 = fir1(filter_order1, cutoff1, 'low') * up_factor1;
```

Kluczowe parametry filtru:

- Częstotliwość odcięcia (`cutoff1`) jest ustawiona na połowę oryginalnej częstotliwości próbkowania znormalizowaną do nowej częstotliwości próbkowania
- Mnożenie przez `up_factor1` zapewnia zachowanie energii sygnału po interpolacji

3. **Filtracja interpolacyjna i kompensacja opóźnienia**:

```matlab
x1_manual = filter(b_interp1, 1, x1_upsampled);
delay1 = filter_order1 / 2;
x1_manual = x1_manual(delay1+1:end);
x1_manual = [x1_manual; zeros(delay1, 1)];
```

Dla bardziej złożonego przypadku konwersji z 32000 Hz do 48000 Hz, proces wymaga zarówno nadpróbkowania jak i decymacji:

```matlab
up_factor2 = 3;
down_factor2 = 2;
```

Co odpowiada konwersji (32000 × 3) ÷ 2 = 48000 Hz.

#### Porównanie metod i wyniki

Kod oblicza błąd średniokwadratowy (MSE) między sygnałami uzyskanymi obiema metodami a teoretycznym sygnałem wzorcowym:

```matlab
mse_resampled = mean((x4_resampled - x4_expected_1).^2);
mse_manual = mean((x4_manual - x4_expected_2).^2);
```

Zgodnie z wynikami prezentowanymi w kodzie, metoda wykorzystująca funkcję `resample` zwykle daje mniejszy błąd, co wynika z faktu, że funkcja ta implementuje zaawansowane techniki takie jak filtry polifazowe, które są bardziej efektywne obliczeniowo i dokładniejsze.

#### Analiza widmowa

Program wykonuje również analizę widmową najlepszej metody (funkcji `resample`), co pozwala na obserwację składowych częstotliwościowych sygnału wynikowego:

```matlab
NFFT = 2^nextpow2(length(x4_resampled));
X4 = fft(x4_resampled, NFFT) / length(x4_resampled);
```

Widmo pokazuje trzy wyraźne piki odpowiadające częstotliwościom sygnałów składowych (f1 = 1001.2 Hz, f2 = 303.1 Hz i f3 = 2110.4 Hz), co potwierdza poprawność procesu repróbkowania.

#### Praktyczne zastosowanie - miksowanie plików WAV

Zadanie kończy się praktycznym przykładem zastosowania repróbkowania do miksowania rzeczywistych plików dźwiękowych w formacie WAV:

```matlab
[wav1, fs_wav1] = audioread('x1.wav');
[wav2, fs_wav2] = audioread('x2.wav');

[p_wav1, q_wav1] = rat(fs_target/fs_wav1); 
[p_wav2, q_wav2] = rat(fs_target/fs_wav2);

wav1_resampled = resample(wav1, p_wav1, q_wav1);
wav2_resampled = resample(wav2, p_wav2, q_wav2);

wav_mixed = wav1_resampled + wav2_resampled;
```

Ten fragment kodu pokazuje typowy proces w produkcji audio, gdzie sygnały z różnych źródeł muszą być najpierw przekonwertowane do wspólnej częstotliwości próbkowania przed zmiksowaniem.

### Zadanie 4: Filtr różniczkujący, dekodowanie FM

Zadanie 4 dotyczy implementacji i porównania różnych metod demodulacji sygnału FM. Jest to praktyczne zastosowanie teorii modulacji częstotliwości i filtrów różniczkujących.

#### Teoretyczne podstawy demodulacji FM

Demodulacja FM polega na odzyskaniu informacji zawartej w zmianach częstotliwości sygnału nośnego. Sygnał FM można wyrazić jako:

$$x_{FM}(t) = \cos\left(2\pi f_c t + 2\pi K \int_{0}^{t} x(\tau) d\tau\right)$$

gdzie sygnał modulujący x(t) wpływa na fazę chwilową sygnału nośnego. Pochodna fazy chwilowej jest proporcjonalna do sygnału modulującego:

$$\frac{d\phi(t)}{dt} = 2\pi K \cdot x(t)$$

#### Analiza implementacji w MATLAB

W zadaniu 4 zaimplementowano trzy różne metody demodulacji sygnału FM:

##### Metoda 1: Demodulacja przez obliczenie składowych I(n) oraz Q(n)

Ta metoda wykorzystuje demodulację kwadraturową, która polega na przemnożeniu sygnału FM przez sygnały nośne przesunięte w fazie o 90 stopni:

```matlab
I = lowpass(x .* cos(2*pi*fc*t'), fc+5e3, fs); % Składowa I (w fazie)
Q = lowpass(-x .* sin(2*pi*fc*t'), fc+5e3, fs); % Składowa Q (w kwadraturze)

y = I + 1j*Q; % Utworzenie sygnału analitycznego

y_shifted = [y(2:end); 0]; % Przesunięcie o jedną próbkę
phase_diff = angle(y .* conj(y_shifted)); % Różnica faz między próbkami
```

Kluczowe kroki:

1. Przemnożenie sygnału przez funkcje sinus i cosinus o częstotliwości nośnej
2. Filtracja dolnoprzepustowa dla wyodrębnienia składowych baseband
3. Utworzenie sygnału analitycznego (zespolonego)
4. Obliczenie różnicy faz między kolejnymi próbkami, która jest proporcjonalna do sygnału modulującego

##### Metoda 2: Demodulacja przez różniczkowanie i obwiednię

Ta metoda bazuje na teorii, że obwiednia zróżniczkowanego sygnału FM zawiera informację o sygnale modulującym:

```matlab
% Projektowanie filtru różniczkującego
h_diff_ideal = zeros(size(n_diff));
for i = 1:length(n_diff)
    if n_diff(i) == 0
        h_diff_ideal(i) = 0;
    else
        h_diff_ideal(i) = (sin(pi*n_diff(i)) - pi*n_diff(i)*cos(pi*n_diff(i))) / (pi*n_diff(i)^2);
    end
end

% Zastosowanie okna Blackmana
window_diff = blackman(Ndiff)';
h_diff = h_diff_ideal .* window_diff;

% Projektowanie filtru pasmowo-przepustowego
h_bp = fir1(Nbp-1, [f_low, f_high], 'bandpass', blackman(Nbp));

% Kaskada filtrów
h_diff_bp = conv(h_diff, h_bp);

% Filtracja sygnału FM
x_filtered = filter(h_diff_bp, 1, x);

% Odzyskanie obwiedni
x_squared = x_filtered.^2;
x_env = filter(h_lp, 1, x_squared);
x_demod_method2 = sqrt(abs(x_env));
```

Kluczowe kroki:

1. Projektowanie filtru różniczkującego jako filtru FIR z oknem Blackmana
2. Projektowanie filtru pasmowo-przepustowego dla wyodrębnienia sygnału wokół częstotliwości nośnej
3. Połączenie obu filtrów (splot odpowiedzi impulsowych)
4. Filtracja sygnału FM i obliczenie obwiedni przez kwadraturowanie i filtrację dolnoprzepustową

##### Metoda 3: Demodulacja przez kaskadę oddzielnych filtrów BP i DIFF

Ta metoda jest modyfikacją metody 2, gdzie zamiast jednego złożonego filtru stosuje się oddzielne filtry pasmowo-przepustowy i różniczkujący:

```matlab
% Filtracja pasmowo-przepustowa
x_bp_filtered = filter(h_bp, 1, x);

% Prosty filtr różniczkujący
b_diff_simple = [-1, 1];
a_diff_simple = 1;

% Różniczkowanie
x_diff_filtered = filter(b_diff_simple, a_diff_simple, x_bp_filtered);

% Kwadraturowanie i obwiednia
x_squared_method3 = x_diff_filtered.^2;
x_env_method3 = filter(h_lp, 1, x_squared_method3);
x_demod_method3 = sqrt(abs(x_env_method3));
```

Kluczowe różnice:

1. Najpierw stosuje się filtr pasmowo-przepustowy
2. Następnie stosuje się prosty filtr różniczkujący (y[n] = x[n] - x[n-1])
3. Dalszy proces jest identyczny jak w metodzie 2

#### Metoda bonusowa: Implementacja filtru różniczkującego pasmowo-przepustowego za pomocą firls

Dodatkowo kod zawiera implementację filtru różniczkującego pasmowo-przepustowego z wykorzystaniem funkcji `firls`, która projektuje filtr metodą najmniejszych kwadratów:

```matlab
f = [0, f_low-0.01, f_low, f_high, f_high+0.01, 1]; % Punkty częstotliwości
a = [0, 0, 1j, -1j, 0, 0]; % Odpowiedź różniczkująca

h_bp_diff = firls(Nbp_diff-1, f, abs(a));
h_bp_diff = real(h_bp_diff);
```

Ten podejście umożliwia zaprojektowanie filtru o złożonej charakterystyce bezpośrednio, bez konieczności kaskadowego łączenia dwóch oddzielnych filtrów.

#### Porównanie metod i wyniki

Kod generuje wykresy porównujące zdemodulowane sygnały uzyskane różnymi metodami:

```matlab
subplot(3,1,1);
plot
```