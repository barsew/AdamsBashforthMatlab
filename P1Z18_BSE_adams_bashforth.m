function [y, h, x] = P1Z18_BSE_adams_bashforth(alfa, beta, n, y_alfa, f, c, b, a)
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Rozwiązywanie liniowych równań różniczkowych pierwszego i drugiego rzędu
% metodą Adamsa-Bashfortha rzedu 4-go
% (pierwsze trzy kroki metodą Rungego-Kutty rzędu 4-go (wzór Gilla))
% równania postaci a(x)*y'' + b(x)*y' + c(x)*y = f(x)
% Wejście:
%   [alfa, beta] - liczby rzeczywiste, tworzą przedział, na którym
%                  przybliżamy szukaną funkcję
%              n - liczba naturalna dodatnia (większa od 2), ilość
%                  podziałów przedziału
%         y_alfa - warunek początkowy (wektor warunków) y(alfa) = y_alfa
%                  (dla równań 2-go rzędu y(alfa) = y_alfa(1), y'(alfa) =
%                  y_alfa(2)
%              f - uchwyt do funkcji zmiennej x
%              c - uchwyt do funkcji zmiennej x
%              b - uchwyt do funkcji zmiennej x (dla równań 1-go rzędu
%                  b(x) ~= 0, dla x z przedziału [alfa, beta])
%              a - uchwyt do funkcji zmiennej x, domyślnie 0 (dla równań
%                  1-go rzędu)(dla równań 2-go rzędu a(x) ~= 0, dla x z 
%                  przedziału [alfa, beta])
% Wyjście:
%              y - wektor 1x(n+1), przybliżenia szukanej funkcji w węzłach
%              h - liczba rzeczywista, krok całkowania
%              x - wektor 1x(n+1), węzły równoodległe

rzad_rownania = 2;
if (nargin < 8)
    rzad_rownania = 1; % sprawdzenie rzędu równania
end

% Dla równań pierwszego rzędu
if (rzad_rownania == 1)
    [y, h, x] = AB_1(alfa, beta, n, y_alfa, f, c, b);
    return
end

% Równania drugiego rzędu
h = (beta - alfa) / n ; % krok całkowania
x = alfa:h:beta; % węzły

Y = zeros(2, n + 1);
k = 4; % liczba węzłów, w których wartość będzie policzona metodą RK


F = @(x,Y) [Y(2, :); (f(x) - c(x) .* Y(1, :) - b(x) .* Y(2, :)) ./ a(x)];

% Metoda RK  4-go rzędu do obliczenia Y2, Y3, Y4
Y(:, 1:k) = RK_4TH_GILL(h, k, y_alfa, F, x(1:k));

%AB dla rownan 2 rzędu
for i = 4:n
    Y(:, i + 1) = Y(:, i) + (h / 24) * (55 * F(x(i), Y(:, i)) - 59 * ...
        F(x(i - 1), Y(:, i - 1)) + 37 * F(x(i - 2), Y(:, i - 2)) - 9 ...
        * F(x(i - 3), Y(:, i - 3)));
end

y = Y(1, :);

end % function

