function Y = RK_4TH_GILL(h, k, y_alfa, F, x)
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Funkcja rozwiązująca równania różniczkowe rzędu 2-go, sprowadzone do
% układu równań postaci:
% Y = [y; y'], Y' = F(x,Y) metodą Rungego-Kutty rzędu 4-go (wzór Gilla)
% Wejście:
%       h - liczba rzeczywista, krok całkowania
%       k - liczba naturaln, ilość węzłów
%  y_alfa - wektor rozmiaru 2x1, warunki początkowe równania,
%           y(alfa) = y_alfa(1), y'(alfa) = y_alfa(2)
%       F - uchwyt do funkcji dwóch zmiennych, musi ona zwracać wektor
%       rozmiaru 2x1
%       x - wektor rozmiaru 1xk, węzły
% Wyjście:
%       Y - macierz rozmiaru 2xk, wartości y (pierwszy wiersz) i y'
%           (drugi wiersz) w wezłach x

% alokacja pamięci
Y = zeros(2, k);
Y(:, 1) = y_alfa;

% współczynniki wykorzystwane w metodzie
s1 = (-1 + sqrt(2)) / 2;
s2 = 1 - sqrt(2) / 2;
s3 = -sqrt(2) / 2;
s4 = 1 + sqrt(2)/2;
s5 = 2 - sqrt(2);
s6 = 2 + sqrt(2);

for i = 1:(k - 1)
    k1 = h .* F(x(i), Y(:, i));
    k2 = h .* F(x(i) + h / 2, Y(:, i) + k1 / 2);
    k3 = h .* F(x(i) + h / 2, Y(:, i) + s1 * k1 + s2 * k2);
    k4 = h .* F(x(i) + h, Y(:, i) + s3 * k2 + s4 * k3);

    Y(:, i + 1) = Y(:, i) + (k1 + s5 * k2 + s6 * k3 + k4) ./ 6;
end

end % function

