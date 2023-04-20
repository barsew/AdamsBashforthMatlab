function [y, h, x] = AB_1(alfa, beta, n, y_alfa, f, b, a)
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Funkcja rozwiązująca liniowe równania różniczkowe rzędu 1-go metodą 
% Adamsa-Bashfortha rzedu 4-go
% (pierwsze trzy kroki metodą Rungego-Kutty rzędu 4-go (wzór Gilla))
% równania postaci a(x)y' + b(x)y = f(x)
% Wejście:
%   [alfa, beta] - liczby rzeczywiste, tworzą przedział, na którym
%                  przybliżamy szukaną funkcję
%              n - liczba naturalna dodatnia (większa od 2), ilość
%                  podziałów przedziału
%         y_alfa - warunek początkowy (wektor warunków) y(alfa) = y_alfa
%              f - uchwyt do funkcji zmiennej x
%              b - uchwyt do funkcji zmiennej x
%              a - uchwyt do funkcji zmiennej x (dla równań 1-go rzędu
%                  a(x) ~= 0, dla x z przedziału [alfa, beta])
% Wyjście:
%              y - wektor 1x(n+1), przybliżenia szukanej funkcji w węzłach
%              h - liczba rzeczywista, krok całkowania
%              x - wektor 1x(n+1), węzły równoodległe

h = (beta - alfa) / n;
x = alfa:h:beta;

% alokacja pamieci na y
y = zeros(1, n + 1);
y(1) = y_alfa; % warunek poczatkowy

F = @(x,y) (f(x) - b(x) .* y) ./ a(x);

% metoda RK  4-go rzędu do obliczenia y1, y2, y3:
k = 4;

% współczynniki wykorzystwane w metodzie RK
s1 = (-1 + sqrt(2)) / 2;
s2 = 1 - sqrt(2) / 2;
s3 = -sqrt(2) / 2;
s4 = 1 + sqrt(2) / 2;
s5 = 2 - sqrt(2);
s6 = 2 + sqrt(2);

for i = 1:(k - 1)
    k1 = h .* F(x(i), y(i));
    k2 = h .* F(x(i) + h / 2, y(i) + k1 / 2);
    k3 = h .* F(x(i) + h / 2, y(i) + s1 * k1 + s2 * k2);
    k4 = h .* F(x(i) + h, y(i) + s3 * k2 + s4 * k3);

    y(i + 1) = y(i) + (k1 + s5 * k2 + s6 * k3 + k4) ./ 6;
end

% AB dla rownan 1 rzędu
for i = 4:n
    y(i + 1) = y(i) + (h / 24) * (55 * F(x(i), y(i)) - ...
             59 * F(x(i - 1), y(i - 1)) + 37 * F(x(i - ...
             2), y(i - 2)) - 9 * F(x(i - 3), y(i - 3)));
end

end % function

