function [] = test2()
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Funckja testująca zależność błędu globalnego w zależności od kroku
% całkowania

fprintf("Funkcja sprawdza zależność błędu globalnego od kroku całko" + ...
    "wania\ndla wybranych liniowych równań różniczkowych rzędu 1-go" + ...
    " i 2-go\nna przedziale [alfa, beta], dla n = 100, 200. Dl" + ...
    "a metody\nAdamsa-Bashfortha rzędu 4-go błąd globalny powinien " + ...
    "być około\nrzędu h^4, dla dostatecznie małych h.\n");
pause;

n = [100, 200];
h = zeros(1, length(n));
blad = zeros(1, length(n));

% Równanie 1

alfa = -1;
beta = 1;
fprintf("\nPrzedział: [%d, %d]\n", alfa, beta);
fprintf("Równanie 1: y' + 2y = x, y(-1) = 0\n");
fprintf("Dokładne rozwiązanie: y = (2x + 3e^(-2x - 2) - 1) / 4\n");
pause;
for i = 1:length(n)
    [y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), 0, @(x) ...
                   x, @(x) 2, @(x) 1);
    fun = (2 .* x + 3 .* exp(-2 .* x - 2) - 1) ./ 4;
    blad(i) = max(abs(fun - y));
    fprintf("\nn = %d, h%d = %d\n", n(i), i, h(i));
    fprintf("Błąd %d: %d\n", i, blad(i));
    pause;
end

% Równanie 2

alfa = 1;
beta = 2;
fprintf("\nPrzedział: [%d, %d]\n", alfa, beta);
fprintf("Równanie 2: y'' + y = sinx, y(1) = 0, y'(1) = 0\n");
fprintf("Dokładne rozwiązanie: y = (-xcosx + cosx - cos(1)" + ...
        "sin(1 - x)) / 2\n");
pause;
for i = 1:length(n)
    [y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), [0; 0], ...
                   @(x) sin(x), @(x) 1, @(x) 0, @(x) 1);
    fun = (-x .* cos(x) + cos(x) - cos(1) .* sin(1 - x)) ./ 2;
    blad(i) = max(abs(fun - y));
    fprintf("\nn = %d, h%d = %d\n", n(i), i, h(i));
    fprintf("Błąd %d: %d\n", i, blad(i));
    pause;
end

end % function

