function [] = test4()
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Funckja testująca błąd globalny w metodzie Rungego-Kutty 4-go rzędu
% (wzór Gilla) i porównująca go z błędem globalnym w metodzie
% Adamsa-Bashfortha 4-go rzędu.

fprintf("Funkcja sprawdza zależność błędu globalnego od kroku " + ...
    "całkowania\ndla liniowego równania różniczkowego rzędu 2-" + ...
    "go na przedziale\n[alfa, beta], dla n = 100, 200 w metodz" + ...
    "ie Rungego-Kutty 4-go\nrzędu (wzór Gilla) i porównuje do " + ...
    "błędu globalnego w metodzie\nAdamsa-Bashfortha 4-go rzędu" + ...
    " dla tego samego równania. Dla me-\ntody Rungego-Kutty bł" + ...
    "ąd globalny powinien być około rzędu h^4,\ndla dostateczn" + ...
    "ie małych h.\n");
pause;

alfa = 1;
beta = 2;
n = [100, 200];
h = (beta - alfa) ./ n;
a = @(x) 1;
b = @(x) 0;
c = @(x) 2;
f = @(x) x.^3;
y_alfa = [0; 0];

F = @(x,Y) [Y(2, :); (f(x) - c(x) .* Y(1, :) - b(x) .* Y(2, :)) ./ a(x)];
% RK
fprintf("\nPrzedział: [%d, %d]\n", alfa, beta);
fprintf("Równanie 1: y'' + 2y = x^3, y(1) = 0, y'(1) = 0\n");
fprintf("Dokładne rozwiązanie: y = (x^3 - 3x + 2cos(sqrt(2)" + ...
    "*(x - 1))) / 2\n");
pause;
fprintf("\nMetoda Rungego-Kutty 4-go rzędu:\n");
for i = 1:length(n)
    x = alfa:h(i):beta;
    Y = RK_4TH_GILL(h(i), n(i) + 1, y_alfa, F, x);
    fun = (x.^3 - 3 .* x + 2 .* cos(sqrt(2) .* (x - 1))) ./ 2;
    blad = max(abs(fun - Y(1, :)));
    fprintf("\nn = %d, h%d = %d\n", n(i), i, h(i));
    fprintf("Błąd %d: %d\n", i, blad);
    pause;
end

% AB
fprintf("\nMetoda Adamsa-Bashfortha 4-go rzędu:\n");
for i = 1:length(n)
    [y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, ...
        n(i), y_alfa, f, c, b, a);
    fun = (x.^3 - 3 .* x + 2 .* cos(sqrt(2) .* (x - 1))) ./ 2;
    fprintf("\nn = %d, h%d = %d\n", n(i), i, h(i));
    blad(i) = max(abs(fun - y));
    fprintf("Błąd %d: %d\n", i, blad(i));
    pause;
end

end % function
    
