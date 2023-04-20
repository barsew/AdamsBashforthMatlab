function [] = test1()
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Funckja testująca zależność błędu od kroku całkowania

% Przedział na którym przybliżamy funkjcę, [alfa, beta]
alfa = 0;
beta = 1;
n = [100, 200, 400];
h = zeros(1, length(n));
blad = zeros(1, length(n));

fprintf(['Sprawdzenie zależności pomiędzy błędem globalnym metody,' ...
    ' a krokiem \ncałkowania, dla wybranych równań 1-go i 2-go rzę' ...
    'du na przedziale\n[0, 1] oraz dla następujących wartości n: 1' ...
    '00, 200, 400. W metodzie\nAdamsa-Bashfortha rzędu 4-go x-krot' ...
    'ne zmniejszenie kroku całkowania\npowinno skutkować około x^4' ...
    '-krotnym zmniejszeniem wartości błędu \nglobalnego dla dostat' ...
    'ecznie małych h.\n\n']);

% Równanie 1

pause;
fprintf("Równanie 1: y' - y = x, y(0) = 0\n");
fprintf("Dokładne rozwiązanie: y = (2x + e^(-2x) - 1) / 4\n\n");
pause;
for i = 1:length(n)
    [y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), 0, @(x) ...
                   x, @(x) 2, @(x) 1);
    fun = (2 .* x + exp(-2 .* x) - 1) ./ 4;
    fprintf("n = %d, h%d = %d\n", n(i), i, h(i));
    blad(i) = max(abs(fun - y));
    fprintf("Błąd %d: %d\n\n", i, blad(i));
    pause;
end
fprintf("h1 / h2 = %d, blad1 / blad2 = %d\n", h(1) / h(2), ...
    blad(1) / blad(2));
fprintf("h1 / h3 = %d, blad1 / blad3 = %d\n", h(1) / h(3), ...
    blad(1) / blad(3));
fprintf("h2 / h3 = %d, blad2 / blad3 = %d\n", h(2) / h(3), ...
    blad(2) / blad(3));

% Równanie 2

pause;
fprintf("\nRównanie 2: y'' + y' = x, y(0) = 1, y'(0) = 1\n");
fprintf("Dokładne rozwiązanie: y = x^2/2 - x - 2e^(-x) + 3\n\n");
pause;
for i = 1:length(n)
    [y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), [1; 1], ...
                   @(x) x, @(x) 0, @(x) 1, @(x) 1);
    fun = x .* x ./ 2 - x - 2 .* exp(-x) + 3;
    fprintf("n = %d, h%d = %d\n", n(i), i, h(i));
    blad(i) = max(abs(fun - y));
    fprintf("Błąd %d: %d\n\n", i, blad(i));
    pause;
end
fprintf("h1 / h2 = %d, blad1 / blad2 = %d\n", h(1) / h(2), ...
    blad(1) / blad(2));
fprintf("h1 / h3 = %d, blad1 / blad3 = %d\n", h(1) / h(3), ...
    blad(1) / blad(3));
fprintf("h2 / h3 = %d, blad2 / blad3 = %d\n", h(2) / h(3), ...
    blad(2) / blad(3));

% Równanie 3

pause;
fprintf("\nRównanie 3: y'' + y = xsinx, y(0) = 0, y'(0) = 0\n");
fprintf("Dokładne rozwiązanie: y = -x(xcosx - sinx) / 4\n\n");
pause;
for i = 1:length(n)
[y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), [0; 0], @(x) ...
               x .* sin(x), @(x) 1, @(x) 0, @(x) 1);
    fun = -x .* (x .* cos(x) - sin(x)) ./ 4;
    fprintf("n = %d, h%d = %d\n", n(i), i, h(i));
    blad(i) = max(abs(fun - y));
    fprintf("Błąd %d: %d\n\n", i, blad(i));
    pause;
end
fprintf("h1 / h2 = %d, blad1 / blad2 = %d\n", h(1) / h(2), ...
    blad(1) / blad(2));
fprintf("h1 / h3 = %d, blad1 / blad3 = %d\n", h(1) / h(3), ...
    blad(1) / blad(3));
fprintf("h2 / h3 = %d, blad2 / blad3 = %d\n", h(2) / h(3), ...
    blad(2) / blad(3));

% Równanie 4

pause;
fprintf("\nRównanie 4: y'' - y' = 1, y(0) = 1, y'(0) = 1\n");
fprintf("Dokładne rozwiązanie: y = -x + 2*e^x - 1\n\n");
pause;
for i = 1:length(n)
    [y, h(i), x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), [1; 1], ...
                   @(x) 1, @(x) 0, @(x) -1, @(x) 1);
    fun = -x + 2 .* exp(x) - 1;
    fprintf("n = %d, h%d = %d\n", n(i), i, h(i));
    blad(i) = max(abs(fun - y));
    fprintf("Błąd %d: %d\n\n", i, blad(i));
    pause;
end
fprintf("h1 / h2 = %d, blad1 / blad2 = %d\n", h(1) / h(2), ...
    blad(1) / blad(2));
fprintf("h1 / h3 = %d, blad1 / blad3 = %d\n", h(1) / h(3), ...
    blad(1) / blad(3));
fprintf("h2 / h3 = %d, blad2 / blad3 = %d\n", h(2) / h(3), ...
    blad(2) / blad(3));

end % function

