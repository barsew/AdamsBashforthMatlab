function [] = num_test1()
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Test numeryczny, wartość błędu globalnego dla bardzo małych h.

alfa = 1;
beta = 1.01;
n = [100, 1000, 10000];
y_alfa = [0; 0];

fprintf("\nRównanie: y'' + y = x^2, y(1) = 0, y'(1) = 0\n");
fprintf("Dokładne rozwiązanie: y = x^2 + 2sin(1 - x) + cos(1 - x) - 2\n");
fprintf("Przedział: [%d, %d]\n\n", alfa, beta);
pause; 

for i = 1:length(n)
    [y, h, x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), y_alfa, ...
        @(x) x.^2, @(x) 1, @(x) 0, @(x) 1);
    fun = x.^2 + 2 .* sin(1 - x) + cos(1 - x) - 2;
    blad = max(abs(y - fun));
    fprintf("\nn = %d, h = %d\n", n(i), h);
    fprintf("Błąd %d: %d\n", i, blad);
    pause;
end

end % function

