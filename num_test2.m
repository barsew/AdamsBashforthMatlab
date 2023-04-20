function [] = num_test2()
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Test numeryczny, bład globalny dla równań których rozwiązaniem jest
% funkcja, która ma w sobie e^x, wtedy dla dużych x, błędy metody powinny
% być duże.

alfa = 0;
beta = 20;
n = [1000, 2000, 20000];
y_alfa = [1; 1];

fprintf("\nRównanie: y'' - y' = 0, y(0) = 1, y'(0) = 1\n");
fprintf("Dokładne rozwiązanie: y = e^x\n");
fprintf("Przedział: [%d, %d]\n", alfa, beta);
pause;

for i = 1:length(n)
    [y, h, x] = P1Z18_BSE_adams_bashforth(alfa, beta, n(i), y_alfa, ...
        @(x) 0, @(x) 0, @(x) -1, @(x) 1);
    fun = exp(x);
    blad = max(abs(y - fun));
    fprintf("\nn = %d, h = %d\n", n(i), h);
    fprintf("Błąd %d: %d\n", i, blad);
    pause;
end

end % function

