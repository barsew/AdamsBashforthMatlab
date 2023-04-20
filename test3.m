function [] = test3()
% Projekt 1, zadanie 18
% Bartosz Seweryn, 320733
%
% Funckja sprawdzająca błąd globalny dla równań rózniczkowych, których
% rozwiązaniami są wielomiany różnych stopni.
 
fprintf("Funckja sprawdzająca błąd globalny dla równań różniczkowych" + ...
    ",\nktórych rozwiązaniami są wielomiany różnych stopni. Metoda\n" + ...
    "Adamsa-Bashfortha rzędu 4-go powinna być dokładna dla równań,\n" + ...
    "których rozwizązaniami są wielomiany stopnia mniejszego, bądź\n" + ...
    "równego 4.\n");
pause;

n = 100;
alfa = 0;
beta = 1;
[~, h, x] = P1Z18_BSE_adams_bashforth(alfa, beta, n, 1, @(x) 0, @(x) 0, ...
            @(x) 1);
rownania = ["y' = 0, y(0) = 1", "y' = 1, y(0) = 1", "y' = x, y(0) = 1", ...
    "y' = x^2, y(0) = 1", "y' = x^3, y(0) = 1", "y' = x^4, y(0) = 1"];
rozwiazania = ["y = 1", "y = x + 1", "y = (x^2 + 2) / 2", "(x^3 + 3)" + ...
    " / 3", "(x^4 + 4) / 4", "(x.^5 + 5) / 5"];
f = {@(x) 0, @(x) 1, @(x) x, @(x) x.^2, @(x) x.^3, @(x) x.^4};
funkcje = [ones(1, length(x)); x + 1; (x .* x + 2) ./ 2; (x.^3 + 3) ... 
./ 3; (x.^4 + 4) ./ 4; (x.^5 + 5) ./ 5]; 

for i = 1:length(rownania)
    fprintf("\nRównanie %d: %s\n", i, rownania(i));
    fprintf("Dokładne rozwiązanie: %s\n", rozwiazania(i));
    fprintf("Przedział: [%d, %d]\n", alfa, beta);
    [y, ~, ~] = P1Z18_BSE_adams_bashforth(alfa, beta, n, 1, ...
                cell2mat(f(i)), @(x) 0, @(x) 1);
    blad = max(abs(y - funkcje(i, :)));
    fprintf("\nn = %d, h = %d\n", n, h);
    fprintf("Błąd %d: %d\n", i, blad);
    pause;
end

end % function

