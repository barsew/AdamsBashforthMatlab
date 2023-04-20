% testy 
clearvars;

n = 10000000; 
a = @(x) 1;
b = @(x) 0;
c = @(x) -1;
f = @(x) x;
y_alfa = [1, 1];
% y_alfa = 1;
alfa = 0;
beta = 10;

% h = (beta - alfa) / n;
% x = alfa:h:beta;
% F = @(x,y) [a(x); f(x) - b(x)];
% F([3, 2], 3)

[y, h ,x] = P1Z18_BSE_adams_bashforth(alfa, beta, n, y_alfa, f, c, b, a);
%[y, h, x] = AB_test(alfa, beta, n, a, b, f, y_alfa);
%[y, h] = rk_3_8(alfa, beta, n, a, b, f, y_alfa);
%[y, h] = RK_3_8_NORMAL(alfa, beta, n, F, y_alfa);
h

%fun = (1/2)*(3*exp(-x) + sin(x) - cos(x));
%fun = exp(x);
%fun = x + cos(x);
%fun = -1 + x + 2.*exp(-x/2).*cos((sqrt(3).*x)/2) + (2.*exp(-x/2).*sin((sqrt(3).*x)/2))/sqrt(3);
%fun = (sqrt(3)*sin(sqrt(3).*x/2)+cos(sqrt(3).*x./2))./exp(x/2);
%fun = (x.^4)./12 + x + 1;
%fun = exp(x);
fun = -x - exp(-x) ./ 2 + 3 .* exp(x) ./ 2;
blad = (abs(y - fun));

figure
hold on
plot(x, fun);
 plot(x, y);
legend('analitical', 'my-fun');
hold off
 
figure
plot(x, blad);

