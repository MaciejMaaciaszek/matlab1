% Przedmiot: Techniki Obliczeniowe 
% Kierunek studiów: Mechatronika 
% Semestr: 2
% Rok akademicki: 2019/2020
% Data (dzień-miesiąc-rok): <<22.06.2020>>
%
% Imię:             <<Maciej>>
% Nazwisko:         <<Maciaszek>>
% Numer albumu ZUT: <<46759>>
f = @(x,y, dx, dy) 1 ./ (1 + (x-dx).^2 + (y-dy).^2);
g = @(x,y) f(x,y, 0, 0) + 0.5 * f(x, y, 4, 0);

x = -6 : 0.1 : 6;
y = -6 : 0.1 : 6;

[xx, yy] = meshgrid(x,y);

zz = g(xx, yy);

subplot(3,3, [1,2]);
surf(xx, yy, zz);

contour (xx, yy, zz);
