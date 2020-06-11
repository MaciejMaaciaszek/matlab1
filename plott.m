f = @(x,y, dx, dy) 1 ./ (1 + (x-dx).^2 + (y-dy).^2);
g = @(x,y) f(x,y, 0, 0) + 0.5 * f(x, y, 4, 0);

x = -6 : 0.1 : 6;
y = -6 : 0.1 : 6;

[xx, yy] = meshgrid(x,y);

zz = g(xx, yy);

subplot(3,3, [1,2]);
surf(xx, yy, zz);

contour (xx, yy, zz);
