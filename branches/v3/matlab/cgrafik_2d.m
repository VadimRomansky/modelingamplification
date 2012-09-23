clear;
load density2d.dat;
load velocity2d.dat;

N1 = size(density2d,1);
N2 = size(density2d,2);
[X,Y]=meshgrid(1:1:N1, 1:1:N2);   

figure(1);
surf(X,Y,velocity2d);

figure(2);
surf(X,Y,density2d);