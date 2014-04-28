clear;
load temp_grid.dat;
N1=1;
N2=1000;
N3=6;
e = size(temp_grid,1)/N2 - 1;
d = fix(3*e/4);
c = fix(e/2);
b = fix(e/4);
a = fix(e/8);

figure(1);
plot (temp_grid(1:N2,1),temp_grid(a*N2 + (1:N2),2),'cyan',temp_grid(1:N2,1),temp_grid(b*N2 + (1:N2),2),'green',temp_grid(1:N2,1),temp_grid(c*N2 + (1:N2),2),'blue',temp_grid(1:N2,1),temp_grid(d*N2 + (1:N2),2),'black',temp_grid(1:N2,1),temp_grid(e*N2 + (1:N2),2),'red');
title ('r');
xlabel ('i');
ylabel ('r cm');
grid ;

figure(2);
plot (temp_grid(1:N2,2), temp_grid(1:N2,1),'cyan',temp_grid(b*N2 + (1:N2),2),temp_grid(1:N2,1),'green',temp_grid(c*N2 + (1:N2),2),temp_grid(1:N2,1),'blue',temp_grid(d*N2 + (1:N2),2),temp_grid(1:N2,1),'black',temp_grid(e*N2 + (1:N2),2),temp_grid(1:N2,1),'red');
title ('r');
xlabel ('r cm');
ylabel ('i');
grid ;