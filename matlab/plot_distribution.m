clear;
load distribution.dat;
N1=1;
N2=1000;
N3=1000;
e = N2 - 5;
d = fix(3*N2/4);
c = fix(N2/2);
b = fix(N2/4);
a = 5;

figure(1);
plot (distribution(1:N3,1),distribution(1:N3,2),'cyan',distribution(a*N3 + 1:N3,1),distribution(a*N3 + 1:N3,2),'green',distribution(b*N3 + 1:N3,1),distribution(b*N3 + 1:N3,2),'blue',distribution(c*N3 + 1:N3,1),distribution(c*N3 + 1:N3,2),'black',distribution(d*N3 + 1:N3,1),distribution(d*N3 + 1:N3,2),'yellow',distribution(e*N3 + 1:N3,1),distribution(e*N3 + 1:N3,2),'red');
title ('u(x)');
xlabel ('r cm');
ylabel ('U cm/s');
grid ;