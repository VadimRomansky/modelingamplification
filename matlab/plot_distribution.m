clear;
load distribution.dat;
load fullDistribution.dat;
load coordinateDistribution.dat;
N1=1;
N2=1000;
N3=100;
e = N2 - 1;
d = fix(3*N2/4);
c = fix(N2/2);
b = fix(N2/4);
a = 5;

figure(1);
plot (distribution(1:N3,1),distribution(1:N3,2),'cyan',distribution(a*N3 + 1:N3,1),distribution(a*N3 + 1:N3,2),'green',distribution(b*N3 + 1:N3,1),distribution(b*N3 + 1:N3,2),'blue',distribution(c*N3 + 1:N3,1),distribution(c*N3 + 1:N3,2),'black',distribution(d*N3 + 1:N3,1),distribution(d*N3 + 1:N3,2),'yellow',distribution(e*N3 + 1:N3,1),distribution(e*N3 + 1:N3,2),'red');
title ('f(p)');
xlabel ('p g*cm/s');
ylabel ('f');
grid ;

figure(2);
plot (fullDistribution(1:N3,1),fullDistribution(1:N3,2),'red');
title ('f(p)');
xlabel ('p g*cm/s');
ylabel ('f');
grid ;

figure(3);
plot (coordinateDistribution(1:N2,1),coordinateDistribution(1:N2,2),'red');
title ('f(p)');
xlabel ('r cm');
ylabel ('f');
grid ;