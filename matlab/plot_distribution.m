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

N4 = size(fullDistribution,1);
a1 = 1;
a2 = fix(N4/(2*N3));
a3 = fix(N4/N3) - 1;

figure(2);
plot (fullDistribution(1:N3,1),fullDistribution(1:N3,2),'red',fullDistribution(1:N3,1),fullDistribution(a1*N3 + (1:N3),2),'green',fullDistribution(1:N3,1),fullDistribution(a2*N3 + (1:N3),2),'black',fullDistribution(1:N3,1),fullDistribution(a3*N3 + (1:N3),2),'blue');
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