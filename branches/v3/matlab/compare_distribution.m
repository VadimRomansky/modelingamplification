clear;
load distribution.dat;
load fullDistribution.dat;
load coordinateDistribution.dat;
N1=1;
N2=1000;
N3=100;
N5=6;
N7 = 6;

m = 1.67*10^-24;
c = 2.997*10^10;

eta = 0.0003;
rho = 1.6*10^-24;
u1=2900000000;
u2 = u1/4;
p0=distribution(N5,1);
sigma = u1/u2;
gamma = -3*sigma/(sigma - 1);
Q=eta*rho*u1*u1/(c*m*p0*p0);
%Q = 1;
electron = 4.84*10^-10;
B=3*10^-6;


N4 = size(fullDistribution,1);
a3 = fix(N4/N3) - 1;

N6 = size(coordinateDistribution,1);
b3 = fix(N6/N2) - 1;
b2 = fix(0.5*N6/N2) - 1;
b1 = fix(0.25*N6/N2) - 1;

funLocal(1:N3,1:3) = 0;
funCoorDinate(1:N3,1:5) = 0;

x0=0;

for i = 1: N3,
    v = distribution(i,1)/(sqrt(m*m + distribution(i,1)*distribution(i,1)/(c*c)));
    D = distribution(i,1)*v*c/(electron*B);
    funLocal(i,1) = distribution(i,1)/(m*c);
    funLocal(i,2) = distribution(a3*N3 + i,2);
    funLocal(i,3) = (3*Q/(p0*(u1-u2)))*((distribution(i,1)/p0)^gamma)*exp(u1*x0/D);
end

v = distribution(N7,1)/(sqrt(m*m + distribution(N7,1)*distribution(N7,1)/(c*c)));
D = distribution(N7,1)*v*c/(electron*B);
for i = 1: N2,
    funCoorDinate(i,1) = coordinateDistribution(i,1);
    funCoorDinate(i,2) = coordinateDistribution(b3*N2 + i,2);
    funCoorDinate(i,3) = coordinateDistribution(b2*N2 + i,2);
    funCoorDinate(i,4) = coordinateDistribution(b1*N2 + i,2);
    if coordinateDistribution(i,1) < 0
        funCoorDinate(i,5) = (3*Q/(p0*(u1-u2)))*((distribution(N7,1)/p0)^gamma)*exp(u1*coordinateDistribution(i,1)/D);
    else 
        funCoorDinate(i,5) = (3*Q/(p0*(u1-u2)))*((distribution(N7,1)/p0)^gamma);
    end;
end

figure(1);
plot (funLocal(1:N3,1), funLocal(1:N3,2),'blue', funLocal(1:N3,1), funLocal(1:N3,3),'red');
title ('f(p)');
xlabel ('p/mc');
ylabel ('f');
legend('experiment','analityc',4);
grid ;

figure(2);
plot (funCoorDinate(1:N2,1), funCoorDinate(1:N2,2),'blue',funCoorDinate(1:N2,1), funCoorDinate(1:N2,3),'green',funCoorDinate(1:N2,1), funCoorDinate(1:N2,4),'black', funCoorDinate(1:N2,1), funCoorDinate(1:N2,5),'red');
title ('f(p)');
xlabel ('r cm');
ylabel ('f');
legend('numericaly t3','numericaly t2','numericaly t1','analityc',4);
grid ;
