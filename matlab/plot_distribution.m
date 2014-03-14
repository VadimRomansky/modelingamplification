clear;
load distribution.dat;
load fullDistribution.dat;
load coordinateDistribution.dat;
N1=1;
N2=1000;
N3=100;

m = 1.67*10^-24;
c = 3*10^10;


N4 = size(fullDistribution,1);
a1 = fix(N4/(4*N3));
a2 = fix(N4/(2*N3));
a3 = fix(N4/N3) - 1;

N5 = size(coordinateDistribution,1);
b1 = fix(N5/(4*N2));
b2 = fix(N5/(2*N2));
b3 = fix(N5/N2) - 1;

fun(1:N3,1:4) = 0;
funLocal(1:N3,1:4) = 0;

for i = 1: N3,
    fun(i,1) = fullDistribution(i,1)/(m*c);
    fun(i,2) = fullDistribution(a1*N3 + i,2)*(fun(i,1)^4);
    fun(i,3) = fullDistribution(a2*N3 + i,2)*(fun(i,1)^4);
    fun(i,4) = fullDistribution(a3*N3 + i,2)*(fun(i,1)^4);
    
    funLocal(i,1) = distribution(i,1)/(m*c);
    funLocal(i,2) = distribution(a1*N3 + i,2)*(fun(i,1)^4);
    funLocal(i,3) = distribution(a2*N3 + i,2)*(fun(i,1)^4);
    funLocal(i,4) = distribution(a3*N3 + i,2)*(fun(i,1)^4);
end

figure(1);
plot (fun(1:N3,1), fun(1:N3,2),'blue', fun(1:N3,1), fun(1:N3,3),'green', fun(1:N3,1), fun(1:N3,4),'red');
title ('f(p)*p^4');
xlabel ('p/mc');
ylabel ('f');
grid ;


figure(2);
plot (coordinateDistribution(1:N3,1)/(m*c),fullDistribution(1:N3,2),'red',fullDistribution(a1*N3 + (1:N3),1)/(m*c),fullDistribution(a1*N3 + (1:N3),2),'green',fullDistribution(a2*N3 + (1:N3),1)/(m*c),fullDistribution(a2*N3 + (1:N3),2),'black',fullDistribution(a3*N3 + (1:N3),1)/(m*c),fullDistribution(a3*N3 + (1:N3),2),'blue');
title ('f(p)');
xlabel ('p/mc');
ylabel ('f');
grid ;

figure(3);
plot (fun(1:N3,1), fun(1:N3,4),'red');
title ('f(p)*p^4');
xlabel ('p/mc');
ylabel ('f');
grid ;

figure(4);
plot (funLocal(1:N3,1), funLocal(1:N3,2),'blue', funLocal(1:N3,1), funLocal(1:N3,3),'green', funLocal(1:N3,1), funLocal(1:N3,4),'red');
title ('f(p)*p^4');
xlabel ('p/mc');
ylabel ('f');
grid ;

figure(5);
plot (coordinateDistribution(b1*N2 + (1:N2),1), coordinateDistribution(b1*N2 + (1:N2),2),'blue', coordinateDistribution(b2*N2 + (1:N2),1), coordinateDistribution(b2*N2 + (1:N2),2),'green', coordinateDistribution(b3*N2 + (1:N2),1), coordinateDistribution(b3*N2 + (1:N2),2),'red');
title ('f(p = const, r)');
xlabel ('r');
ylabel ('f');
grid ;