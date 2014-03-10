clear;
load distribution.dat;
load fullDistribution.dat;
load coordinateDistribution.dat;
N1=1;
N2=1000;
N3=100;

m = 6.24*10^-24;
c = 3*10^10;


N4 = size(fullDistribution,1);
a1 = fix(N4/(4*N3));
a2 = fix(N4/(2*N3));
a3 = fix(N4/N3) - 1;

fun(1:N3,1:4) = 0;

for i = 1: N3,
    fun(i,1) = fullDistribution(a3*N3 + i,1)/(m*c);
    fun(i,2) = fullDistribution(a1*N3 + i,2)*(fun(i,1)^4);
    fun(i,3) = fullDistribution(a2*N3 + i,2)*(fun(i,1)^4);
    fun(i,4) = fullDistribution(a3*N3 + i,2)*(fun(i,1)^4);
end

figure(1);
plot (fun(1:N3,1), fun(1:N3,2),'blue', fun(1:N3,1), fun(1:N3,3),'green', fun(1:N3,1), fun(1:N3,4),'red');
title ('f(p)*p^4');
xlabel ('p/mc');
ylabel ('f');
grid ;


figure(2);
plot (fullDistribution(1:N3,1)/(m*c),fullDistribution(1:N3,2),'red',fullDistribution(a1*N3 + (1:N3),1)/(m*c),fullDistribution(a1*N3 + (1:N3),2),'green',fullDistribution(a2*N3 + (1:N3),1)/(m*c),fullDistribution(a2*N3 + (1:N3),2),'black',fullDistribution(a3*N3 + (1:N3),1)/(m*c),fullDistribution(a3*N3 + (1:N3),2),'blue');
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
