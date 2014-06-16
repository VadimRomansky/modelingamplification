clear;
load distribution.dat;
load fullDistribution.dat;
load coordinateDistribution.dat;
load xfile.dat;
load pfile.dat;
N1=1;
N2=size(xfile,1);
N3=size(pfile,1);

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
hzfun(1:N3,1:4) = 0;
J(1:N3,1:4) = 0;

for i = 1: N3,
    fun(i,1) = fullDistribution(i,1)/(m*c);
    fun(i,2) = fullDistribution(a1*N3 + i,2)*(fullDistribution(i,1)^4);
    fun(i,3) = fullDistribution(a2*N3 + i,2)*(fullDistribution(i,1)^4);
    fun(i,4) = fullDistribution(a3*N3 + i,2)*(fullDistribution(i,1)^4);
    
    funLocal(i,1) = distribution(i,1)/(m*c);
    funLocal(i,2) = distribution(a1*N3 + i,2)*(distribution(i,1)^4);
    funLocal(i,3) = distribution(a2*N3 + i,2)*(distribution(i,1)^4);
    funLocal(i,4) = distribution((a3-1)*N3 + i,2)*(distribution(i,1)^4);
    
    hzfun(i,1) = fullDistribution(i,1)/(m*c);
    hzfun(i,2) = fullDistribution(a1*N3 + i,2)*(fullDistribution(i,1)^7);
    hzfun(i,3) = fullDistribution(a2*N3 + i,2)*(fullDistribution(i,1)^7);
    hzfun(i,4) = fullDistribution(a3*N3 + i,2)*(fullDistribution(i,1)^7); 
    
    J(i,1) = distribution(i,1)/(m*c);
    J(i,2) = distribution(a1*N3 + i,3);
    J(i,3) = distribution(a2*N3 + i,3);
    J(i,4) = distribution((a3-1)*N3 + i,3);
end

figure(1);
plot (fun(1:N3,1), fun(1:N3,2),'blue', fun(1:N3,1), fun(1:N3,3),'green', fun(1:N3,1), fun(1:N3,4),'red');
title ('f(p)*p^4');
xlabel ('p/mc');
ylabel ('f');
grid ;


%figure(2);
%plot (fullDistribution(1:N3,1)/(m*c),fullDistribution(1:N3,2),'cyan',fullDistribution(a1*N3 + (1:N3),1)/(m*c),fullDistribution(a1*N3 + (1:N3),2),'blue',fullDistribution(a2*N3 + (1:N3),1)/(m*c),fullDistribution(a2*N3 + (1:N3),2),'green',fullDistribution(a3*N3 + (1:N3),1)/(m*c),fullDistribution(a3*N3 + (1:N3),2),'red');
%title ('f(p)');
%xlabel ('p/mc');
%ylabel ('f');
%grid ;

%figure(3);
%plot (fun(1:N3,1), fun(1:N3,4),'red');
%title ('f(p)*p^4');
%xlabel ('p/mc');
%ylabel ('f');
%grid ;

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

%figure(6);
%plot (fun(8:N3,1), fun(8:N3,2),'blue', fun(8:N3,1), fun(8:N3,3),'green', fun(8:N3,1), fun(8:N3,4),'red');
%title ('f(p)*p^4');
%xlabel ('p/mc');
%ylabel ('f');
%grid ;

%figure(7);
%plot (hzfun(1:N3,1), hzfun(1:N3,2),'blue', hzfun(1:N3,1), hzfun(1:N3,3),'green', hzfun(1:N3,1), hzfun(1:N3,4),'red');
%title ('f(p)*p^7');
%xlabel ('p/mc');
%ylabel ('f');
%grid ;

figure(8);
plot (coordinateDistribution(b1*N2 + (1:N2),1), coordinateDistribution(b1*N2 + (1:N2),3),'blue', coordinateDistribution(b2*N2 + (1:N2),1), coordinateDistribution(b2*N2 + (1:N2),3),'green', coordinateDistribution(b3*N2 + (1:N2),1), coordinateDistribution(b3*N2 + (1:N2),3),'red');
title ('J(p = const, r)');
xlabel ('r');
ylabel ('J');
grid ;

figure(9);
plot (coordinateDistribution(b1*N2 + (1:N2),1), coordinateDistribution(b1*N2 + (1:N2),4),'blue', coordinateDistribution(b2*N2 + (1:N2),1), coordinateDistribution(b2*N2 + (1:N2),4),'green', coordinateDistribution(b3*N2 + (1:N2),1), coordinateDistribution(b3*N2 + (1:N2),4),'red');
title ('full J(r)');
xlabel ('r');
ylabel ('J');
grid ;

figure(10);
plot (J(1:N3,1), J(1:N3,2),'blue', J(1:N3,1), J(1:N3,3),'green', J(1:N3,1), J(1:N3,4),'red');
title ('J(p)');
xlabel ('p/mc');
ylabel ('J');
grid ;