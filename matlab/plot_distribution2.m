clear;
load distribution2.dat;
load xfile.dat;
load pfile.dat;
N1=1;
N2=size(xfile,1);
N3=size(pfile,1);
m = 1.67*10^-24;
c = 3*10^10;


figure(1);
plot (distribution2(1:N3,1), distribution2(1:N3,2),'blue', distribution2(1:N3,1), distribution2(1:N3,3),'red');
title ('f(p)*p^4');
xlabel ('p/mc');
ylabel ('f');
legend('без усиления','с усилением',1);
grid ;
