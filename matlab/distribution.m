clear;
load tamc_distribution.dat;
N=100;
N2 = 6;
e = size(tamc_distribution,1)/N - 1;
d = fix(3*e/4);
c = fix(e/2);
b = fix(e/4);
a = 5;
pdf(1:N,1:N2 + 1)=0;
for j=1:N,
    pdf(j,1)=tamc_distribution(j,1);
    pdf(j,2)=tamc_distribution(j,2);
    pdf(j,3)=tamc_distribution(a*N + j,2);
    pdf(j,4)=tamc_distribution(b*N + j,2);
    pdf(j,5)=tamc_distribution(c*N + j,2);
    pdf(j,6)=tamc_distribution(d*N + j,2);
    pdf(j,7)=tamc_distribution(e*N + j,2);
end;
figure(1);
plot (pdf(1:N,1)/10^4,pdf(1:N,2),'cyan',pdf(1:N,1)/10^4,pdf(1:N,3),'green',pdf(1:N,1)/10^4,pdf(1:N,4),'blue',pdf(1:N,1)/10^4,pdf(1:N,5),'black',pdf(1:N,1)/10^4,pdf(1:N,6),'yellow',pdf(1:N,1)/10^4,pdf(1:N,7),'red');
title ('pdf');
xlabel ('p  10^-16 g*cm/s');
ylabel ('dN/N');
grid ;