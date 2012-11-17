clear;
load tamc_pdf0.dat;
N=100;
pdf(1:N,1:3)=0;
for j=1:N,
    pdf(j,1)=tamc_pdf(j,1);
    pdf(j,2)=tamc_pdf(j,2);
    pdf(j,3)=tamc_pdf(j,3);
end;
figure(1);
plot (pdf(1:N,1),pdf(1:N,2),'blue',pdf(1:N,1),pdf(1:N,3),'red');
title ('pdf');
xlabel ('p');
ylabel ('dN/N');
grid ;
