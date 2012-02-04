clear;
load tamc_pdf.dat;
N=50;
pdf(1:N,1:3)=0;
for j=1:N,
    pdf(j,1)=tamc_pdf(j,1);
    pdf(j,2)=tamc_pdf(j,2);
    pdf(j,3)=tamc_pdf(j,3);
end;
figure(1);
plot (pdfm(1:N,1),pdfm(1:N,2),'blue',pdfm(1:N,1),pdfm(1:N,3),'red');
title ('pdf');
xlabel ('p');
ylabel ('dN/N');
grid ;
