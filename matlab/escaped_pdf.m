clear;
load escaped_pdf.dat;
N=size(escaped_pdf,1);
pdf(1:N,1:3)=0;
for j=1:N,
    pdf(j,1)=escaped_pdf(j,1);
    pdf(j,2)=escaped_pdf(j,2);
    pdf(j,3)=escaped_pdf(j,3);
end;
figure(1);
plot (pdf(1:N,1)/10^4,pdf(1:N,2),'blue',pdf(1:N,1)/10^4,pdf(1:N,3),'red');
title ('pdf');
xlabel ('p  10^-16 g*cm/s');
ylabel ('dN/N');
grid ;
