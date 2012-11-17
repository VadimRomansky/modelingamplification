clear;
load pdf_mid.dat;
N=100;
pdfm(1:N,1:6)=0;
for j=1:N,
    pdfm(j,1)=pdf_mid(2*j-1,3);
    pdfm(j,2)=pdf_mid(5*2*N+2*j-1,4);
    pdfm(j,3)=pdf_mid(6*2*N+2*j-1,4);
    pdfm(j,4)=pdf_mid(7*2*N+2*j-1,4);
    pdfm(j,5)=pdf_mid(8*2*N+2*j-1,4);
    pdfm(j,6)=pdf_mid(15*2*N+2*j-1,4);
end;
figure(1);
plot (pdfm(1:N,1),pdfm(1:N,6),'blue');
title ('pdf mib');
xlabel ('p');
ylabel ('');
grid ;