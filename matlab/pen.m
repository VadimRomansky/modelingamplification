clear;
load u.dat;
N1=1;
N2 = size(u,1);

figure(1);
plot (u(1:N2,1),u(1:N2,2),'red');
title ('u(x)');
xlabel ('r cm');
ylabel ('U cm/s');
grid ;