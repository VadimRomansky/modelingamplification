clear;
load tamc_field.dat;
N=100;
figure(1);
plot (tamc_field(1:N,1),tamc_field(1:N,2),'blue');
title ('effective field');
xlabel ('x');
ylabel ('B(x)');
grid ;