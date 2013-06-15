clear;
load iterations.dat;
N = size(iterations,1);
figure(1);
plot (iterations(1:N,1),iterations(1:N,3),'red');
title ('upstream1 mass');
xlabel ('ieration');
ylabel ('m g');
grid ;