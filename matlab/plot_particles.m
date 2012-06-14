clear;
load particles.dat;
N=size(particles,1);

figure(1);
plot (particles(1:N,1),particles(1:N,2),'blue');
title ('particles');
xlabel ('x cm');
ylabel ('p g*cm/s');
grid ;