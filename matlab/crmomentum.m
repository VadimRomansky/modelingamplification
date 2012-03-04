clear;
load tamc_cosmic_ray_momentum.dat;
N1=1;
N2=100;
figure(1);
plot (tamc_cosmic_ray_momentum(1:N2,1),tamc_cosmic_ray_momentum(1:N2,2),'blue');
title ('cosmic ray bound momentum');
xlabel ('x');
ylabel ('p');
grid 
figure(2);
plot (tamc_cosmic_ray_momentum(1:N2,1),tamc_cosmic_ray_momentum(1:N2,3),'blue');
title ('cosmic ray flux');
xlabel ('x');
ylabel ('j');
grid ;