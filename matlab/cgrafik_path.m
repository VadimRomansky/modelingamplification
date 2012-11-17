clear;
load tamc_cosmic_ray_path.dat;
load tamc_not_cosmic_ray_path.dat;
N1 = tamc_cosmic_ray_path(1);
N2 = tamc_not_cosmic_ray_path(1);
cosmic_ray_path(1:N1,1:2) = 0;
not_cosmic_ray_path(1:N2,1:2) = 0;
for i = 1:N1,
    cosmic_ray_path(i,1) = i;
    cosmic_ray_path(i,2) = tamc_cosmic_ray_path(i+1);
end;

for i = 1:N2,
    not_cosmic_ray_path(i,1) = i;
    not_cosmic_ray_path(i,2) = tamc_not_cosmic_ray_path(i+1);
end;
figure(1);
plot (cosmic_ray_path(1:N1,2),cosmic_ray_path(1:N1,1),'blue');
title ('cosmic ray path');
xlabel ('x');
ylabel ('step');
grid ;

figure(2);
plot (not_cosmic_ray_path(1:N2,2),not_cosmic_ray_path(1:N2,1),'blue');
title ('not cosmic ray path');
xlabel ('x');
ylabel ('step');
grid ;