clear;
load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

N = Nx*Ny*Nz;
Nt = size(concentrations, 1)/N;

a = 0;
b = fix(Nt/2)-2;
c = Nt - 1;


electronConcentration(1:Nx, 1:3) = 0;
protonConcentration(1:Nx, 1:3) = 0;
chargeDensity(1:Nx, 1:3) = 0;


for i=1:Nx,   
   electronConcentration(i, 1) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + a*N, 1);
   electronConcentration(i, 2) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + b*N, 1);
   electronConcentration(i, 3) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + c*N, 1);
   protonConcentration(i, 1) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + a*N, 2);
   protonConcentration(i, 2) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + b*N, 2);
   protonConcentration(i, 3) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + c*N, 2);
   chargeDensity(i, 1) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + a*N, 3);
   chargeDensity(i, 2) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + b*N, 3);
   chargeDensity(i, 3) = concentrations((Nz*Ny*(i-1) + Nz*3 + 1) + c*N, 3);
end;
figure(1);
plot (Xfile(1:Nx,1),electronConcentration(1:Nx,1), 'red',Xfile(1:Nx,1),electronConcentration(1:Nx,2), 'green',Xfile(1:Nx,1),electronConcentration(1:Nx,3), 'blue');
title ('N_e');
xlabel ('x/r_g');
ylabel ('N_e cm^-3');
grid ;

figure(2);
plot (Xfile(1:Nx,1),protonConcentration(1:Nx, 1), 'red', Xfile(1:Nx,1), protonConcentration(1:Nx, 2), 'green',Xfile(1:Nx,1),protonConcentration(1:Nx, 3), 'blue');
title ('N_p');
xlabel ('x/r_g');
ylabel ('N_p cm^-3');
grid ;

figure(3);
plot (Xfile(1:Nx,1),chargeDensity(1:Nx, 1), 'red', Xfile(1:Nx,1), chargeDensity(1:Nx, 2), 'green', Xfile(1:Nx,1), chargeDensity(1:Nx, 3), 'blue');
title ('rho');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;