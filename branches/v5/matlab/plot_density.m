clear;
load concentrations.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;
load divergence_error.dat;
load flux.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

N = Nx*Ny*Nz;
Nt = size(concentrations, 1)/N;

ynumber = 2;
znumber = 2;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

electronConcentration(1:Nx, 1:3) = 0;
protonConcentration(1:Nx, 1:3) = 0;
chargeDensity(1:Nx, 1:3) = 0;
shiftChargeDensity(1:Nx, 1:3) = 0;
divergenceError(1:Nx, 1:3) = 0;
electricFluxX(1:Nx, 1:3) = 0;
electricFluxY(1:Nx, 1:3) = 0;
electricFluxZ(1:Nx, 1:3) = 0;

for i=1:Nx,   
   electronConcentration(i, 1) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 1);
   electronConcentration(i, 2) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 1);
   electronConcentration(i, 3) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 1);
   protonConcentration(i, 1) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 2);
   protonConcentration(i, 2) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 2);
   protonConcentration(i, 3) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 2);
   chargeDensity(i, 1) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 3);
   chargeDensity(i, 2) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 3);
   chargeDensity(i, 3) = concentrations((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 3);
   divergenceError(i, 1) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 2);
   divergenceError(i, 2) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 2);
   divergenceError(i, 3) = divergence_error((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 2);
   electricFluxX(i, 1) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 1);
   electricFluxX(i, 2) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 1);
   electricFluxX(i, 3) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 1);
   electricFluxY(i, 1) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 2);
   electricFluxY(i, 2) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 2);
   electricFluxY(i, 3) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 2);
   electricFluxZ(i, 1) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + a*N, 3);
   electricFluxZ(i, 2) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + b*N, 3);
   electricFluxZ(i, 3) = flux((Nz*Ny*(i-1) + Nz*(ynumber-1) + znumber) + c*N, 3);
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

figure(4);
plot (Xfile(1:Nx,1),divergenceError(1:Nx, 1), 'red', Xfile(1:Nx,1), divergenceError(1:Nx, 2), 'green', Xfile(1:Nx,1), divergenceError(1:Nx, 3), 'blue');
title ('divergence error');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;

figure(5);
plot (Xfile(2:Nx,1),divergenceError(2:Nx, 1), 'red', Xfile(2:Nx,1), divergenceError(2:Nx, 2), 'green', Xfile(2:Nx,1), divergenceError(2:Nx, 3), 'blue');
title ('divergence error');
xlabel ('x/r_g');
ylabel ('rho sgs*cm^-3');
grid ;

figure(6);
plot (Xfile(1:Nx,1), electricFluxX(1:Nx, 1), 'red', Xfile(1:Nx,1), electricFluxX(1:Nx, 2), 'green', Xfile(1:Nx,1), electricFluxX(1:Nx, 3), 'blue');
title ('electric flux');
xlabel ('x/r_g');
ylabel ('Jx');
grid 

figure(7);
plot (Xfile(1:Nx,1), electricFluxY(1:Nx, 1), 'red', Xfile(1:Nx,1), electricFluxY(1:Nx, 2), 'green', Xfile(1:Nx,1), electricFluxY(1:Nx, 3), 'blue');
title ('electric flux');
xlabel ('x/r_g');
ylabel ('Jy');
grid ;

figure(8);
plot (Xfile(1:Nx,1), electricFluxZ(1:Nx, 1), 'red', Xfile(1:Nx,1), electricFluxZ(1:Nx, 2), 'green', Xfile(1:Nx,1), electricFluxZ(1:Nx, 3), 'blue');
title ('electric flux');
xlabel ('x/r_g');
ylabel ('Jz');
grid ;