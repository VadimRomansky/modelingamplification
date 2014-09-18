clear;
load Efield.dat;
load Bfield.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

Ex(1:Nx) = 0;
Ey(1:Nx) = 0;
Ez(1:Nx) = 0;

Bx(1:Nx) = 0;
By(1:Nx) = 0;
Bz(1:Nx) = 0;


for i=1:Nx,
   Ex(i) = Efield((Nz+1)*(Ny+1)*(i-1) + (Nz+1)*3 + 1, 1);
   Ey(i) = Efield((Nz+1)*(Ny+1)*(i-1) + (Nz+1)*3 + 1, 2);
   Ez(i) = Efield((Nz+1)*(Ny+1)*(i-1) + (Nz+1)*3 + 1, 3);
   
   Bx(i) = Bfield(Nz*Ny*(i-1) + Nz*3 + 1, 1);
   By(i) = Bfield(Nz*Ny*(i-1) + Nz*3 + 1, 2);
   Bz(i) = Bfield(Nz*Ny*(i-1) + Nz*3 + 1, 3);
end;
figure(1);
plot (Xfile(1:Nx,1),Ex(1:Nx), 'red');
title ('Ex');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Ey(1:Nx), 'red');
title ('Ey');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Ez(1:Nx), 'red');
title ('Ez');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(4);
plot (Xfile(1:Nx,1),Bx(1:Nx), 'red');
title ('Bx');
xlabel ('x/r_g');
ylabel ('B gauss');
grid ;

figure(5);
plot (Xfile(1:Nx,1),By(1:Nx), 'red');
title ('By');
xlabel ('x/r_g');
ylabel ('B gauss');
grid ;

figure(6);
plot (Xfile(1:Nx,1),Bz(1:Nx), 'red');
title ('Bz');
xlabel ('x/r_g');
ylabel ('B gauss');
grid ;

