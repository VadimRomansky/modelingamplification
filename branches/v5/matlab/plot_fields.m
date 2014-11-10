clear;
load Efield.dat;
load Bfield.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

NE = Nx*Ny*Nz;
NB = Nx*Ny*Nz;
Nt = size(Efield, 1)/NE;

ynumber = 2;
znumber = 2;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Ex(1:Nx, 1:3) = 0;
Ey(1:Nx, 1:3) = 0;
Ez(1:Nx, 1:3) = 0;

Bx(1:Nx, 1:3) = 0;
By(1:Nx, 1:3) = 0;
Bz(1:Nx, 1:3) = 0;


for i=1:Nx,
   Ex(i,1) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + a*NE, 1);
   Ex(i,2) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + b*NE, 1);
   Ex(i,3) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + c*NE, 1);
   Ey(i,1) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + a*NE, 2);
   Ey(i,2) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + b*NE, 2);
   Ey(i,3) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + c*NE, 2);
   Ez(i,1) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + a*NE, 3);
   Ez(i,2) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + b*NE, 3);
   Ez(i,3) = Efield((Nz)*(Ny)*(i-1) + (Nz)*ynumber + znumber + c*NE, 3);
   
   Bx(i, 1) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + a*NB, 1);
   Bx(i, 2) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + b*NB, 1);
   Bx(i, 3) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + c*NB, 1);
   By(i, 1) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + a*NB, 2);
   By(i, 2) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + b*NB, 2);
   By(i, 3) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + c*NB, 2);
   Bz(i, 1) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + a*NB, 3);
   Bz(i, 2) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + b*NB, 3);
   Bz(i, 3) = Bfield((Nz*Ny*(i-1) + Nz*ynumber + znumber) + c*NB, 3);
end;
figure(1);
plot (Xfile(1:Nx,1),Ex(1:Nx,1), 'red',Xfile(1:Nx,1),Ex(1:Nx,2), 'green',Xfile(1:Nx,1),Ex(1:Nx,3), 'blue');
title ('Ex');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Ey(1:Nx, 1), 'red', Xfile(1:Nx,1), Ey(1:Nx, 2), 'green',Xfile(1:Nx,1),Ey(1:Nx, 3), 'blue');
title ('Ey');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Ez(1:Nx, 1), 'red', Xfile(1:Nx,1), Ez(1:Nx, 2), 'green', Xfile(1:Nx,1), Ez(1:Nx, 3), 'blue');
title ('Ez');
xlabel ('x/r_g');
ylabel ('E gauss');
grid ;

figure(4);
plot (Xfile(1:Nx,1),Bx(1:Nx, 1), 'red', Xfile(1:Nx,1),Bx(1:Nx, 2), 'green', Xfile(1:Nx,1),Bx(1:Nx, 3), 'blue');
title ('Bx');
xlabel ('x/r_g');
ylabel ('B gauss');
grid ;

figure(5);
plot (Xfile(1:Nx,1),By(1:Nx, 1), 'red', Xfile(1:Nx,1),By(1:Nx, 2), 'green', Xfile(1:Nx,1),By(1:Nx, 3), 'blue');
title ('By');
xlabel ('x/r_g');
ylabel ('B gauss');
grid ;

figure(6);
plot (Xfile(1:Nx,1),Bz(1:Nx, 1), 'red', Xfile(1:Nx,1),Bz(1:Nx, 2), 'green', Xfile(1:Nx,1),Bz(1:Nx, 3), 'blue');
title ('Bz');
xlabel ('x/r_g');
ylabel ('B gauss');
grid ;
