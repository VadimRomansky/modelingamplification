clear;
load velocity.dat;
load Xfile.dat;
load Yfile.dat;
load Zfile.dat;

Nx = size(Xfile, 1);
Ny = size(Yfile, 1);
Nz = size(Zfile, 1);

N = Nx*Ny*Nz;
Nt = size(velocity, 1)/N;
ynumber = 2;
znumber = 2;

a = 0;
b = fix(Nt/2);
c = Nt - 1;

Vx(1:Nx, 1:3) = 0;
Vy(1:Nx, 1:3) = 0;
Vz(1:Nx, 1:3) = 0;


for i=1:Nx,
   Vx(i,1) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 1);
   Vx(i,2) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 1);
   Vx(i,3) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 1);
   Vy(i,1) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 2);
   Vy(i,2) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 2);
   Vy(i,3) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 2);
   Vz(i,1) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + a*N, 3);
   Vz(i,2) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + b*N, 3);
   Vz(i,3) = velocity((Nz)*(Ny)*(i-1) + (Nz)*(ynumber-1) + znumber + c*N, 3);
end;
figure(1);
plot (Xfile(1:Nx,1),Vx(1:Nx,1), 'red',Xfile(1:Nx,1),Vx(1:Nx,2), 'green',Xfile(1:Nx,1),Vx(1:Nx,3), 'blue');
title ('Vx');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(2);
plot (Xfile(1:Nx,1),Vy(1:Nx, 1), 'red', Xfile(1:Nx,1), Vy(1:Nx, 2), 'green',Xfile(1:Nx,1),Vy(1:Nx, 3), 'blue');
title ('Vy');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

figure(3);
plot (Xfile(1:Nx,1),Vz(1:Nx, 1), 'red', Xfile(1:Nx,1), Vz(1:Nx, 2), 'green', Xfile(1:Nx,1), Vz(1:Nx, 3), 'blue');
title ('Vz');
xlabel ('x/r_g');
ylabel ('V cm/s');
grid ;

