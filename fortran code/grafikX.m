clear;
load gfull.dat;
load Nx.dat;
load Np.dat;
load x.dat
g0(1:Nx,1:(Np-1))=0;
fteor(1:Nx)=0;
j=0;
for i=1:Nx,
    for k=1:(Np-1), 
        j=k+(i-1)*(Np-1);
        g0(i,k)=gfull(j);
    end
end
for i=1:Nx/2,
    fteor(i)=g0(Nx/2,50)*exp(10*x(i));
end
for i=(Nx/2+1):Nx,
    fteor(i)=g0(Nx/2,50);
end
figure(1);
plot(x(1:Nx),g0(1:Nx,50,1),'g',x(1:Nx),fteor(1:Nx),'r');
title('f');
xlabel('x');
ylabel('f');
grid;
