clear;
load tamc_radial_profile.dat;
N1=1;
N2=100;
N3=6;
e = size(tamc_radial_profile,1)/N2 - 1;
d = fix(3*e/4);
c = fix(e/2);
b = fix(e/4);
a = 5;
r(1:N2,1:N3) = 0;
vel(1:N2,1:N3 + 1) = 0;
rho(1:N2,1:N3 + 1)= 0;
temp(1:N2,1:N3 + 1) = 0;
pressure(1:N2,1:N3 + 1) = 0;
crpressure(1:N2,1:N3 + 1) = 0;
flux(1:N2,1:N3 + 1) = 0;
for j=1:N2,
    r(j,1) = tamc_radial_profile(j,1);
    r(j,2) = tamc_radial_profile(a*N2 + j,1);
    r(j,3) = tamc_radial_profile(b*N2 + j,1);
    r(j,4) = tamc_radial_profile(c*N2 + j,1);
    r(j,5) = tamc_radial_profile(d*N2 + j,1);
    r(j,6) = tamc_radial_profile(e*N2 + j,1);
    vel(j,1)=tamc_radial_profile(j,1);
    vel(j,2)=tamc_radial_profile(j,2);
    vel(j,3)=tamc_radial_profile(a*N2 + j,2);
    vel(j,4)=tamc_radial_profile(b*N2 + j,2);
    vel(j,5)=tamc_radial_profile(c*N2 + j,2);
    vel(j,6)=tamc_radial_profile(d*N2 + j,2);
    vel(j,7)=tamc_radial_profile(e*N2 + j,2);
    rho(j,1)=tamc_radial_profile(j,1);
    rho(j,2)=tamc_radial_profile(j,3);
    rho(j,3)=tamc_radial_profile(a*N2 + j,3);
    rho(j,4)=tamc_radial_profile(b*N2 + j,3);
    rho(j,5)=tamc_radial_profile(c*N2 + j,3);
    rho(j,6)=tamc_radial_profile(d*N2 + j,3);
    rho(j,7)=tamc_radial_profile(e*N2 + j,3);
    pressure(j,1)=tamc_radial_profile(j,1);
    pressure(j,2)=tamc_radial_profile(j,4);
    pressure(j,3)=tamc_radial_profile(a*N2 + j,4);
    pressure(j,4)=tamc_radial_profile(b*N2 + j,4);
    pressure(j,5)=tamc_radial_profile(c*N2 + j,4);
    pressure(j,6)=tamc_radial_profile(d*N2 + j,4);
    pressure(j,7)=tamc_radial_profile(e*N2 + j,4);       
    crpressure(j,1)=tamc_radial_profile(j,1);
    crpressure(j,2)=tamc_radial_profile(j,5);
    crpressure(j,3)=tamc_radial_profile(a*N2 + j,5);
    crpressure(j,4)=tamc_radial_profile(b*N2 + j,5);
    crpressure(j,5)=tamc_radial_profile(c*N2 + j,5);
    crpressure(j,6)=tamc_radial_profile(d*N2 + j,5);
    crpressure(j,7)=tamc_radial_profile(e*N2 + j,5);      
    temp(j,1)=tamc_radial_profile(j,1);
    temp(j,2)=tamc_radial_profile(j,6);
    temp(j,3)=tamc_radial_profile(a*N2 + j,6);
    temp(j,4)=tamc_radial_profile(b*N2 + j,6);
    temp(j,5)=tamc_radial_profile(c*N2 + j,6);
    temp(j,6)=tamc_radial_profile(d*N2 + j,6);
    temp(j,7)=tamc_radial_profile(e*N2 + j,6);
    flux(j,1)=tamc_radial_profile(j,1);
    flux(j,2)=tamc_radial_profile(j,2)*tamc_radial_profile(j,3)*flux(j,1)*flux(j,1);
    flux(j,3)=tamc_radial_profile(a*N2 + j,2)*tamc_radial_profile(a*N2 + j,3)*flux(j,1)*flux(j,1);
    flux(j,4)=tamc_radial_profile(b*N2 + j,2)*tamc_radial_profile(b*N2 + j,3)*flux(j,1)*flux(j,1);
    flux(j,5)=tamc_radial_profile(c*N2 + j,2)*tamc_radial_profile(c*N2 + j,3)*flux(j,1)*flux(j,1);
    flux(j,6)=tamc_radial_profile(d*N2 + j,2)*tamc_radial_profile(d*N2 + j,3)*flux(j,1)*flux(j,1);
    flux(j,7)=tamc_radial_profile(e*N2 + j,2)*tamc_radial_profile(e*N2 + j,3)*flux(j,1)*flux(j,1);  
end;
figure(1);
plot (r(1:N2,1),vel(1:N2,2),'cyan',r(1:N2,2),vel(1:N2,3),'green',r(1:N2,3),vel(1:N2,4),'blue',r(1:N2,4),vel(1:N2,5),'black',r(1:N2,5),vel(1:N2,6),'yellow',r(1:N2,6),vel(1:N2,7),'red');
title ('u(x)');
xlabel ('r cm');
ylabel ('U cm/s');
grid ;
figure(2);
plot (r(1:N2,1),rho(1:N2,2),'cyan',r(1:N2,2),rho(1:N2,3),'green',r(1:N2,3),rho(1:N2,4),'blue',r(1:N2,4),rho(1:N2,5),'black',r(1:N2,5),rho(1:N2,6),'yellow',r(1:N2,6),rho(1:N2,7),'red');
title ('density');
xlabel ('r cm');
ylabel ('rho 10^-5 g/cm^3');
grid ;
figure(3);
plot (r(1:N2,1),pressure(1:N2,2),'cyan',r(1:N2,2),pressure(1:N2,3),'green',r(1:N2,3),pressure(1:N2,4),'blue',r(1:N2,4),pressure(1:N2,5),'black',r(1:N2,5),pressure(1:N2,6),'yellow',r(1:N2,6),pressure(1:N2,7),'red');
title ('pressure');
xlabel ('r cm');
ylabel ('P din/cm^2');
grid ;
figure(4);
plot (r(1:N2,1),crpressure(1:N2,2),'cyan',r(1:N2,2),crpressure(1:N2,3),'green',r(1:N2,3),crpressure(1:N2,4),'blue',r(1:N2,4),crpressure(1:N2,5),'black',r(1:N2,5),crpressure(1:N2,6),'yellow',r(1:N2,6),crpressure(1:N2,7),'red');
title ('cosnic ray pressure');
xlabel ('r cm');
ylabel ('P din/cm^2');
grid;
figure(5);
plot (r(1:N2,1),temp(1:N2,2),'cyan',r(1:N2,2),temp(1:N2,3),'green',r(1:N2,3),temp(1:N2,4),'blue',r(1:N2,4),temp(1:N2,5),'black',r(1:N2,5),temp(1:N2,6),'yellow',r(1:N2,6),temp(1:N2,7),'red');
title ('tempertature');
xlabel ('r cm');
ylabel ('T');
grid;
%figure(6);
%plot (r(1:N2,1),flux(1:N2,2),'cyan',r(1:N2,2),flux(1:N2,3),'green',r(1:N2,3),flux(1:N2,4),'blue',r(1:N2,4),flux(1:N2,5),'black',r(1:N2,5),flux(1:N2,6),'yellow',r(1:N2,6),flux(1:N2,7),'red');
%title ('mass flux');
%xlabel ('r cm');
%ylabel ('g*cm/s');
%grid;