clear;
load tamc_radial_profile.dat;
N1=1;
N2=250;
N3=6;
e = size(tamc_radial_profile,1)/N2 - 1;
d = fix(3*e/4);
c = fix(e/2);
b = fix(e/4);
a = 5;
vel(1:N2,1:N3 + 1)=0;
rho(1:N2,1:N3 + 1) =0;
temp(1:N2,1:N3 + 1) = 0;
avervel(1:N2,1:N3 + 1)=0;
massvel(1:N2,1:N3 + 1) = 0;
crflux(1:N2,1:N3 + 1) = 0;
for j=1:N2,
    vel(j,1)=tamc_radial_profile(j,1);
    vel(j,2)=tamc_radial_profile(j,2);
    vel(j,3)=tamc_radial_profile(a*N2 + j,2);
    vel(j,4)=tamc_radial_profile(b*N2 + j,2);
    vel(j,5)=tamc_radial_profile(c*N2 + j,2);
    vel(j,6)=tamc_radial_profile(d*N2 + j,2);
    vel(j,7)=tamc_radial_profile(e*N2 + j,2);
    avervel(j,1)=tamc_radial_profile(j,1);
    avervel(j,2)=tamc_radial_profile(j,3);
    avervel(j,3)=tamc_radial_profile(a*N2 + j,3);
    avervel(j,4)=tamc_radial_profile(b*N2 + j,3);
    avervel(j,5)=tamc_radial_profile(c*N2 + j,3);
    avervel(j,6)=tamc_radial_profile(d*N2 + j,3);
    avervel(j,7)=tamc_radial_profile(e*N2 + j,3);  
    massvel(j,1)=tamc_radial_profile(j,1);   
    rho(j,1)=tamc_radial_profile(j,1);
    rho(j,2)=tamc_radial_profile(j,4);
    rho(j,3)=tamc_radial_profile(a*N2 + j,4);
    rho(j,4)=tamc_radial_profile(b*N2 + j,4);
    rho(j,5)=tamc_radial_profile(c*N2 + j,4);
    rho(j,6)=tamc_radial_profile(d*N2 + j,4);
    rho(j,7)=tamc_radial_profile(e*N2 + j,4);
    temp(j,1)=tamc_radial_profile(j,1);
    temp(j,2)=tamc_radial_profile(j,5);
    temp(j,3)=tamc_radial_profile(a*N2 + j,5);
    temp(j,4)=tamc_radial_profile(b*N2 + j,5);
    temp(j,5)=tamc_radial_profile(c*N2 + j,5);
    temp(j,6)=tamc_radial_profile(d*N2 + j,5);
    temp(j,7)=tamc_radial_profile(e*N2 + j,5); 
    crflux(j,1)=tamc_radial_profile(j,1);
    crflux(j,2)=tamc_radial_profile(N2 + j,6);
    crflux(j,3)=tamc_radial_profile(a*N2 + j,6);
    crflux(j,4)=tamc_radial_profile(b*N2 + j,6);
    crflux(j,5)=tamc_radial_profile(c*N2 + j,6);
    crflux(j,6)=tamc_radial_profile(d*N2 + j,6);
    crflux(j,7)=tamc_radial_profile(e*N2 + j,6);     
end;
figure(1);
plot (vel(1:N2,1),vel(1:N2,2),'cyan',vel(1:N2,1),vel(1:N2,3),'green',vel(1:N2,1),vel(1:N2,4),'blue',vel(1:N2,1),vel(1:N2,5),'black',vel(1:N2,1),vel(1:N2,6),'yellow',vel(1:N2,1),vel(1:N2,7),'red');
title ('u(x)');
xlabel ('r cm');
ylabel ('U cm/s');
grid ;
figure(2);
plot (rho(1:N2,1),rho(1:N2,2),'cyan',rho(1:N2,1),rho(1:N2,3),'green',rho(1:N2,1),rho(1:N2,4),'blue',rho(1:N2,1),rho(1:N2,5),'black',rho(1:N2,1),rho(1:N2,6),'yellow',rho(1:N2,1),rho(1:N2,7),'red');
title ('density');
xlabel ('r cm');
ylabel ('rho 10^-5 g/cm^3');
grid ;
figure(3);
plot (avervel(1:N2,1),avervel(1:N2,2),'cyan',avervel(1:N2,1),avervel(1:N2,3),'green',avervel(1:N2,1),avervel(1:N2,4),'blue',avervel(1:N2,1),avervel(1:N2,5),'black',avervel(1:N2,1),avervel(1:N2,6),'yellow',avervel(1:N2,1),avervel(1:N2,7),'red');
title ('average u(x)');
xlabel ('r cm');
ylabel ('U cm/s');
grid ;
figure(4);
plot (temp(1:N2,1),temp(1:N2,2),'cyan',temp(1:N2,1),temp(1:N2,3),'green',temp(1:N2,1),temp(1:N2,4),'blue',temp(1:N2,1),temp(1:N2,5),'black',temp(1:N2,1),temp(1:N2,6),'yellow',temp(1:N2,1),temp(1:N2,7),'red');
title ('tempertature');
xlabel ('r cm');
ylabel ('T');
grid ;
figure(5);
plot (crflux(1:N2,1),crflux(1:N2,2),'cyan',crflux(1:N2,1),crflux(1:N2,3),'green',crflux(1:N2,1),crflux(1:N2,4),'blue',crflux(1:N2,1),crflux(1:N2,5),'black',crflux(1:N2,1),crflux(1:N2,6),'yellow',crflux(1:N2,1),crflux(1:N2,7),'red');
title ('cosmic ray flux');
xlabel ('r cm');
ylabel ('j A/cm^2');
grid ;