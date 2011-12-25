clear;
load tamc_radial_profile.dat;
N1=1;
N2=100;
N3=6;
a = 50;
b = 100;
c = 100;
d = 100;
e = 100;
vel(1:N2,1:N3 + 1)=0;
rho(1:N2,1:N3 + 1) =0;
temp(1:N2,1:N3 + 1) = 0;
avervel(1:N2,1:N3 + 1)=0;
for j=1:N2,
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
    temp(j,1)=tamc_radial_profile(j,1);
    temp(j,2)=tamc_radial_profile(j,4);
    temp(j,3)=tamc_radial_profile(a*N2 + j,4);
    temp(j,4)=tamc_radial_profile(b*N2 + j,4);
    temp(j,5)=tamc_radial_profile(c*N2 + j,4);
    temp(j,6)=tamc_radial_profile(e*N2 + j,4);
    temp(j,7)=tamc_radial_profile(d*N2 + j,4);
    avervel(j,1)=tamc_radial_profile(j,1);
    avervel(j,2)=tamc_radial_profile(j,5);
    avervel(j,3)=tamc_radial_profile(a*N2 + j,5);
    avervel(j,4)=tamc_radial_profile(b*N2 + j,5);
    avervel(j,5)=tamc_radial_profile(c*N2 + j,5);
    avervel(j,6)=tamc_radial_profile(d*N2 + j,5);
    avervel(j,7)=tamc_radial_profile(e*N2 + j,5);    
end;
figure(1);
plot (vel(1:N2,1),vel(1:N2,2),'red',vel(1:N2,1),vel(1:N2,3),'green',vel(1:N2,1),vel(1:N2,4),'blue',vel(1:N2,1),vel(1:N2,5),'black',vel(1:N2,1),vel(1:N2,6),'yellow',vel(1:N2,1),vel(1:N2,7),'cyan');
title ('mass u(x)');
xlabel ('r');
ylabel ('U');
grid ;
figure(2);
plot (rho(1:N2,1),rho(1:N2,2),'red',rho(1:N2,1),rho(1:N2,3),'green',rho(1:N2,1),rho(1:N2,4),'blue',rho(1:N2,1),rho(1:N2,5),'black',rho(1:N2,1),rho(1:N2,6),'yellow',rho(1:N2,1),rho(1:N2,7),'cyan');
title ('density');
xlabel ('r');
ylabel ('rho');
grid ;
figure(3);
plot (temp(1:N2,1),temp(1:N2,2),'red',temp(1:N2,1),temp(1:N2,3),'green',temp(1:N2,1),temp(1:N2,4),'blue',temp(1:N2,1),temp(1:N2,5),'black',temp(1:N2,1),temp(1:N2,6),'yellow',temp(1:N2,1),temp(1:N2,7),'cyan');
title ('temperature');
xlabel ('r');
ylabel ('T');
grid ;
figure(4);
plot (avervel(1:N2,1),avervel(1:N2,2),'red',avervel(1:N2,1),avervel(1:N2,3),'green',avervel(1:N2,1),avervel(1:N2,4),'blue',avervel(1:N2,1),avervel(1:N2,5),'black',avervel(1:N2,1),avervel(1:N2,6),'yellow',avervel(1:N2,1),avervel(1:N2,7),'cyan');
title ('average u(x)');
xlabel ('r');
ylabel ('U');
grid ;