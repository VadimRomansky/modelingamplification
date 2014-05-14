clear;
load tamc_radial_profile.dat;
N1=1;
N2=1000;
e = size(tamc_radial_profile,1)/N2 - 1;

r(1:N2) = 0;
vel(1:N2) = 0;
rho(1:N2,1:N3 + 1)= 0;
temp(1:N2) = 0;
pressure(1:N2) = 0;

maxV = 0;
maxPressure = 0;
maxDensity = 0;

for j=1:N2
    if(tamc_radial_profile(e*N2 + j,2) > maxV)
        maxV = tamc_radial_profile(e*N2 + j,2);
    end;
    if(tamc_radial_profile(e*N2 + j,3) > maxDensity)
        maxDnsity = tamc_radial_profile(e*N2 + j,3);
    end;
    if(tamc_radial_profile(e*N2 + j,4) > maxPressure)
        maxPressure = tamc_radial_profile(e*N2 + j,4);
    end;
end;

for j=1:N2,
    r(j) = tamc_radial_profile(e*N2 + j,1);
    vel(j)=tamc_radial_profile(e*N2 + j,2)/maxV;
    rho(j)=tamc_radial_profile(e*N2 + j,3)/maxDensity;
    pressure(j)=tamc_radial_profile(e*N2 + j,4)/maxPressure;
end;
figure(1);
plot (r(1:N2),vel(1:N2),'red',r(1:N2),rho(1:N2),'green',r(1:N2),pressure(1:N2),'blue');
title ('u, rho, p');
xlabel ('r cm');
ylabel ('U  rho P');
legend('скорость','плотность','давление',4)
grid ;
