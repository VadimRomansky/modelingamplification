clear;
load tamc_energy_pdf0.dat;
load tamc_energy_pdf1.dat;
load tamc_energy_pdf2.dat;
load tamc_energy_pdf3.dat;
load tamc_energy_pdf4.dat;
N=size(tamc_energy_pdf0,1);

figure(1);
plot (tamc_energy_pdf0(1:N,1),tamc_energy_pdf0(1:N,2),'blue',tamc_energy_pdf0(1:N,1),tamc_energy_pdf0(1:N,3),'red');
title ('pdf');
xlabel ('E');
ylabel ('dN/N');
grid ;

figure(2);
plot (tamc_energy_pdf1(1:N,1),tamc_energy_pdf1(1:N,2),'blue',tamc_energy_pdf1(1:N,1),tamc_energy_pdf1(1:N,3),'red');
title ('pdf');
xlabel ('E');
ylabel ('dN/N');
grid ;

figure(3);
plot (tamc_energy_pdf2(1:N,1),tamc_energy_pdf2(1:N,2),'blue',tamc_energy_pdf2(1:N,1),tamc_energy_pdf2(1:N,3),'red');
title ('pdf');
xlabel ('E');
ylabel ('dN/N');
grid ;

figure(4);
plot (tamc_energy_pdf3(1:N,1),tamc_energy_pdf3(1:N,2),'blue',tamc_energy_pdf3(1:N,1),tamc_energy_pdf3(1:N,3),'red');
title ('pdf');
xlabel ('E');
ylabel ('dN/N');
grid ;

figure(5);
plot (tamc_energy_pdf4(1:N,1),tamc_energy_pdf4(1:N,2),'blue',tamc_energy_pdf4(1:N,1),tamc_energy_pdf4(1:N,3),'red');
title ('pdf');
xlabel ('E');
ylabel ('dN/N');
grid ;