clear;
load zpdf1.dat;
load zpdf2.dat;
load zpdf3.dat;
load zpdf4.dat;
load zpdf5.dat;
load zpdf6.dat;
load zpdf7.dat;
load zpdf8.dat;
load zpdf9.dat;
load zpdf0.dat;
N=size(zpdf1,1);

figure(1);
plot (zpdf0(1:N,1),zpdf0(1:N,2),'blue',zpdf0(1:N,1),zpdf0(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(2);
plot (zpdf1(1:N,1),zpdf1(1:N,2),'blue',zpdf1(1:N,1),zpdf1(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(3);
plot (zpdf2(1:N,1),zpdf2(1:N,2),'blue',zpdf2(1:N,1),zpdf2(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20');
ylabel ('dN/N');
grid ;

figure(4);
plot (zpdf3(1:N,1),zpdf3(1:N,2),'blue',zpdf3(1:N,1),zpdf3(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(5);
plot (zpdf4(1:N,1),zpdf4(1:N,2),'blue',zpdf4(1:N,1),zpdf4(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(6);
plot (zpdf5(1:N,1),zpdf5(1:N,2),'blue',zpdf5(1:N,1),zpdf5(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(7);
plot (zpdf6(1:N,1),zpdf6(1:N,2),'blue',zpdf6(1:N,1),zpdf6(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(8);
plot (zpdf7(1:N,1),zpdf7(1:N,2),'blue',zpdf7(1:N,1),zpdf7(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(9);
plot (zpdf8(1:N,1),zpdf8(1:N,2),'blue',zpdf8(1:N,1),zpdf8(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;

figure(10);
plot (zpdf9(1:N,1),zpdf9(1:N,2),'blue',zpdf9(1:N,1),zpdf9(1:N,3),'red');
title ('pdf');
xlabel ('pz 10^-20 g*cm/s');
ylabel ('dN/N');
grid ;