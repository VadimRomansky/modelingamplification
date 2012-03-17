clear;
load zpdf1.dat;
load zpdf2.dat;
load zpdf3.dat;
load zpdf4.dat;
load zpdf0.dat;
N=500;

figure(1);
plot (zpdf0(1:N,1),zpdf0(1:N,2),'blue',zpdf0(1:N,1),zpdf0(1:N,3),'red');
title ('pdf');
xlabel ('pz');
ylabel ('dN/N');
grid ;

figure(2);
plot (zpdf1(1:N,1),zpdf1(1:N,2),'blue',zpdf1(1:N,1),zpdf1(1:N,3),'red');
title ('pdf');
xlabel ('pz');
ylabel ('dN/N');
grid ;

figure(3);
plot (zpdf2(1:N,1),zpdf2(1:N,2),'blue',zpdf2(1:N,1),zpdf2(1:N,3),'red');
title ('pdf');
xlabel ('pz');
ylabel ('dN/N');
grid ;

figure(4);
plot (zpdf3(1:N,1),zpdf3(1:N,2),'blue',zpdf3(1:N,1),zpdf3(1:N,3),'red');
title ('pdf');
xlabel ('pz');
ylabel ('dN/N');
grid ;

figure(5);
plot (zpdf4(1:N,1),zpdf4(1:N,2),'blue',zpdf4(1:N,1),zpdf4(1:N,3),'red');
title ('pdf');
xlabel ('pz');
ylabel ('dN/N');
grid ;