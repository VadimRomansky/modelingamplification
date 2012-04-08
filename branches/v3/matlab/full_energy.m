clear;
load tamc_iteration.dat;
N=70;

figure(1);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,2),'blue',tamc_iteration(1:N,1),tamc_iteration(1:N,3),'red');
title ('energy');
xlabel ('iteration');
ylabel ('energy');
grid ;

figure(2);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,4),'blue',tamc_iteration(1:N,1),tamc_iteration(1:N,5),'red');
title ('momentumZ');
xlabel ('iteration');
ylabel ('momentum');
grid ;

figure(3);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,6),'blue',tamc_iteration(1:N,1),tamc_iteration(1:N,7),'red');
title ('momentumX');
xlabel ('iteration');
ylabel ('momentum');
grid ;

figure(4);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,8),'blue',tamc_iteration(1:N,1),tamc_iteration(1:N,9),'red');
title ('momentumY');
xlabel ('iteration');
ylabel ('momentum');
grid ;

figure(5);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,10),'blue');
title ('particles count');
xlabel ('iteration');
ylabel ('count');
grid ;

figure(6);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,11),'blue');
title ('particles weight');
xlabel ('iteration');
ylabel ('weight');
grid ;