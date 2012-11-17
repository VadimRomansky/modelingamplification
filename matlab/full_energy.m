clear;
load tamc_iteration.dat;
N=size(tamc_iteration,1);

figure(1);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,2),'blue',tamc_iteration(1:N,1),tamc_iteration(1:N,3),'red');
title ('energy');
xlabel ('iteration');
ylabel ('energy erg');
grid ;

figure(2);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,4),'blue',tamc_iteration(1:N,1),tamc_iteration(1:N,5),'red');
title ('momentumZ');
xlabel ('iteration');
ylabel ('momentum g*cm/s');
grid ;

figure(5);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,6),'blue');
title ('particles count');
xlabel ('iteration');
ylabel ('count');
grid ;

figure(6);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,7),'blue');
title ('particles weight');
xlabel ('iteration');
ylabel ('weight');
grid ;

figure(7);
plot (tamc_iteration(2:N,1),tamc_iteration(2:N,8),'blue');
title ('shock wave point');
xlabel ('iteration');
ylabel ('r cm');
grid ;

figure(8);
plot (tamc_iteration(1:N,1),tamc_iteration(1:N,8),'blue', tamc_iteration(1:N,1),tamc_iteration(N,8)*exp((2/5)*log(tamc_iteration(1:N,1)/tamc_iteration(N,1))),'red');
title ('shock wave point');
xlabel ('iteration');
ylabel ('r cm');
grid ;

load EM_cur.dat;
load EM_in.dat;
load EM_out.dat;
load M_wall.dat;
N=size(EM_cur,1);

figure(9);
plot (EM_cur(1:N,1),EM_cur(1:N,2),'blue',EM_cur(1:N,1),EM_in(1:N,2) + EM_out(1:N,2),'red');
title ('energy');
xlabel ('iteration');
ylabel ('energy erg');
grid ;

figure(10);
plot (EM_cur(1:N,1),EM_cur(1:N,3),'blue',EM_cur(1:N,1),EM_in(1:N,3) + EM_out(1:N,3) + M_wall(1:N,2),'red');
title ('momentumZ');
xlabel ('iteration');
ylabel ('momentum g*cm/s');
grid ;

clear;