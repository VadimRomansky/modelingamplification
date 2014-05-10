clear;
load extra_iterations.dat;
N = size(extra_iterations,1);
figure(1);
plot (extra_iterations(1:N,1),extra_iterations(1:N,3),'red');
title ('mass');
xlabel ('ieration');
ylabel ('m g');
grid ;
figure(2);
plot (extra_iterations(1:N,1),extra_iterations(1:N,4),'red');
title ('momentum');
xlabel ('ieration');
ylabel ('p g*cm/s');
grid ;
figure(3);
plot (extra_iterations(1:N,1),extra_iterations(1:N,5),'red', extra_iterations(1:N,1),extra_iterations(1:N,6),'black', extra_iterations(1:N,1),extra_iterations(1:N,7),'green', extra_iterations(1:N,1),extra_iterations(1:N,8),'yellow', extra_iterations(1:N,1),extra_iterations(1:N,9),'blue');
title ('energy');
xlabel ('ieration');
ylabel ('e erg');
legend('full energy', 'kinetic', 'temal', 'particles', 'magnetic',4);
grid ;