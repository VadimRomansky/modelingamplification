clear;
load traectory.dat;
N1=1;
N2=size(traectory,1);


figure(1);
plot (traectory(1:N2,1),traectory(1:N2,2),'red',traectory(1:N2,1),traectory(1:N2,3),'green',traectory(1:N2,1),traectory(1:N2,4),'blue');
title ('r');
xlabel ('t');
ylabel (')');
legend(4, 'x', 'y','z');
grid ;

figure(2);
plot (traectory(1:N2,1),traectory(1:N2,5),'red',traectory(1:N2,1),traectory(1:N2,6),'green',traectory(1:N2,1),traectory(1:N2,7),'blue', traectory(1:N2,1),traectory(1:N2,8),'black');
title ('p');
xlabel ('t');
ylabel ('p');
legend(4, 'px', 'py','pz', 'p');
grid ;

figure(3);
plot (traectory(1:N2,3),traectory(1:N2,4));
title ('r');
xlabel ('y');
ylabel ('z');
grid ;
