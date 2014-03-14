clear;
load grid.dat;
N1=1;
N2=96;
vel(1:96,1:5)=0;
Fm(1:96,1:5)=0;
Fen(1:96,1:5)=0;
jd(1:96,1:5)=0;
for j=1:87,
    vel(j,1)=grid(j,4);
    vel(j,2)=grid(j,5);
    vel(j,3)=grid(96+j,5);
    vel(j,4)=grid(192+j,5);
    vel(j,5)=grid(288+j,5);
    Fm(j,2)=grid(j,6);
    Fm(j,3)=grid(96+j,6);
    Fm(j,4)=grid(192+j,6);
    Fm(j,5)=grid(288+j,6);
    Fen(j,2)=grid(j,7);
    Fen(j,3)=grid(96+j,7);
    Fen(j,4)=grid(192+j,7);
    Fen(j,5)=grid(288+j,7);
    jd(j,2)=grid(j,17);
    jd(j,3)=grid(96+j,17);
    jd(j,4)=grid(192+j,17);
    jd(j,5)=grid(288+j,17);
end;
for j=88:96,
    vel(j,1)=4+grid(j,3);
    vel(j,2)=grid(j,5);
    vel(j,3)=grid(96+j,5);
    vel(j,4)=grid(192+j,5);
    vel(j,5)=grid(288+j,5);
    Fm(j,2)=grid(j,6);
    Fm(j,3)=grid(96+j,6);
    Fm(j,4)=grid(192+j,6);
    Fm(j,5)=grid(288+j,6);
    Fen(j,2)=grid(j,7);
    Fen(j,3)=grid(96+j,7);
    Fen(j,4)=grid(192+j,7);
    Fen(j,5)=grid(288+j,7);
    jd(j,2)=grid(j,17);
    jd(j,3)=grid(96+j,17);
    jd(j,4)=grid(192+j,17);
    jd(j,5)=grid(288+j,17);
end;
figure(1);
plot (vel(1:96,1),vel(1:96,5),'blue');
title ('u(x)');
xlabel ('');
ylabel ('');
grid ;
figure(2);
plot (vel(1:96,1),Fm(1:96,5),'blue');
title ('Fm');
xlabel ('');
ylabel ('');
grid ;
figure(3);
plot (vel(1:96,1),Fen(1:96,5),'blue');
title ('Fen');
xlabel ('');
ylabel ('');
grid ;
figure(4);
plot (vel(1:96,1),jd(1:96,5),'blue');
title ('jd');
xlabel ('');
ylabel ('');
grid ;