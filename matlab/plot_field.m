clear;
load field.dat;

N1=1;
N2=1000;
N3=100;


N4 = size(field,1);
a1 = fix(N4/(4*N3));
a2 = fix(N4/(2*N3));
a3 = fix(N4/N3) - 1;

funLocal(1:N3,1:4) = 0;
growthRate(1:N3,1:4) = 0;

for i = 1: N3,
    
    funLocal(i,1) = field(i,1);
    funLocal(i,2) = field(a1*N3 + i,2);
    funLocal(i,3) = field(a2*N3 + i,2);
    funLocal(i,4) = field(a3*N3 + i,2);
    
    growthRate(i,1) = field(i,1);
    growthRate(i,2) = field(a1*N3 + i,3);
    growthRate(i,3) = field(a2*N3 + i,3);
    growthRate(i,4) = field(a3*N3 + i,3);

end

figure(1);
plot (funLocal(1:N3,1), funLocal(1:N3,2),'blue', funLocal(1:N3,1), funLocal(1:N3,3),'green', funLocal(1:N3,1), funLocal(1:N3,4),'red');
title ('W');
xlabel ('k cm^-1');
ylabel ('W');
grid ;

figure(2);
plot (growthRate(1:N3,1), growthRate(1:N3,2),'blue', growthRate(1:N3,1), growthRate(1:N3,3),'green', growthRate(1:N3,1), growthRate(1:N3,4),'red');
title ('G');
xlabel ('k cm^-1');
ylabel ('G');
grid ;