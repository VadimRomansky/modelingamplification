clear;
load field.dat;
load full_field.dat;
load diff_coef.dat;
load kfile.dat;
load xfile.dat;

N1=1;
N2=size(xfile,1);
N3=size(kfile,1);
Nmomentum = 100;


N4 = size(field,1);
a1 = fix(N4/(4*N3));
a2 = fix(N4/(2*N3));
a3 = fix(N4/N3) - 1;

funLocal(1:N3,1:4) = 0;
fun(1:N3,1:4) = 0;
growthRate(1:N3,1:4) = 0;
growthRate2(1:N3,1:4) = 0;
coef(1:Nmomentum,1:10) = 0;

for i = 1 : Nmomentum
    coef(i,1) = diff_coef(i,1);
    coef(i,2) = diff_coef(a1*Nmomentum + i,2);
    coef(i,3) = diff_coef(a2*Nmomentum + i,2);
    coef(i,4) = diff_coef(a3*Nmomentum + i,2);
    coef(i,5) = diff_coef(a1*Nmomentum + i,3);
    coef(i,6) = diff_coef(a2*Nmomentum + i,3);
    coef(i,7) = diff_coef(a3*Nmomentum + i,3);
    coef(i,8) = diff_coef(a1*Nmomentum + i,4);
    coef(i,9) = diff_coef(a2*Nmomentum + i,4);
    coef(i,10) = diff_coef(a3*Nmomentum + i,4);
end;

for i = 1: N3,
    
    funLocal(i,1) = field(i,1);
    funLocal(i,2) = field(a1*N3 + i,2);
    funLocal(i,3) = field(a2*N3 + i,2);
    funLocal(i,4) = field(a3*N3 + i,2);
    
    fun(i,1) = field(i,1);
    fun(i,2) = field(a1*N3 + i,4);
    fun(i,3) = field(a2*N3 + i,4);
    fun(i,4) = field(a3*N3 + i,4);
    
    growthRate(i,1) = field(i,1);
    growthRate(i,2) = field(a1*N3 + i,3);
    growthRate(i,3) = field(a2*N3 + i,3);
    growthRate(i,4) = field(a3*N3 + i,3);

    growthRate2(i,1) = field(i,1);
    growthRate2(i,2) = field(a1*N3 + i,3)*funLocal(i,2);
    growthRate2(i,3) = field(a2*N3 + i,3)*funLocal(i,3);
    growthRate2(i,4) = field(a3*N3 + i,3)*funLocal(i,4);
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

figure(3);
plot (fun(1:N3,1), fun(1:N3,2),'blue', fun(1:N3,1), fun(1:N3,3),'green', fun(1:N3,1), fun(1:N3,4),'red');
title ('W');
xlabel ('k cm^-1');
ylabel ('W');
grid ;

figure(4);
plot (growthRate2(1:N3,1), growthRate2(1:N3,2),'blue', growthRate2(1:N3,1), growthRate2(1:N3,3),'green', growthRate2(1:N3,1), growthRate2(1:N3,4),'red');
title ('G');
xlabel ('k cm^-1');
ylabel ('G');
grid ;

figure(5);
plot (coef(1:Nmomentum,1), coef(1:Nmomentum,2),'blue',coef(1:Nmomentum,1), coef(1:Nmomentum,5),'blue',coef(1:Nmomentum,1), coef(1:Nmomentum,8),'blue', coef(1:Nmomentum,1), coef(1:Nmomentum,3),'green', coef(1:Nmomentum,1), coef(1:Nmomentum,4),'red', coef(1:Nmomentum,1), coef(1:Nmomentum,7),'red', coef(1:Nmomentum,1), coef(1:Nmomentum,10),'red');
title ('D');
xlabel ('p g*cm/s');
ylabel ('D cm^2/s');
grid ;

figure(6);
[X Y] = meshgrid(kfile, xfile);
mesh(X, Y, full_field);
grid;

%xgrid(1:N3)=0;
%for i = 1:N2
   % xgrid(i) = -xfile(i);
%end;

%figure(6);
%[X Y] = meshgrid(kfile, xgrid);
%mesh(X, Y, full_field);
%grid;