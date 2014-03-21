clear;
load fp4.dat;
load Np.dat;
figure(1);
plot (fp4(1:(Np-1),2),fp4(1:(Np-1),3),'g');
title ('f');
xlabel ('p');
ylabel ('f');
grid ;