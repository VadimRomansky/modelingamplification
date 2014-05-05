clear;
load shock_wave.dat;
N1=1;
N2=size(shock_wave,1);
year = 3.15*10^7;
t(1:N2,1:7) = 0;
alpha = 1;
E = alpha*10^49;
rho = 1.6*10^-24;
for i = 1:N2
    t(i,1) = shock_wave(i,2)/year;
    t(i,2) = shock_wave(N2,4)*(shock_wave(i,2)/shock_wave(N2,2))^(2/5);
    t(i,3) = (E*shock_wave(i,2)*shock_wave(i,2)/rho)^(1/5);
    t(i,4) = shock_wave(i,5);
    t(i,5) = shock_wave(i,6);
    t(i,6) = shock_wave(N2,5)*(shock_wave(i,2)/shock_wave(N2,2))^(-3/5);
    t(i,7) = 0.4*((E/rho)^(1/5))*(shock_wave(i,2)^(-3/5));
end
t(1,6)=0;
figure(1);
plot (t(1:N2,1),shock_wave(1:N2,4),'blue', t(1:N2,1), t(1:N2,2), 'red', t(1:N2,1), t(1:N2,3), 'green');
title ('shock wawe r');
xlabel ('time year');
ylabel ('r');
legend('experiment','sedov','real sedov',4);
grid ;

figure(2);
plot (t(1:N2,1),t(1:N2,4),'blue', t(1:N2,1), t(1:N2,5), 'yellow', t(1:N2,1), t(1:N2,6), 'red', t(1:N2,1), t(1:N2,7), 'green');
title ('shock wawe V');
xlabel ('time year');
ylabel ('V');
legend('shockWave','gas','-3/5', 'sedov', 4);
grid ;