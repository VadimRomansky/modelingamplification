clear;
load shock_wave.dat;
N1=1;
N2=size(shock_wave,1);
year = 3.15*10^7;
t(1:N2,1:2) = 0;
for i = 1:N2
    t(i,1) = shock_wave(i,2)/year;
    t(i,2) = shock_wave(N2,4)*(shock_wave(i,2)/shock_wave(N2,2))^(2/5);
end
figure(1);
plot (t(1:N2,1),shock_wave(1:N2,4),'blue', t(1:N2,1), t(1:N2,2), 'red');
title ('shock wawe r');
xlabel ('time year');
ylabel ('r');
legend('experiment','sedov',4);
grid ;