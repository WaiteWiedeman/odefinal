clc;clear;close all;

[t,y] = ode45(@vdp1,[0 50],[2; 0]);

figure;
plot(t,y(:,1),'-',t,y(:,2),'--')
title('Solution of van der Pol Equation (\mu = 1) with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

figure;
plot(y(:,1),y(:,2))
title('Phase plot of van der Pol Equation (\mu = 1) with ODE45');
xlabel('y1');
ylabel('y2');
