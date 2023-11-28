clc; clear; close all;

syms x0(t) x1(t) x2(t) t a b
alph = 0.01;

x0(t) = a*cos(t) + b*sin(t);

x1(t) = ((-1/32)*b^3 + (1/32)*a^2*b + (1/8)*b)*cos(t) - (1/32)*a*(5*b^2 - 7*a^2 + 16)*sin(t) ...
    + t*((1/2)*a*(1-(a^2+b^2)/4)*cos(t)) + t*((1/2)*b*(1 - (b^2 - a^2)/4)*sin(t))...
    + (1/8)*b*((b^2 - a^2)/4 - 1)*cos(3*t) + (1/8)*a*((3*b^2 - a^2)/4)*sin(3*t);

x_1T = x0 + alph*x1;
x_1T = subs(x_1T, [a b], [1.5 sqrt(4-1.5^2)]);
v_1T = diff(x_1T,t);

figure;
fplot(x_1T,[0 400])
hold on 
fplot(v_1T,[0 400])
hold off
title('Solution of van der Pol Equation (\mu = 0.01) using Perturbation approximation');
xlabel('t');
ylabel('Solution');
legend('x','v')

figure;
fplot(x_1T,v_1T)
title('Phase plot of van der Pol Equation (\mu = 0.01) with Perturbation approximation');
xlabel('x');
ylabel('v');

eqn = diff(x2(t),t,t) + x2(t) == diff(x1(t),t) - x0(t)^2*diff(x1(t),t) - 2*x0(t)*diff(x0(t),t)*x1(t);
Dx2 = diff(x2,t);
cond = [x2(0)==a, Dx2(0)==b];

S_x2 = dsolve(eqn,cond);
S_2T = x0 + alph*x1 + alph^2*S_x2;
x_2T = subs(S_2T, [a b], [1.5 sqrt(4-1.5^2)]);
v_2T = diff(x_2T,t);

figure;
fplot(x_2T,[0 50])
hold on 
fplot(v_2T,[0 50])
hold off
title('Solution of van der Pol Equation (\mu = 0.01) using Perturbation approximation');
xlabel('t');
ylabel('Solution');
legend('x','v')

figure;
fplot(x_2T,v_2T)
title('Phase plot of van der Pol Equation (\mu = 0.01) with Perturbation approximation');
xlabel('x');
ylabel('v');

%% 
clc;clear;close all;

[t,y] = ode45(@vdp1,[0 50],[1.5 sqrt(4-1.5^2)]);

figure;
plot(t,y(:,1),'-',t,y(:,2),'--')
title('Solution of van der Pol Equation (\mu = 0.01) with ODE45');
xlabel('Time t');
ylabel('Solution');
legend('x','v')

figure;
plot(y(:,1),y(:,2))
title('Phase plot of van der Pol Equation (\mu = 0.01) with ODE45');
xlabel('x');
ylabel('v');

function dydt = vdp1(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

dydt = [y(2); 0.01*(1-y(1)^2)*y(2)-y(1)];
end

