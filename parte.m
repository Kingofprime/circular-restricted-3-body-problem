%% This code will take 2-3 minutes to run
clc;
clear;
clear all;
load 'EM_L4-304P1.mat'
t = 0:0.001:T;
omega =1;
%% ode45
options=odeset('reltol',1e-15,'abstol',1e-15);
[t,xnp]= ode45(@(t,xnp) out(xnp,MU1),t,x0,options);
[t,xp]=ode45(@(t,xp) out(xp,MU1),t,x0+perturbation,options); %perturb
 
%% rotation about 3rd axis
for i=1:length(t)
    theta= omega*t(i);
    dcm3=[cos(theta),-sin(theta);sin(theta),cos(theta)];
    xin(i,:)=dcm3*xnp(i,1:2)';
    xinp(i,:)=dcm3*xp(i,1:2)'; %perturb
    deltax(i,:)=xp(i,:)-xnp(i,:); 
end
%% Orbit graphs
figure 
plot(xnp(:,1),xnp(:,2)); % xnt vs ynt
xlabel("x position");
ylabel("y position");
title('Nominal Orbit');

figure
plot (xin(:,1),xin(:,2)) %Xnt vs Ynt
xlabel("x position");
ylabel("y position");
title('Inertial Frame Orbit');
%% non -linear graphs
figure
subplot(2,1,1)
plot(t,deltax(:,1))
hold on
plot(t,deltax(:,2))
hold off
xlabel("time");
ylabel("position");
title('Non-linear Position vs Time');
legend('X','Y');

subplot(2,1,2)
plot(t,deltax(:,3))
hold on
plot(t,deltax(:,4))
xlabel("time");
ylabel("velocity");
title('Non-linear velocity vs Time');
legend('X','Y');
hold off

%% linearization , similar to circularhoop.m code
syms x y mu
uxx=(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)) + 1; %from derivates.m
uyy=(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1; %from derivates.m
uxy=(3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2)^(5/2); %from derivates.m

a = [0, 0, 1, 0; 0, 0, 0, 1;uxx, uxy, 0, 2;uxy, uyy, -2, 0];
na = double(subs(a,[x y mu],[0.5-MU1 sin(deg2rad(60)) MU1]));

[vec,val]=eig(na);
val=diag(val);
coeff=vec\(perturbation); % matrix left division
for i=1:length(t)
    deltalinear(i,:)=(coeff(1)*vec(:,1)*exp(val(1)*t(i)))+(coeff(2)*vec(:,2)*exp(val(2)*t(i)))+(coeff(3)*vec(:,3)*exp(val(3)*t(i)))+(coeff(4)*vec(:,4)*exp(val(4)*t(i)));
    if imag(deltalinear(i,:))<=(1e-10)
        deltalinear(i,:)=real(deltalinear(i,:)); % to convert imag to real
    end
end
%% linearized graph
figure
subplot(2,1,1)
plot(t,(deltalinear(:,1)))
hold on
plot(t,(deltalinear(:,2)))
xlabel("time");
ylabel("position");
title('linear position vs Time');
legend('X','Y');
hold off

subplot(2,1,2)
plot(t,(deltalinear(:,3)))
hold on
plot(t,(deltalinear(:,4)))
xlabel("time");
ylabel("velocity");
title('linear velocity vs Time');
legend('X','Y');
hold off

%% comparing linear to non-linear graphs
figure
subplot(2,1,1)
plot(t,deltax(:,1))
hold on
plot(t,deltalinear(:,1))
hold off
xlabel("time");
ylabel("position");
title('non-linear and linear vs Time for X');
legend('Non-Linear','Linear')

subplot(2,1,2)
plot(t,deltax(:,2))
hold on
plot(t,deltalinear(:,2))
hold off
xlabel("time");
ylabel("position");
title('non-linear and linear vs Time for Y');
legend('Non-Linear','Linear')

figure
subplot(2,1,1)
plot(t,deltax(:,3))
hold on
plot(t,deltalinear(:,3))
hold off
xlabel("time");
ylabel("velocity");
title('non-linear and linear vs Time for X');
legend('Non-Linear','Linear')

subplot(2,1,2)
plot(t,deltax(:,4))
hold on
plot(t,deltalinear(:,4))
hold off
xlabel("time");
ylabel("velocity");
title('non-linear and linear vs Time for Y');
legend('Non-Linear','Linear')

%% Intial ode45 conditions 
function out = out(x0,mu)
x=x0(1);
y=x0(2);
x1=x0(3);
y1=x0(4);
ux = -((1-mu)*(x+mu)/((((x+mu)^2)+(y^2))^0.5)^3)-(mu*(x-1+mu)/((((x-1+mu)^2)+(y^2))^0.5)^3)+x; % from derivatives.m
uy = -((1-mu)*y/((((x+mu)^2)+(y^2))^0.5)^3) - (mu*y/((((x-1+mu)^2)+(y^2))^0.5)^3) + y; %derivatives.m
x2 = ux + (2*y1);
y2 = uy - (2*x1);
out =[x1;y1;x2;y2];
end