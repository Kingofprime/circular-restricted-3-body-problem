%304 Project 1 - Groupmembers: Thomas Waltz, Ankit Gupta, and Kameron
%Metcalf

clc
syms x mu y;
p1 = (((x+mu)^2)+(y^2))^0.5; %from part a
p2 = (((x-1+mu)^2)+(y^2))^0.5; %from part a

ux = -((1-mu)*(x+mu)/p1^3)-(mu*(x-1+mu)/p2^3)+x; %from part a
uy = -((1-mu)*y/p1^3) - (mu*y/p2^3) + y;  %from part a

uxx = diff(ux,x)
uyy = diff(uy,y)
uxy = diff(ux,y)
