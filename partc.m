%304 Project 1 - Groupmembers: Thomas Waltz, Ankit Gupta, and Kameron
%Metcalf

clc, clear

y1=sqrt(3)/2; %from part b


%sun-earth

muse = 3.0039*10^-7;
lsex = [0.995363, 1.004637, -1.00001, 0.5-muse, 0.5-muse];
lsey = [0,0,0,y1,-y1];
for i= 1:5
 ees(:,i) = ev(muse, lsex(1,i) ,lsey(1,i));
end

 %Earth-Moon
muem = 1.2151*10^-2;
lemx = [0.836915, 1.15568, -1.00506, 0.5-muem, 0.5-muem];
lemy = [0,0,0,y1,-y1];
for i= 1:5
 eem(:,i) = ev(muem, lemx(1,i) ,lemy(1,i));
end

%Saturn-Titan
must = 2.366*10^-4;
lstx = [0.9575, 1.0425, -1.0001, 0.5-must, 0.5-must];
lsty = [0,0,0,y1,-y1];
for i= 1:5
 est(:,i) = ev(must, lstx(1,i) ,lsty(1,i));
end
ees
eem
est


%function ev gets us U's equation required for calculation
function ev = ev(mu,x,y) 

Uxx=(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)) + 1; %from derivates.m
Uyy=(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1; %from derivates.m
Uxy=(3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2)^(5/2); %from derivates.m


A = [0, 0, 1, 0; 0, 0, 0, 1;Uxx, Uxy, 0, 2;Uxy, Uyy, -2, 0];
ev = eig(A);
end
