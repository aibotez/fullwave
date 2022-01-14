
clear
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
c = 3*10^8;
x = [0.25:0.01:0.5];

% ne=Caune(x);
% figure
% plot(x,ne)

x1=2.35 - x+0.25;
k = -0.6*8/(2.35-1.9);
b = -2.35*k;
ne0 = (k.*x1+b).^1;
ne0 = ne0.*10^19;
% figure
% plot(x,ne0)

f = [30:1:40].*10^9;
omg = 2*pi.*f;

B = (-k*e^2)/(epsilon0*me);
for i = 1:length(f)
   [ta(i) xc(i)] = test3(omg(i),0.25); 
end
x=2.35 - xc+0.25;
tac = (4*omg./c/B).*(omg.^2-e.^2.*(k.*x+b)./epsilon0./me).^0.5;
x1=2.35 - 0.25+0.25;
tac0 = (4*omg./c/B).*(omg.^2-e.^2.*(k.*x1+b)./epsilon0./me).^0.5;
figure
plot(f,tac-tac0)


% [ta0 xc0] = test3(omg(1),0.25);

for i = 1:length(f)
    i
    tac(i) = tatest(omg(i),0.25);
%     [ta(i) xc] = test3(omg(i),xc0);
%     dx(i) = intt(omg(1),omg(i),ta(i))/1.38+xc0;
%     ne(i) = fidcutne(omg(i));
end
figure
plot(f,tac)
% figure
% plot(f,dx)
% figure
% plot(dx,ne,'r',x,ne0,'b.')




function neg=Caune(y)
stpot = 0.25;
idx = find(y>=stpot);
ped_pos=2.25;  %pedestal position    
h=2.2; %height
w=0.06; %width
slope1=0.2; %slope of the core part  
slope2=-0.5;% slope of edge part
a0=[h/2 h/2 slope1 ped_pos w slope2];
neg = [1:length(y)].*0;
for j = idx
    dist = 2.35-y(j)+stpot;
    ne0 = MTANH(a0,dist);
    neg(j)=ne0;
end

end

function y=MTANH(a,x)  % MTANH function to generate the density profile for H-mode plasma
A=a(1);
B=a(2);
alpha=a(3);
x0=a(4);
w=a(5);
z=(x0-x)./w;
mtanh=((1+alpha*z).*exp(z)-(1-(0.04)*z).*exp(-z))./(exp(z)+exp(-z));
y=(A*mtanh+B)*1.e19;
if y<0
   y=0 
end
end
function ne1 = fidcutne(omeg)
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
    ne1 = omeg^2*epsilon0*me/e/e;
end
function pha = test2(omeg)
c = 3*10^8;
%      x = 0.25:0.001:0.5;
%      x1=2.35 - x;
     k = -9/(2.35-1.9);
     b = -2.35*k;
%      ne = (k.*x1+b).*10^19;
     ra = 0.25;
     rc = fidcut(omeg);
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
%     omegp = (ne.*e^2./epsilon0./me).^0.5;
%     Kr = ((omeg^2-omegp.^2)./c).^0.5;
% Kr = (omeg^2/c-(k*(2.35 - 0.25+0.25)+b)*10^19*e^2/c/epsilon0/me)^0.5
    syms x
    Kr = (omeg^2/c-(k*(2.35 - x+0.25)+b)*10^19*e^2/c/epsilon0/me)^0.5;
    pha = 2*int(Kr,ra,rc)-pi/2;
end
function [ta] = tatest(omeg,xc0)
     ra = xc0;
     rc = fidcut(omeg);
     k = -0.6*8/(2.35-1.9);
     b = -2.35*k;
c = 3*10^8;
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
    A = omeg^2-(b*e^2)/(epsilon0*me);
    B = (-k*e^2)/(epsilon0*me);
    x=2.35 - rc+0.25;
    tac = (4*omeg/c/B)*(omeg^2-e^2*(k*x+b)/epsilon0/me)^0.5;
    x=2.35 - ra+0.25;
    taa = (4*omeg/c/B)*(omeg^2-e^2*(k*x+b)/epsilon0/me)^0.5;
    ta = tac-taa;
end
function [ta,rc] = test3(omeg,xc0)
c = 3*10^8;
%      x = 0.25:0.001:0.5;
%      x1=2.35 - x+0.25;
     k = -0.6*8/(2.35-1.9);
     b = -2.35*k;
%      ne = (k.*x1+b).*10^19;
     ra = xc0;
     rc = fidcut(omeg);
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
%     omegp = (ne.*e^2./epsilon0./me).^0.5;
%     Kr = ((omeg^2-omegp.^2)./c).^0.5;
    syms x
    Kr = omeg/(omeg^2-(k*(2.35 - x+0.25)+b)^1*10^19*e^2/epsilon0/me)^0.5;
%     ta = 2*4*pi*pi*int(Kr,ra,rc)/c;
    ta = 2*1*int(Kr,ra,rc)/c;
  
end

function ne = fidne(xc)
x1=2.35 - xc+0.25;
k = -9/(2.35-1.9);
b = -2.35*k;
ne = (k.*x1+b).*10^19;
end
function x = fidcut(omeg)
x = 0.25:0.001:0.5;
x1=2.35 - x+0.25;
     k = -0.6*8/(2.35-1.9);
     b = -2.35*k;
     ne = (k.*x1+b).^1;
     ne = ne.*10^19;
         e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
omegp = (ne.*e^2./epsilon0./me).^0.5;
idx1 = find(omegp>=omeg);
idx2 = find(omegp<omeg);
if idx1(1)>=idx2(end)
    idx = idx1(1);
else
    idx = idx2(1);
end
x = x(idx);
end
function dx = intt(omegao,omegap,to)
c = 3*10^8;
% fomg = @(omega) to./((omegap^2-omega.^2).^0.5);
% y = integral(fomg,omegao,omegap);
syms omega
fomg = to/((omegap^2-omega^2)^0.5);
y = int(fomg,omegao,omegap);
dx = y*c/(1*pi)/(1);
end