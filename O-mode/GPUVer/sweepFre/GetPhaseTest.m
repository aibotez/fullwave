clear;
% t=1:1:2000;
% N=2000;
% fs=1000;
% for k=1:1:2000
%     x(k)=cos(2*pi*10*k/1000+pi*20*(k/1000)^2-pi/6);
% end
% % figure
% % plot(t,x);
% y=fft(x,2000);
% mag=abs(y);
% n=0:N-1;
% f=n*fs/N;
% 
% ph=2*angle(y);
% figure
% plot(f(1:N/2),ph(1:N/2));
% 
% ph=2*angle(y(1:N/2));
% 
% 
% f=n*fs/N;
% % plot(f(1:N/2),mag(1:N/2)*2/N);
% figure
% plot(f(1:N/2),ph(1:N/2));


% t=linspace(0,5,160);
% x=cos(2*pi*t-pi);
% figure
% plot(t,x)
% phase = GetPhase(x)
% X=fft(x);
% figure
% stem(abs(X));
% figure;
% stem(angle(X)/pi*180);



fs = 1000;
t = [0:1/fs:10];
t0 = 1;
ns = length(t);
N=ns;
dt = t(2)-t(1);
[ft,ylau] = fres(dt,ns);
yrefl = GeneReflS(ylau,t0,fs,dt);

% 
% [tn,phL] = GetPhases(t,ylau,ft,dt);
% [tn,phR] = GetPhases(t,yrefl,ft,dt);

[f,ph1] = GetP(ylau,fs);
[f,ph2] = GetP(yrefl,fs);

% y1=fft(ylau);
% mag=abs(y1);
% n=0:N-1;
% f=n*fs/N;
% ph1=2*angle(y1);
% 
% y2=fft(yrefl);
% mag=abs(y2);
% n=0:N-1;
% f=n*fs/N;
% ph2=2*angle(y2);

% 
% figure
% plot(f(1:N/2),ph2(1:N/2)-ph1(1:N/2));
% figure
% plot(f,mag)

dph = diff(ph2-ph1)./diff(f);
figure
plot(f(2:end),dph,'.')


figure
plot(tn,phR-phL)

figure
plot(t,yrefl,'b',t,ylau,'r')

function [ft,signal] = fres(dt,ns)
f = 5;
signal = [];
B = 1*f;
tm = ns*dt;
k=2*B/tm;
t = [0:ns-1].*dt;
signal = cos(2 * pi * f .* t+pi*k.*(t).^2);
ft = f+k.*t./2;
m=0;
end
function sig = GeneReflS(y,t0,fs,dt)
widx = fix(t0*fs);
N = length(y);
sig = [1:N].*0;
ns = N-widx+1;
[ft,signal] = fres(dt,ns);
sig(widx:end) = signal;
m=0;
end
function phase = GetPhase(s)
P = fft(s);
Pmax = max(P);
idx = find(P==Pmax);
if length(idx)>1
    idx = idx(1);
end
ph = angle(P);
phase = ph(idx);
end
function [tn,ph] = GetPhases(t,y,f,dt)
t0 = 1/(f(1));
fftpoint = fix(t0/dt);
N = length(y);
overlap = 0.5;
curidx=1;
tn=[];
% F=[];
ph=[];
while curidx+fftpoint <= N
    idx1 = curidx;
    idx2 = idx1+fftpoint;
    ys = y(idx1:idx2);
    ts = t(idx1:idx2);
    tn=[tn mean(ts)];
    phase = GetPhase(ys);
    ph = [ph phase];
    curidx = curidx+fix(overlap*fftpoint);
%     figure
%     plot(f,P)
    m=0;
end
end
function [f,ph] = GetP(y,fs)
N = length(y);
P=fft(y);
mag=abs(P);
n=0:N-1;
f1=n*fs/N;
ph1=2*angle(P);
f = f1(1:N/2);
ph = ph1(1:N/2);
end