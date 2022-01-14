



clear;
dt=8.3333e-13;
y0 = 0.025;
c = 3*10^8;
t0 = 2.2e-9;
dfdt = 3.6000e+17;
file='result/FullData.mat';
[R,Z,ez,t,LaunchData,ReceiveData,phase,phaset,ReflI,ReflQ]=readdata(file);
fs = 1/(t(2)-t(1));
yn = LaunchData.*ReceiveData;

% figure
% plot(t,ReceiveData)
% idx = find(t>t0);
% y1 = LaunchData(idx-min(idx)+1);
% y2 = ReceiveData(idx);
% t1 = t(idx-min(idx)+1);
% yn1 = y1.*y2;



% f = [30:1:45].*10^9;
% omg = 2*pi.*f;
% for i = 1:length(f)
%     i
%     [ta(i) xc] = test3(omg(i));
%     dx(i) = intt(0,omg(i),ta(i))/2+0.25;
%     ne(i) = fidne(xc);
% end
% figure
% plot(f,ta)
% figure
% plot(dx,ne)

% f = 1*30*10^9;
% omg = 2*pi*f;
% ta = test3(omg)
% to=7.4538e-10;
% dx = intt(0,omg,to)
% pha = test0(omg);

% [F,P]=Fft(t,yn);
% figure
% plot(F,P)
y3 = bandpass(t,yn);


fftpoint=2048;
[f,P,time]=STFT(t,ReceiveData,fftpoint);
% figure
% imagesc(time,f,P')

c = 3*10^8;
% figure
% plot(t,y3)
[F,P]=Fft(t,y3);
% figure
% plot(F,P)

to0 = 9.441533121560042e-09;
to0 = 1.5e-9;
% to0 = 1.8e-08;
% [F1,F2] = fdpek(t,y3,ReceiveData);
[F1,F2] = fdpek(t,LaunchData,ReceiveData);
figure
plot(F2,F1)

Fb = F1;
fftpoint=1024;
% [F1,F2,Fb] = sft(t,LaunchData,ReceiveData,fftpoint);
% figure
% plot(F2,Fb)
% figure
% plot(F2,F1)
% F1 = smooth(F1,20);
F1=Fb;
f = 35*10^9;
lamda = c/f;
dt=lamda/40/c/1;
B = 0.3*f;
tm = 30000*dt;
k=2*B/tm;
dfdt = k;
% dfdt=dfdt;
tao = 1*(F1./dfdt)-to0;

% tao = smooth(tao,10);
% figure
% plot(F2,smooth(tao,10))
% tao = F1./dfdt;
L = c*tao./2;
% figure
% plot(F2,L,'.')
idx = find(F2>=2.1*10^10);


omega = 2*pi*F2(idx);
to = tao(idx);
% to = to-to(1);
figure
plot(omega./2./pi,to)
L=0;
% yls=[];
yls(1)=0;
xc = fidcut(omega(1));
for i =1:length(omega)
    i
    dx = intt(0,omega(i),to(i));
    ne1(i) = fidcutne(omega(i));
    yls(i) = dx+0.25;
%     L = L+dx; 
%     yls = [yls dx];
end
% figure 
% plot(1:length(yls),yls)
% figure
% plot(omega./2./pi,yls)
figure
plot(yls,ne1,'.')


fftpoint=2048;
[F0,F1]=STFT1(t,y3,fftpoint,ReceiveData);

figure
plot(F0,F1,'.')

[F,P]=Fft(t,y3);
figure
plot(F,P)
figure
plot(t,y3)
% figure
% plot(t,y3)
% [f,ph] = GetP(y3,fs);
% dph = diff(ph);
% figure
% plot(f(2:end),dph)


% fftpoint = 4096;
% [f,P,time]=STFT(t,yn,fftpoint);
% figure
% imagesc(time,f,P')
function x = fidcut(omeg)
x = 0.25:0.001:0.5;
x1=2.35 - x+0.25;
     k = -8/(2.35-1.9);
     b = -2.35*k;
     ne = (k.*x1+b).*10^19;
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
function pha = test0(omeg)
c = 3*10^8;
k = 8*10^19/(0.5-0.225);
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
    W2 = -c*e^2*k/epsilon0/me;
    pha = (4/(3*W2))*omeg^2/c-pi/2;
    pha = 2*(4/(3*W2))*omeg/c;
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
    syms x
    Kr = (omeg^2/c-(k*(2.35 - x)+b)*10^19*e^2/c/epsilon0/me)^0.5;
    pha = 2*int(Kr,ra,rc)-pi/2;
end
function [ta,rc] = test3(omeg)
c = 3*10^8;
%      x = 0.25:0.001:0.5;
%      x1=2.35 - x+0.25;
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
    syms x
    Kr = omeg/(omeg^2-(k*(2.35 - x+0.25)+b)*10^19*e^2/epsilon0/me)^0.5;
    ta = 2*int(Kr,ra,rc)/c;
  
end

function ne = fidne(xc)
x1=2.35 - xc+0.25;
k = -9/(2.35-1.9);
b = -2.35*k;
ne = (k.*x1+b).*10^19;
end
function ne1 = fidcutne(omeg)
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
    ne1 = omeg^2*epsilon0*me/e/e;
end
function dx = intt(omegao,omegap,to)
c = 3*10^8;
% fomg = @(omega) to./((omegap^2-omega.^2).^0.5);
% y = integral(fomg,omegao,omegap);
syms omega
fomg = to/((omegap^2-omega^2)^0.5);
y = int(fomg,omegao,omegap);
dx = y*c/(1*pi)/(2);
end
function [f,P,time]=STFT(t,y,fftpoint)
fs = 1/(t(2)-t(1));
N = length(y);
curidx=1;
time=[];
overlap=0.1;
% F=[];
P=[];
while curidx+fftpoint <= N
    idx1 = curidx;
    idx2 = idx1+fftpoint;
    ys = y(idx1:idx2);
    ts = t(idx1:idx2);
    time=[time mean(ts)];
    [f,p]=Fft(ts,ys);
    P = [P;p];
    curidx = curidx+overlap*fftpoint;
%     figure
%     plot(f,p)
    m=0;
end
end
function [F,P]=Fft(t,y)
dt = t(3)-t(2);
t1 = t.*1.e9;
idx = find(t1>=0);%2.73
FS=1/dt;
sample=y(idx);
t=t(idx);
N=length(sample);
Y = fft(sample);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
F = FS*((0):(N/2))/N;
P=P1;
end
function [R,Z,ez,t,LaunchData,ReceiveData,phase,phaset,ReflI,ReflQ]=readdata(file)
    FullDatas = load(file);
    ez = FullDatas.ez;
    LaunchData = FullDatas.LaunchData;
    ReceiveData = FullDatas.ReceiveData;
    R = FullDatas.R;
    Z = FullDatas.Z;
    t = FullDatas.t;
    phase = FullDatas.phase;
    phaset = FullDatas.tphase;
    ReflI = FullDatas.ReflI;
    ReflQ = FullDatas.ReflQ;
end

function y3 = bandpass(t,y)
%滤波
fp=5*10^9;%滤掉频率超过25Hz的信号
Fs=1/(t(2)-t(1));
fs = Fs/4;
wp=2*fp/Fs;ws=2*fs/Fs;rp=0.1;rs=60;   %DF指标（低通滤波器的通、阻带边界频）
[N,wp]=ellipord(wp,ws,rp,rs); %调用ellipord计算椭圆DF阶数N和通带截止频率wp
[B,A]=ellip(N,rp,rs,wp);      %调用ellip计算椭圆带通DF系统函数系数向量B和A
% y2=filter(B,A,y1); %滤波
y3=filtfilt(B,A,y);
end

function[F1,F2,Fb] = sft(t,y1,y2,fftpoint)
fftpoint=2^10;
fs = 1/(t(2)-t(1));
N = length(t);
curidx=1;
time=[];
overlap=0.9;
% F=[];
P=[];
j=1;
while curidx+fftpoint <= N
    idx1 = curidx;
    idx2 = idx1+fftpoint;
    ys1 = y1(idx1:idx2);
    ts = t(idx1:idx2);
    [f1,p1]=Fft(ts,ys1);
    idx = find(p1==max(p1));
    if length(idx)>1
       idx = idx(1); 
    end
    F1(j) = f1(idx);
    
    ys2 = y2(idx1:idx2);
    ts = t(idx1:idx2);
    [f2,p2]=Fft(ts,ys2);
    idx = find(p2==max(p2));
    if length(idx)>1
       idx = idx(1); 
    end
    F2(j) = f2(idx);

    curidx = curidx+overlap*fftpoint;
    j=j+1;
%     figure
%     plot(f,p)
    m=0;
end
F1s = smooth(F1,20);
F2s = smooth(F2,20);
Fb = F1s-F2s;
figure
plot(1:length(F1),Fb)
% figure
% plot(1:length(F1),Fb)
m=0;
end

function [F1,F2] = fdpek(t,y1,y2)
idx = find(t>2.4e-9);
y1 = y1(idx);
y2 = y2(idx);
t = t(idx);
% y1 = abs(y1);
% y2 = abs(y2);
[s,idx1] = findpeaks(y1);
[s,idx2] = findpeaks(y2);
y11 = y1(idx1);
y22 = y2(idx2);
t1 = t(idx1);
t2 = t(idx2);
% figure
% plot(t1,y11,'r.')
F1 = [];
F20 = [];
% F2 = [];
% for i =2:length(y11)
%     idx1 = i-1;
%     idx2 = i;
%     dt = t1(idx2)-t1(idx1);
%     F1 = [F1 1/dt];
% end
F1 = 1./diff(t1);
F1=smooth(F1,1);
% figure
% plot(t,y1,'b',t1,y11,'r.')

for i=2:length(t1)
    idx10 = find(t>=t1(i-1));
    idx1 = idx10(1);
    idx20 = find(t>=t1(i));
    idx2 = idx20(1);
    ti = t(idx1:idx2);
    yi = y2(idx1:idx2);
    [F,P]=Fft(ti,yi);
    idx = find(P==max(P));
    if length(idx)>1
        idx = idx(1);
    end
    F20 = [F20 F(idx)];
end

% F20 = 1./diff(t2);
% dfdt = 3.6000e+17;
% tao = F1./dfdt;
% figure
% plot(t1(2:end),tao)
% c = 3*10^8;
% L = c*tao./2;
% figure
% plot(t1(2:end),L)

% for i =2:length(y22)
%     idx1 = i-1;
%     idx2 = i;
%     dt = t2(idx2)-t2(idx1);
%     F20 = [F20 1/dt];
% end
% F20=smooth(F20,10);
ts1 = t1(2:end);
idx = find(ts1>1*10^-9);
ts = ts1(idx);
F20 = F20(idx);
F1 = F1(idx);

F2=smooth(F20,20);
% figure
% plot(ts,F20,'.')
% figure
% plot(t2(2:end),smooth(F20,30),'.')
% for i =1:length(F1)
%     ti = t1(i);
%     idx = find(t2>=ti);
%     idxx = idx(1);
%     F2(i) = F20(idxx);
% end

% figure
% plot(F2,F1,'.');
m=0;
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
function [F0,F1]=STFT1(t,y,fftpoint,y2)
fs = 1/(t(2)-t(1));
N = length(y);
curidx=1;
time=[];
overlap = 0.3;
% F=[];
P=[];
F0=[];
F1=[];
while curidx+fftpoint <= N
    idx1 = curidx;
    idx2 = idx1+fftpoint;
    ys = y(idx1:idx2);
    y2s = y2(idx1:idx2);
    ts = t(idx1:idx2);
    
%     xi = ts(1):1/(1*fs):ts(end)+1*(ts(end)-ts(1));
%     ys = interp1(ts,ys,xi, 'spline');
%     y2s = interp1(ts,y2s,xi, 'spline');
%     ts = xi;
    
    time=[time mean(ts)];
    
    
    [f1,p1]=Fft(ts,ys);
    [f2,p2]=Fft(ts,y2s);
    
%     figure(100)
%     subplot(2,1,1)
%     plot(ts,ys)
%     subplot(2,1,2)
%     plot(f1,p1)
    
    idxx = find(p1<0.0001);
    p1(idxx) = 0;
    idxx = find(p2<0.0001);
    p2(idxx) = 0;
    
    dinffidx = find(p1 == max(p1));
    if length(dinffidx)>1
       dinffidx=1; 
    end
    F1 = [F1 f1(dinffidx)];
    
    dinffidx = find(p2 == max(p2));
    if length(dinffidx)>1
       dinffidx=1; 
    end
    F0 = [F0 f2(dinffidx)];
    
    curidx = curidx+overlap*fftpoint;
%     figure
%     plot(f,p)
    m=0;
end
end
