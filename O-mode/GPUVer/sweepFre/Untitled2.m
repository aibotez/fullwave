% [f2,P2]=Fft(time,ReflWave);
% figure
% plot(f2,P2)

% clc,close all,clear all;
% f0      = 24e9; %% 中心频率24GHZ
% fm      = 100;      %% 三角调制波频率100Hz
% tm      = 1/fm;     %% 三角调制波周期，单位s
% N       = 1024;     %% FFT变换点数
% B       = 250e6;    %% 调制带宽250MHZ
% c       = 3e8;      %% 光速, unit: m/s
% k=2*B/tm;            %% 三角调制波斜率  
% Fs      = 2000e3;    %% AD 采样速率, unit: Hz
% Ts=1/Fs; %%采样周期
% SNR=30;                                       %信噪比
% n=0:N-1;
% t=n/Fs;  
% %t=n/N;
% r=120;                  %实际距离,可自己修改
% t0=2*r/c;
% f=n*Fs/N;    %频率序列
% S_t=cos(2*pi*f0*t+pi*k*t.^2);            %发射信号
% S_r=cos(2*pi*f0*(t-t0)+pi*k*(t-t0).^2);  %接收信号
% S_r=awgn(S_r,SNR);
% Hu=S_t.*S_r+randn(size(t));    
% figure
% plot(t,S_t)
% clear
clear;
% f = 2;
% t = 0:0.01:10;
% B = 1*f;
% tm = 10;
% k=2*B/tm;
% % k=0;
% y = sin(2*pi*f*t+pi*k*t.^2);
% [F,P]=Fft(t,y);
% figure
% plot(t,y)
% figure
% plot(F,P)




clear
nstep=30000;
dt = 8.3333e-13;
dfdt = 3.6000e+17;
t0 = 2.23e-9;%1.9e-9
idxoffet = fix(t0/dt);
file='result/FullData.mat';
[R,Z,ez2,t,LaunchData,ReceiveData,phase,phaset,ReflI,ReflQ]=readdata(file);
[f2,P2]=Fft(t,ReceiveData);
% figure
% plot(f2,P2)


f = 30*10^9;
[ft,Launhsignal] = fres(t0,dt,nstep);
% [F,P]=Fft(t,ReceiveData);
% figure
% plot(F,P)

fs = 1/(t(2)-t(1));
[f,ph1] = GetP(Launhsignal,fs);
[f,ph2] = GetP(ReceiveData,fs);
dph = abs(diff(ph2)./diff(f));

fidx = find(f>=3.16*10^10 & f<=3.6*10^10);
fi = f(fidx);
dph1 = dph(fidx);
dpidx = find(dph1<1.5*10^-7);
dph2 = dph1(dpidx);
fsi = [];
for i =1:length(dph2)
    idx = find(dph2==dph2(i));
    fsi = [fsi fi(idx)];
end
omega = 2*pi*fsi;
dph2 = dph2./(2*pi);
figure
plot(fsi,dph2)
yls=[];
L=0;
to1 = dph2;
to0 = 9.441533121560042e-09;
to = to1-to0;
figure
plot(omega,to)
for i =2:length(omega)
    i
    dx = intt(omega(1),omega(i),to(i));
    L = L+dx;
    yls = [yls dx];
end
figure
plot(1:length(yls),yls)


figure
plot(f(2:end),dph,'.')
figure
plot(f,ph2)

% figure
% plot(t,Launhsignal)
fftpoint=1024;
[f1,P,time]=STFT(t,ReceiveData,fftpoint);
figure
imagesc(time,f1,P')
df = 9.6*10^5;

s1 = ReceiveData(idxoffet:end);
s2 = Launhsignal(1:end-idxoffet+1);
ft1 = ft(1:end-idxoffet+1);
sn = s1+s2;
ti = t(1:end-idxoffet);
[tn,ph,fn] = cauphase1(ti,s1,s2,sn,ft1);
% [F,P]=Fft(t,ReceiveData);
figure
% plot(t,Launhsignal,'b',t,ReceiveData,'r')
plot(tn,ph)
% [tn,ph] = cau(t,y1,y2,yn,f);
% df = abs(ReceiveData0-ReceiveData0);
% figure
% plot(t,df)
phs = smooth(ph,20);
phds = diff(phs);
phd = diff(ph);
fnd = diff(fn);
phds1 = smooth(phd,20);
figure
plot(fn(2:end),abs(phd)./df,'.')
[tim,ys]=meanph(tn(2:end),abs(phd)./df);
figure
plot(tim,ys)


[tim,to]=meanph(tn(2:end),abs(phd)./(2*pi*df));
omega = 2*pi*fn;
% to = abs(phd)./(2*pi*df);
yls=[];
L=0;
for i =6:length(tim)
    i
    dx = intt(omega(i-1),omega(i),to(i-1));
    L = L+dx;
    yls = [yls L];
end

% [f1,P1]=Fft(t0,Erec);
tn = phaset;
% figure
% plot(tn,phase)
idx = find(isnan(phase));
% for i=idx
%     phase(i) = (phase(i+1)+phase(i-1))/2;
% end
y2 = cos(phase)+(-1).^0.5.*sin(phase);
% [f2,P2]=Fft(t,ReceiveData);
% figure
% plot(f2,P2)

figure
subplot(3,1,1)
plot(t,ReceiveData,'r',t,LaunchData,'b')
subplot(3,1,2)
plot(tn,phase)
subplot(3,1,3)
plot(f2,P2)

% figure
% plot(t,LaunchData,'b',t,ReceiveData,'r')
% figure
% plot(f2,P2)
function dx = intt(omegao,omegap,to)
c = 3*10^8;
% fomg = @(omega) to./((omegap^2-omega.^2).^0.5);
% y = integral(fomg,omegao,omegap);
syms omega
fomg = to/((omegap^2-omega^2)^0.5);
y = int(fomg,omegao,omegap);
dx = y*c/(1*pi);
end
function [time,ys]=meanph(tn,phase)
N = length(tn);
width = 8;
overlap=0.5;
curidx=1;
ys=[];
time=[];
while curidx+width <= N
    idx1 = curidx;
    idx2 = idx1+width;
    ys = [ys mean(phase(idx1:idx2))];
    ts = tn(idx1:idx2);
    time = [time mean(ts)];
    curidx = curidx+fix(overlap*width);
end

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
function [f,P,time]=STFT(t,y,fftpoint)
fs = 1/(t(2)-t(1));
N = length(y);
overlap = 0.5;
curidx=1;
time=[];
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
    curidx = curidx+fix(overlap*fftpoint);
%     figure
%     plot(f,P)
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



% n=0:N-1;
% y=fft(sample,N);    %对信号进行快速Fourier变换
% mag1=abs(y)./N;     %求得Fourier变换后的振幅
% mag=smooth(mag1,50);
% f=n*FS/(N-1);    %频率序列
% F = FS*(0:(N/2))/N;
% 
% % F = f(1:N/2);
% P = mag1(1:N/2);
m=0;
end
function [wave]=GetWaveOr(phaset,t)
f = 30*10^9;
dt =8.3333e-13;
df = 0.1*10^9;
n = [1:1:30000];
ft = df.*fix(n*dt/dti)+f;
tft = n.*dt;



dt1 = t(2)-t(1);
dt2 = phaset(2)-phaset(1);
idx = fix(dt2/dt1);
% wave = 
end

function [ft,signal] = fres(t0,dt,nstep)
nt0 = fix(t0/dt);
f = 30*10^9;
% dti = 300/f;
% n = [nt0:1:nstep];
% df = 0.1*10^9;
% ft = ones(1,nt0-1).*f;
% fti = df.*fix(n*dt/dti)+f;
% ft = [ft,fti];
% 
% ns = [1:1:nstep];
% td=150;
signal = [];
% 
% B = 1*f;
% tm = 30000*dt;
% k=2*B/tm;
% k=0;
% for i=1:length(ft)
%     T = i-1;
%     signal = [signal cos(2 * pi * f * dt * T+pi*k*(dt*T).^2)];
% end

B = 0.3*f;
tm = 30000*dt;
k=2*B/tm;
t = [0:30000-1].*dt;
signal = cos(2 * pi * f .* t+pi*k.*(t).^2);
ft = f+k.*t./2;
m=0;
% figure
% plot((1:nstep).*dt,ft./10^10)

end
function [tn,ph,fn] = cauphase1(t,y1,y2,yn,f)
    I1=[];
    I2=[];
    In=[];
    tn=[];
    fn = [];
    N = length(t);
    fs = 1/(t(2)-t(1));
    leni=1;
    Tidx = 1*fix(fs/f(1));
    while leni+Tidx <=N
        idx1 = leni;
        idx2 = idx1+Tidx;
        I1 = [I1 max(abs(y1(idx1:idx2)))^2];
        I2 = [I2 max(abs(y2(idx1:idx2)))^2];
        In = [In max(abs(yn(idx1:idx2)))^2];
        tn=[tn t(idx1)];
        fn = [fn f(idx1)];
        leni = idx2;
        Tidx = 1*fix(fs/f(idx1));
    end
   cosph=(In-I1-I2)./(2.*I1.^0.5.*I2.^0.5);
   ph = acos(cosph);

end
function [tn,ph] = cauphase(t,y1,y2,yn,f)
    fs = 1/(t(2)-t(1));
    T = 1/f;
    widx = 1*fix(T*fs);
    I1=[];
    I2=[];
    In=[];
    tn=[];
    for i =1:widx:length(t)-widx
        I1 = [I1 max(abs(y1(i:i+widx)))^2];
        I2 = [I2 max(abs(y2(i:i+widx)))^2];
        In = [In max(abs(yn(i:i+widx)))^2];
        tn=[tn t(i)];
    end
   cosph=(In-I1-I2)./(2.*I1.^0.5.*I2.^0.5);
   ph = acos(cosph);
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