

clear;
f = 20*10^9;
lamda = 1/f/2;
dt=lamda/40;
[signal0,signal1,t,dfdt] = fres(f,dt);

[tlau,trec,Flau,Frec] = GetFre(t,signal0,signal1,f)

yn = signal0.*signal1;
y3 = bandpass(t,yn);
[F,P]=Fft(t,y3);
c = 3*10^8;
idx = find(P == max(P));
fof = F(idx);
tao = fof/dfdt
L = tao*c/2
% figure
% plot(F,P)
% dfdt




function [signal0,signal1,t,dfdt] = fres(f,dt)
B = 0.3*f;
tm = 30000*dt;
k=2*B/tm;
dfdt = k;
t = [0:30000-1].*dt;
signal0 = cos(2 * pi * f .* t+pi*k.*(t).^2);
c = 3*10^8;
ddt = fix(0.35/c/dt)*2;
tde = t-ddt*dt;
signal1 = cos(2 * pi * f .* tde+pi*k.*(tde).^2);
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
function [tlau,trec,Flau,Frec] = GetFre(t,y1,y2,f)


[s,idx1] = findpeaks(y1);

idx1 = idx1(1:end);
y11 = y1(idx1);
t1 = t(idx1);


for i=2:length(t1)
     idx10 = find(t>=t1(i-1));
    idx1 = idx10(1);
    idx20 = find(t>=t1(i));
    idx2 = idx20(1);
    ti = t(idx1:idx2);
    yi = y2(idx1:idx2);
    if max(abs(yi)) < 0.005
        y2(idx1:idx2) = 0;
    end
    
end


[s,idx2] = findpeaks(y2);
idx2 = idx2(1:end);
t2 = t(idx2);
F1 = 1./diff(t1);
F2 = 1./diff(t2);

figure
plot(t,y2,t2,y2(idx2),'r.')

idx11 = find(F1>=f);
idx22 = find(F2>=f);
t1 = t1(idx11);
t2 = t2(idx22);


F1 = F1(idx11);
F2 = F2(idx22);
% F1=smooth(F1,200);
% F2=smooth(F2,200);
figure
plot(t1,F1,'r',t2,F2)
Frec = F2;
trec = t2;
Flau = F1;
tlau = t1;


% [f2,fb] = findfb(tlau,F1,trec,F2);
% figure
% plot(t1(2:end),F1./10^9,t2(2:end),F2./10^9)
% figure
% plot(f2,fb)

m=0;
end