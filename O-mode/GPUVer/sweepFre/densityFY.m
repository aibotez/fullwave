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

% [t,f1,f2] = GetFre(t,LaunchData,ReceiveData);
% [tlau,trec,Flau,Frec] = GetFre(t,LaunchData,ReceiveData);

% figure
% plot(t,LaunchData)
% y3 = bandpass(t,yn);
% [F,P]=Fft(t,y3);
% figure
% plot(F,P)
% figure
% plot(t,yn,t,y3)

c = 3*10^8;

% [f2,fb] = means(f2,fb,12);


% [F1,F2] = fdpek(t,y3,ReceiveData);
% F1=Fb;
f = 35*10^9;
lamda = c/f;
dt=lamda/40/c/1;
B = 0.3*f;
tm = 30000*dt;
k=2*B/tm;
dfdt = k;
% dfdt=dfdt;
% [tlau,trec,Flau,Frec] = GetFre(t,LaunchData,ReceiveData,f);
[tlau,trec,Flau,Frec] = GetFreSlidFFT(t,LaunchData,ReceiveData,f);
[f2,fb] = findfb(tlau,Flau,trec,Frec);
% [f2,fb] = means(f2,fb,20);
to0 = 9.441533121560042e-09;
to0 = 1.5e-9;
tao = (fb./dfdt)-to0;
figure
plot(f2,tao)
omega = 2*pi*f2;
to = tao;
L=0;
% yls=[];
yls(1)=0;
% xc = fidcut(omega(1));
for i =1:length(omega)
    i
    dx = intt(0,omega(i),to(i));
    ne1(i) = fidcutne(omega(i));
    yls(i) = dx+0.25;
end
figure
plot(yls,ne1)


function [F1,F2] = means(f1,fb,Lidx)
N = length(f1);
overlap = 0.5;
curidx=1;
% F=[];
F2=[];
F1=[];
while curidx+Lidx <= N
    idx1 = curidx;
    idx2 = idx1+Lidx;
    F1 = [F1 mean(f1(idx1:idx2))];
    F2 = [F2 mean(fb(idx1:idx2))];
    curidx = curidx+overlap*Lidx;
end

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

function [f2i,fb] = findfb(t1,f1,t2,f2)

% fb = [];
% for i = 1:length(t2)
% idx = find(t1>=t2(i));
% idx = idx(1);
% fb = [fb f1(idx)-f2(i)];
% end

fb = [];
f2i = [];
len=1;
for i = 1:len:length(t2)-len
idx = find(t1>=t2(i));
if length(idx)==0
    idx = length(t1);
end
idx = idx(1);
fb = [fb mean(f1(idx:idx+len))-mean(f2(i:i+len))];
f2i = [f2i mean(f2(i:i+len))];
end

end
function [tlau,trec,Flau,Frec] = GetFreSlidFFT(t,y1,y2,f)
dt = t(2)-t(1);
tin = t(1):dt/5:t(end);
y1=interp1(t,y1,tin,'spline');
y2=interp1(t,y2,tin,'spline');
t = tin;

[s,idx1] = findpeaks(y1);
[s,idx2] = findpeaks(y2);
y11 = y1(idx1);
y22 = y2(idx2);
t1 = t(idx1);
t2 = t(idx2);
F10 = [];
F20 = [];
t10=[];
t20=[];
% F1 = 1./diff(t1);
% F1=smooth(F1,200);
for i=2:1:length(t1)
    idx10 = find(t>=t1(i-1));
    idx1 = idx10(1);
    idx20 = find(t>=t1(i));
    idx2 = idx20(1);
    ti = t(idx1:idx2);
    yi = y1(idx1:idx2);
    [F,P]=Fft(ti,yi);
    idx = find(P==max(P));
    if length(idx)>1
        idx = idx(1);
    end
    F10 = [F10 F(idx)];
    t10 = [t10 mean(ti)];
end

for i=2:1:length(t2)
    idx10 = find(t>=t2(i-1));
    idx1 = idx10(1);
    idx20 = find(t>=t2(i));
    idx2 = idx20(1);
    ti = t(idx1:idx2);
    yi = y2(idx1:idx2);
    if max(abs(yi)) < 0.005
        yi = yi.*0;
    end
    [F,P]=Fft(ti,yi);
    idx = find(P==max(P));
    if length(idx)>1
        idx = idx(1);
    end
    F20 = [F20 F(idx)];
    t20 = [t20 mean(ti)];
end


idx11 = find(F10>=f);
idx22 = find(F20>=f);
t10 = t10(idx11);
t20 = t20(idx22);
F1 = F10(idx11);
F20 = F20(idx22);
idxmin = find(F20==min(F20));
t2 = t20(idxmin:end);
F2 = F20(idxmin:end);
F1=smooth(F1,200);
F2=smooth(F2,200);
figure
plot(t10,F1,t2,F2)

Frec = F2;
trec = t2;
Flau = F1;
tlau = t10;

m=0;
end
function [tlau,trec,Flau,Frec] = GetFre(t,y1,y2,f)

dt = t(2)-t(1);
tin = t(1):dt/5:t(end);
y1=interp1(t,y1,tin,'spline');
y2=interp1(t,y2,tin,'spline');
t = tin;

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


idx11 = find(F1>=f);
idx22 = find(F2>=f);
t1 = t1(idx11);
t2 = t2(idx22);
F1 = F1(idx11);
F2 = F2(idx22);
idxmin = find(F2==min(F2));
t2 = t2(idxmin:end);
F2 = F2(idxmin:end);


F1=smooth(F1,200);
F2=smooth(F2,200);
% figure
% plot(t1,F1,'r',t2,F2)
Frec = F2;
trec = t2;
Flau = F1;
tlau = t1;

m=0;
end



function [F11,F22] = fdpek(t,y1,y2)
idx = find(t>0.24e-9);
y1 = y1(idx);
y2 = y2(idx);
t = t(idx);
% y1 = abs(y1);
% y2 = abs(y2);
[s,idx1] = findpeaks(y1);
[s,idx2] = findpeaks(y2);
idx1i = find(y1(idx1)>0);
idx1 = idx1(idx1i);
idx2i = find(y2(idx2)>0);
idx2 = idx2(idx2i);
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
% F1=smooth(F1,1);

F11 = [];
F22 = [];
for i=3:length(t1)
    idx10 = find(t>=t1(i-2));
    idx1 = idx10(1);
    idx20 = find(t>=t1(i));
    idx2 = idx20(1);
    ti = t(idx1:idx2);
    yi = y1(idx1:idx2);
    [F,P]=Fft(ti,yi);
    idx = find(P==max(P));
    if length(idx)>1
        idx = idx(1);
    end
    F11 = [F11 F(idx)];
    
    yi = y2(idx1:idx2);
    [F,P]=Fft(ti,yi);
    idx = find(P==max(P));
    if length(idx)>1
        idx = idx(1);
    end
    F22 = [F22 F(idx)];
    
end
F22=smooth(F22,20);

% figure
% plot(1:length(F11),F22)
% 
% figure
% plot(t,y1,'b',t1,y11,'r.')
% figure
% plot(1:length(F1),F1)

% for i=2:length(t1)
%     idx10 = find(t>=t1(i-1));
%     idx1 = idx10(1);
%     idx20 = find(t>=t1(i));
%     idx2 = idx20(1);
%     ti = t(idx1:idx2);
%     yi = y2(idx1:idx2);
%     [F,P]=Fft(ti,yi);
%     idx = find(P==max(P));
%     if length(idx)>1
%         idx = idx(1);
%     end
%     F20 = [F20 F(idx)];
% end
% 
% % F20 = 1./diff(t2);
% % dfdt = 3.6000e+17;
% % tao = F1./dfdt;
% % figure
% % plot(t1(2:end),tao)
% % c = 3*10^8;
% % L = c*tao./2;
% % figure
% % plot(t1(2:end),L)
% 
% % for i =2:length(y22)
% %     idx1 = i-1;
% %     idx2 = i;
% %     dt = t2(idx2)-t2(idx1);
% %     F20 = [F20 1/dt];
% % end
% % F20=smooth(F20,10);
% ts1 = t1(2:end);
% idx = find(ts1>1*10^-9);
% ts = ts1(idx);
% F20 = F20(idx);
% F1 = F1(idx);
% 
% F2=smooth(F20,20);

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
