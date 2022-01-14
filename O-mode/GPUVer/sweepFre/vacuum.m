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

% figure
% plot(t,ReceiveData)
y3 = bandpass(t,ReceiveData);
[F,P]=Fft(t,y3);
figure
plot(F,P)
% figure
% plot(t,yn,t,y3)

c = 3*10^8;
% [F1,F2] = fdpek(t,y3,ReceiveData);
% F1=Fb;
f = 30*10^9;
lamda = c/f/1;
dt=lamda/40/c/1;
B = 0.3*f;
tm = 30000*dt;
k=2*B/tm;
dfdt = k;
[tlau,trec,Flau,Frec] = GetFre(t,LaunchData,ReceiveData,f);
[f2,fb] = findfb(tlau,Flau,trec,Frec);
% dfdt=dfdt;
tao = (fb./dfdt);
tao1 = (fb./dfdt);
figure
plot(f2,tao)
% tao = (F1./dfdt);
L = tao*c/2 +y0;
% idx = find(f2>=f);
% meanL = mean(L(idx));
% dL = meanL-0.2
figure
% plot(f2,tao1);
plot(f2,L);


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

% figure
% plot(t,y2,t2,y2(idx2),'r.')
% figure
% plot(t1(2:end),F1,'r',t2(2:end),F2)

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
% figure
% plot(t2(2:end),diff(F2),'.')
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
