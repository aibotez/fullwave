% [f2,P2]=Fft(time,ReflWave);
% figure
% plot(f2,P2)


% clear
file='result/FullData.mat';
[R,Z,ez,t,LaunchData,ReceiveData1,ReceiveData2,phase1,phase2,phaset,ReflI1,ReflQ1,ReflI2,ReflQ2]=readdata(file);
% [R,Z,ez2,t,LaunchData,ReceiveData,phase,phaset,ReflI,ReflQ]=readdata(file);
dt = t(2)-t(1);
[a,b] = findpeaks(ReceiveData1);
figure
plot(b.*dt,a)
% lev=[0.05,0.1,0.2,0.4,0.8,1.6,3.2];
% figure
% contour(R,Z,ez',lev);
% figure
% imagesc(R,Z,ez')
% set(gca,'YDir','normal')
% df = abs(ReceiveData0-ReceiveData0);
% figure
% plot(t,df)
idx = find(phaset>=1e-8);
[to,po,delay] = TimeDelay(phase1(idx),phase2(idx),phaset);
delay
figure
subplot(3,1,1)
plot(t,ReceiveData1,'r',t,ReceiveData2,'b')
subplot(3,1,2)
plot(phaset,phase1,'r',phaset,phase2,'b')
subplot(3,1,3)
plot(to,po)



% [f1,P1]=Fft(t0,Erec);
tn = phaset;
% figure
% plot(tn,phase)
idx = find(isnan(phase));
% for i=idx
%     phase(i) = (phase(i+1)+phase(i-1))/2;
% end
y2 = cos(phase)+(-1).^0.5.*sin(phase);
[f2,P2]=Fft(tn,y2);
% figure
% plot(f1,P1)

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


function [R,Z,ez,t,LaunchData,ReceiveData1,ReceiveData2,phase1,phase2,phaset,ReflI1,ReflQ1,ReflI2,ReflQ2]=readdata(file)
    FullDatas = load(file);
    ez = FullDatas.ez;
    LaunchData = FullDatas.LaunchData;
    ReceiveData1 = FullDatas.ReceiveData1;
    ReceiveData2 = FullDatas.ReceiveData2;
    R = FullDatas.R;
    Z = FullDatas.Z;
    t = FullDatas.t;
    phase1 = FullDatas.phase1;
    phase2 = FullDatas.phase2;
    phaset = FullDatas.tphase;
    ReflI1 = FullDatas.ReflI1;
    ReflI2 = FullDatas.ReflI2;
    ReflQ1 = FullDatas.ReflQ1;
    ReflQ2 = FullDatas.ReflQ2;
end

function [F,P]=Fft(t,y)
dt = t(3)-t(2);
t1 = t.*1.e9;
idx = find(t1>=2.73);
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

function [to,po,delay] = TimeDelay(y1,y2,t)
y1 = y1-mean(y1);
y2 = y2-mean(y2);
[po,to]=xcorr(y1,y2,'coeff');
dt = t(2)-t(1);
to = to*dt;
[value,index]=max(po);
delay = to(index);
end