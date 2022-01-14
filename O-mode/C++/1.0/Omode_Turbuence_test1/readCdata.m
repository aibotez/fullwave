% [f2,P2]=Fft(time,ReflWave);
% figure
% plot(f2,P2)


% clear
file='Omode_Turbuence_test1/';
[R,Z,ez1,t,LaunchData,ReceiveData,phase,phaset]=readdata(file);


% ezf = abs(ez-ez1);
% figure
% imagesc(R,Z,ezf')
figure
plot(phaset,phase)



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
plot(t,ReceiveData)
subplot(3,1,2)
plot(tn,phase)
subplot(3,1,3)
plot(f2,P2)

% figure
% plot(t,LaunchData,'b',t,ReceiveData,'r')
% figure
% plot(f2,P2)


function [Rd,Zd,ezd,td,LaunchDatad,ReceiveDatad,phased,phasetd]=readdata(file)
    ez = [file 'Ez.mat'];
    ezd=load(ez);
    ezd = ezd.ez;
    LaunchData = [file 'LaunchData.mat'];
    LaunchDatad=load(LaunchData);
    LaunchDatad = LaunchDatad.LaunchData;
    ReceiveData = [file 'ReceiveData.mat'];
    ReceiveDatad=load(ReceiveData);
    ReceiveDatad = ReceiveDatad.ReceiveData;
    R = [file 'R.mat'];
    Rd=load(R);
    Rd = Rd.R;
    Z = [file 'Z.mat'];
    Zd=load(Z);
    Zd = Zd.Z;
    t = [file 't.mat'];
    td=load(t);
    td = td.t;
    
    phase = [file 'phase.mat'];
    phased=load(phase);
    phased = phased.phase;
    phaset = [file 'tphase.mat'];
    phasetd=load(phaset);
    phasetd = phasetd.tphase;
end

function [F,P]=Fft(t,y)
dt = t(3)-t(2);
t1 = t.*1.e9;
idx = find(t1<300 & t1>=2.73);
FS=1/dt;
sample=y(idx);
t=t(idx);
N=length(sample);
n=0:N-1;
y=fft(sample,N);    %对信号进行快速Fourier变换
mag1=abs(y)./1000;     %求得Fourier变换后的振幅
mag=smooth(mag1,50);
f=n*FS/(N-1);    %频率序列
F = f(1:N/2);
P = mag1(1:N/2);
m=0;
end