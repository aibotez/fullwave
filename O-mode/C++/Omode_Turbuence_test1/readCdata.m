% [f2,P2]=Fft(time,ReflWave);
% figure
% plot(f2,P2)
clear;
c= 3*10^8;
v = c/20;

% clear
file='Omode_Turbuence_test1/';
[R,Z,ez1,t,ReceiveData1,ReceiveData2,phase1,phase2,phaset]=readdata(file);
% figure
% plot(t,ReceiveData1,'r',t,ReceiveData2,'b')
% figure
% plot(phaset,phase1,'r',phaset,phase2,'b')
fs = 1/(phaset(2)-phaset(1));
idx = find(phaset>=1e-8);
% y1 = sin(2*pi*20^6*phaset);
% y2 = sin(2*pi*20^6*(phaset+phaset(200)));
phase1 = phase1(idx)-mean(phase1(idx));
phase2 = phase2(idx)-mean(phase2(idx));
[to,po,delay] = TimeDelay(phase1,phase2,phaset(idx));
% muti = phase1(idx).*phase2(idx);
% sum(muti)
figure
plot(to,po)
figure
% plot(phaset,y1,'r',phaset,y2,'b')
plot(phaset(idx),phase1,'r',phaset(idx),phase2,'b')
[Coherence,Crossphase,f,Pxy,Pxx,Pyy]=coherence_me(phase1(idx),phase2(idx),fs);
% ezf = abs(ez-ez1);
% figure
% imagesc(R,Z,ezf')
    figure
    ax(1)=subplot(3,1,1);
         plot(f,log10(Pxx),'b')
         hold on
         plot(f,log10(Pyy),'r')
    ax(2)=subplot(3,1,2);
         plot(f,Coherence,'-')
    ax(3)=subplot(3,1,3);
         plot(f,Crossphase,'');
         linkaxes(ax,'x');



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


function [Rd,Zd,ezd,td,ReceiveData1d,ReceiveData2d,phase1d,phase2d,phasetd]=readdata(file)
    ez = [file 'Ez.mat'];
    ezd=load(ez);
    ezd = ezd.ez;
    ReceiveData1 = [file 'ReceiveData1.mat'];
    ReceiveData1d=load(ReceiveData1);
    ReceiveData1d = ReceiveData1d.ReceiveData1;
    ReceiveData2 = [file 'ReceiveData2.mat'];
    ReceiveData2d=load(ReceiveData2);
    ReceiveData2d = ReceiveData2d.ReceiveData2;
    R = [file 'R.mat'];
    Rd=load(R);
    Rd = Rd.R;
    Z = [file 'Z.mat'];
    Zd=load(Z);
    Zd = Zd.Z;
    t = [file 't.mat'];
    td=load(t);
    td = td.t;
    
    phase1 = [file 'phase1.mat'];
    phase1d=load(phase1);
    phase1d = phase1d.phase1;
    
    phase2 = [file 'phase2.mat'];
    phase2d=load(phase2);
    phase2d = phase2d.phase2;
    
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
function [Coherence,Crossphase,f,Pxy,Pxx,Pyy]=coherence_me(x1,x2,fs)
%========================
%calculate the power spectrum density of signal x
%fs is sampling frequency, fftpoint is fft point.
%if x is vector , the function returns the PSD of x
%========================
a=length(x1);
ffty1=fftshift(fft(x1));
ffty2=fftshift(fft(x2));
Pxy=ffty1.*conj(ffty2)/a;
Pxx=abs(ffty1).^2/a;
Pyy=abs(ffty2).^2/a;
Crossphase=angle(Pxy);
Coherence=abs(Pxy)./sqrt(abs(Pxx).*abs(Pyy));
f=(-a/2:(a/2-1))*(fs/a);
end

function [to,po,delay] = TimeDelay(y1,y2,t)
[po,to]=xcorr(y1,y2);
dt = t(2)-t(1);
to = to*dt;
[value,index]=max(po);
delay = to(index);
end