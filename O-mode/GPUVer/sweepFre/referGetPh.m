% Decompose the amplitude and phase from the sweep signal
% Local fft and smooth are employed.
f1=5;  % the initial freq
s=4;   % sweep rate
fr=50; % Resonant freq
af=[]; % amplitude
pf=[]; % phase
k1=0.02; % damping ratio
df=0.01; % freq interval
for fa=40:df:60
   t1=60/s*log2(fa/f1);
   t2=60/s*log2((fa+df)/f1);
   ta=t1:0.001:t2;
   N=length(ta);
   ft=f1*2.^(s/60*ta);
   A1=sin(2*pi*ft.*ta);
   lamb=ft/fr;
   B1=1./(1-lamb.^2+j*2*k1*lamb); % transfer function
   A2=abs(B1).*sin(2*pi*ft.*ta+angle(B1));
   ffreq=exp(-j*2*pi*(fa-400)*ta);  % freq shift for time domain
   spa=fft(ffreq.*A1);
   spb=fft(ffreq.*A2);
   spr=abs(spb./spa);
   spp=angle(spb./spa);
   k=ceil(N*0.001*400);
   af=[af,spr(k+1)];
   pf=[pf,spp(k+1)];
   figure(100)
   plot(ta,ffreq.*A1)
end
af=smooth(af,7);
pf=smooth(pf,7);  % Key trick
fa=40:df:60;
lamb=fa/fr;
bf=abs(1./(1-lamb.^2+j*2*k1*lamb));
subplot(2,1,1);
plot(fa,af,'r-',fa,bf,'b-.');
legend('Numeric Result','Theoretic Result');
title('Amplitude of Sweep Signal');
xlabel('f');
ylabel('A(f)');
subplot(2,1,2);
bpf=angle(1./(1-lamb.^2+j*2*k1*lamb));
plot(fa,180/pi*pf,'r-',fa,180/pi*bpf,'b-.');
legend('Numeric Result','Theoretic Result');
title('Phase of Sweep Signal');
xlabel('f');
ylabel('\Psi(f)');