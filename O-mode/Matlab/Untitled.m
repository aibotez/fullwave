

t0 = (0:nsteps-1).*dt;
fs = 1/(t0(2)-t0(1));
sigr1 = Erec1;
sign = Erec1+Erec2;
[tn,phase] = cauphase(t0,Erec1,Erec2,sign,f,fs);
figure
plot(tn,phase)


function [tn,ph] = cauphase(t,y1,y2,yn,f,fs)
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