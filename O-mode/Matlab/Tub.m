
clear;

c = 3.e8;   % the speed of light in vacuum
f = 30.e9;  % the frequency of the probing wave in vacuum
theta = 1*pi/180;%天线入射角
Dth = 0.3; %天线方向角
alpha = pi/2-theta; %天线喇叭口平面角度
AtR0=2.5; %天线顶点R坐标
AtZ0=-0.2; %天线顶点Z坐标
AtL = 0.05; %距离天线顶点的距离
rw = AtL;
w=2*AtL*tan(Dth);
w=0.01
% RM = AtR0-AtL*cos(theta);
% ZM = AtZ0+AtL*sin(theta);
lambda = c/f;
k0 = 2*pi/lambda;
ddx = lambda/20;
ddy=1*ddx;
dt=ddx/c/2;
x = 0:ddx:0.4;
y = 0:ddx:0.5;


% r0 = 1.9:0.001:2.6;
% z0 = 0:-0.001:1.9-2.6;
% ne=Caune1(x,y,ddx,ddy,theta);
ne=Caune(x,y);
ne0=ne;
[nx,ny] = size(ne);
w=w/2;
xatp = 0.2;
Ax1=0.245;
Ax2 =0.155;
[Anxl,Anyl]=setLaunchAntenna(xatp,1.5*AtL,w,1.3*AtL,ddx,ddy);
[Anxr1,Anyr1]=setReceiverAntenna(Ax1,1.5*AtL,w,1.3*AtL,ddx,ddy);
[Anxr2,Anyr2]=setReceiverAntenna(Ax2,1.5*AtL,w,1.3*AtL,ddx,ddy);
Anx=[Anxl Anxr1 Anxr2];
Any=[Anyl Anyr1 Anyr2];
% figure
% imagesc(x,y,ne');
% set(gca,'YDir','normal')
% xlabel('R (m)')
% ylabel('Z (m)')  
kstart=1;
kend=nx-100;
omega_pe=omega_pecal(ne);
t0=150.0;
nsteps=6000;
    ez=zeros(nx,ny);
    dz=zeros(nx,ny);
    hx=zeros(nx,ny);
    hy=zeros(nx,ny);
%     ix=zeros(1,nx);
    sx=zeros(nx,ny);
    sx1=zeros(nx,ny);
    sx2=zeros(nx,ny);
    
    ihx = zeros(nx,ny);
    ihy = zeros(nx,ny);
    
    %Calculat the PML parameters
    gi2 = ones(1,nx);
    gi3 = ones(1,nx);
    fi1 = zeros(1,nx);
    fi2 = ones(1,nx);
    fi3 = ones(1,nx);
    
    gj2 = ones(1,ny);
    gj3 = ones(1,ny);
    fj1 = zeros(1,ny);
    fj2 = ones(1,ny);
    fj3 = ones(1,ny);
    coeff_pml = 0.33;
    %Creat the PML 
    npml =8;
    for i = 1:npml
        xnum = npml - i;
        xd = npml;
        xxn = xnum / xd;
        xn = coeff_pml * xxn^3;
        
        gi2(i) = 1/(1+xn);
        gi2(nx - 1 -i) = 1/(1+xn);
        gi3(i) = (1-xn)/(1+xn);
        gi3(nx - 1 -i) = (1-xn)/(1+xn);
        
        gj2(i) = 1/(1+xn);
        gj2(ny - 1 -i) = 1/(1+xn);
        gj3(i) = (1-xn)/(1+xn);
        gj3(ny - 1 -i) = (1-xn)/(1+xn);
        
        xxn = (xnum-0.5)/xd;
        xn = coeff_pml  * xxn^3;
        
        fi1(i) = xn/2;
        fi1(nx - 2 -i) = xn/2;
        fi2(i) = 1/(1+xn);
        fi2(nx - 2 - i) = 1/(1+xn);
        fi3(i) = (1-xn)/(1+xn);
        fi2(nx - 2 - i) = (1-xn)/(1+xn);
        
        fj1(i) = xn/2;
        fj1(ny- 2 - i) = xn/2;
        fj2(i) = 1/(1+xn);
        fj2(ny - 2 - i) = 1/(1+xn);
        fj3(i) = (1-xn)/(1+xn);
        fj3(ny - 2 - i) = (1-xn)/(1+xn);
    end
T=0;
tic
Erec1=[];
Erec2=[];
vt = 2*c/1;
kx=50;
sig0=[];
dn = 4/100;
for n=1:nsteps
    n
    T=T+1;
%     dne = Dne(x,y,ne0,kx,vt,n,dt,ddx,ddy,dn);
%     net = ne0+dne;
%     omega_pe=omega_pecal(net);
    
    
%     figure
%     contour(x,y,net');
%     set(gca,'YDir','normal')
%     xlabel('R (m)')
%     ylabel('Z (m)') 
    
    % calculate the Dz
    for ii=2:nx
        for jj = 2:ny
            dz(ii,jj)=gi3(ii)*gj3(jj)*dz(ii,jj)...
            +0.5*gi2(ii)*gj2(jj)*(hy(ii,jj)-hy(ii-1,jj)...
            -(ddx/ddy)*((hx(ii,jj)-hx(ii,jj-1))));
         end
    end
   % adding the source
   xatp = 0.2;
   yatp=50*ddy;
%     w=w/4;
    sigtmp=[];
   for i = fix((xatp-w/1)/ddx):fix((xatp+w/1)/ddx)
%        is = i-fix((xatp-w/1)/ddx);
       pl = i*ddx-xatp;
%        wwg = 2*w/ddx;
%        signal = 5* cos(-2*pi*f*dt*(t0-T))*sin(is*pi/wwg);
       w0=w/1;
       signal = 5*exp(-pl^2/(w0^2))*cos(-2*pi*f*dt*(t0-T));
       sigtmp = [sigtmp signal];
%        sigx=[sigx i];
%        sig = [sig signal];
       dz(i,fix((yatp)/ddy))=dz(i,fix((yatp)/ddy))+signal;
       
   end
   sig0 = [sig0 mean(sigtmp)];
   iget1 = fix((Ax1-w/1)/ddx):fix((Ax1+w/1)/ddx);
   Erec1 = [Erec1 mean(ez(iget1,fix((yatp)/ddy)))];
   iget2 = fix((Ax2-w/1)/ddx):fix((Ax2+w/1)/ddx);
   Erec2 = [Erec2 mean(ez(iget2,fix((yatp)/ddy)))];
   
   
%    figure(100)
%    plot(sigx,sig)
    for i = 1:length(Anx)
        dz(Anx(i),Any(i)) =0;
    end

%    [dz,Anx,Any]=setLaunchAntenna(xatp,1.5*AtL,w,1.3*AtL,ddx,ddy,dz);
%    dz=setReceiverAntenna(0.23,1.5*AtL,w,1.3*AtL,ddx,ddy,dz);


   
   vc=5e-3;
   for ii=2:nx
       for jj = 2:ny
               ez(ii,jj)=dz(ii,jj)-sx(ii,jj);
               sx(ii,jj)=(1+exp(-vc.*dt)).*sx1(ii,jj)...
                        -exp(-vc.*dt).*sx2(ii,jj)...
                        +((omega_pe(ii,jj)'.^2).*dt/vc).*(1-exp(-vc.*dt)).*ez(ii,jj);
               sx2(ii,jj)=sx1(ii,jj);
               sx1(ii,jj)=sx(ii,jj);
       end
   end
   
   for ii=1:nx-1
       for jj = 1:ny-1
           curl_e = ez(ii,jj) - ez(ii,jj+1);
           ihx(ii,jj) = ihx(ii,jj) + curl_e;
           hx(ii,jj) = fj3(jj)*(ddx/ddy)*hx(ii,jj) + fj2(jj)*0.5*(ddx/ddy)*(curl_e + fi1(ii)*ihx(ii,jj));
       end
    end
   for ii=1:nx-1
       for jj = 1:ny-1
           curl_e = ez(ii,jj) - ez(ii+1,jj);
           ihy(ii,jj) = ihy(ii,jj) + curl_e;
           hy(ii,jj)=fi3(ii)*hy(ii,jj) - fi2(ii)*(0.5*curl_e +fj1(jj)*ihy(ii,jj));
       end
   end
   
   
   lev=[0.05,0.1,0.2,0.4,0.8,1.6,3.2];
   figure(2)
%    subplot(2, 1,1)
    imagesc(x,y,ez');
    colormap('jet');
    hold on
    scatter(x(Anx),y(Any),0.05,'k.')
%     figure
%     plot(1:length(Anx),x(Anx))
%     ylim([0.094 0.5])
%    contour(R(:,1),Z(1,:),ez',lev);
%    imagesc(x,y,ez);
%    contour(R(:,1),Z(1,:),ez');
% xlim([0.176 0.22])
% ylim([0 0.17])
   set(gca,'YDir','normal')
   xlabel('R (m)')
   ylabel('Z (m)')
   hold off
%        subplot(2, 1,2)
%     plot(1:length(Erec),Erec)
%    figure(3)
%     plot(1:length(Erec),Erec)
   pause(0.001)
end

t0 = (0:nsteps-1).*dt;
fs = 1/(t0(2)-t0(1));
sigr1 = Erec1;
sign = Erec1+Erec2;
[tn,phase] = cauphase(t0,Erec1,Erec2,sign,f,fs);
figure
plot(tn,phase)
toc
disp(['运行时间: ',num2str(toc)]);
ss=0;



% [x,y]=coortranser(RM,ZM,theta,r0,z0);

% figure
% imagesc(r0,z0,ne)
% set(gca,'YDir','normal')

% figure
% imagesc(x,y,ne)
% set(gca,'YDir','normal')

function dne = Dne(x,y,ne,kx,vt,ntp,dt,ddx,ddy,dn)
at = ne.*dn;
t = ntp*dt;
s=ddy;

dne=zeros(length(x),length(y));
for i=1:length(x)
    nesm=0;
   for j=1:length(y)
       nesm=0;
       xi = i*ddx;
%        xi=0;
       Atur = at(i,j)*sin(kx*(xi-vt*t));
       nesm = nesm + Atur;
%        Atur = at(i,j)*sin(kx*(xi-vt*t)/2);
%        nesm = nesm + Atur;
       dne(i,j) = nesm;
   end
end
end
function omega_pe = omega_pecal(ne)
    e = 1.602*10^(-19);
    me = 9.1096*10^(-31);
    epsilon0 = 8.85e-12;
    omega_pe = sqrt(ne*e^2/epsilon0/me);
end
function [x,y]=coortranser(ra,za,theta,r0,z0)

Lom = (za^2+ra^2)^0.5;
alpha = atan(za/ra);
beta = alpha-theta;

theta=-(pi/2-theta);
x = r0.*cos(theta)+z0.*sin(theta)+Lom*sin(beta);
y = z0.*cos(theta)+r0.*sin(theta)+Lom*cos(beta);
m=0

end
function neg=Caune(x,y)
stpot = 0.25;
idx = find(y>=stpot);

ped_pos=2.25;  %pedestal position    
h=4; %height
w=0.06; %width
slope1=0.03; %slope of the core part  
slope2=-0.1;% slope of edge part
a0=[h/2 h/2 slope1 ped_pos w slope2];
neg = zeros(length(x),length(y));
for j = idx
    dist = 2.35-y(j)+stpot;
    ne0 = MTANH(a0,dist);
    neg(:,j)=ne0;
end
% MTANH(a0,2.35 - 0.3535 + stpot)
end
function neg=Caune1(x,y,ddx,ddy,theta)
ylj = 0.2+100*ddy;
idx_ylj = fix(ylj/ddy);
neg = zeros(length(x),length(y));
ped_pos=2.25;  %pedestal position    
h=4; %height
w=0.06; %width
slope1=0.03; %slope of the core part  
slope2=-0.1;% slope of edge part
a0=[h/2 h/2 slope1 ped_pos w slope2];
m=0;
for j=idx_ylj:length(y)
    xidx = 1:length(x);
%     yidx = fix((x.*tan(theta)+0.2+100*ddx+m*ddy)/ddx);
    ne0 = MTANH(a0,2.35-m*ddy*cos(theta));
    for i1 = xidx
        j1=fix((x(i1).*tan(theta)+0.2+100*ddy+m*ddy)/ddy);
        if j1<=length(y)
            neg(i1,j1) = ne0;  
        end
    end
%     neg(yidx,xidx) = ne0;
%     figure(1)
% imagesc(x,y,neg)
% set(gca,'YDir','normal')

    m=m+1;
end
m=0;

end
function neg=Caune0(x,y)
ne = 0.*x;
idx = find(x<=2.35);
x0 = x(idx);
ped_pos=2.25;  %pedestal position    
h=4; %height
w=0.06; %width
slope1=0.03; %slope of the core part  
slope2=-0.1;% slope of edge part
 a0=[h/2 h/2 slope1 ped_pos w slope2];
ne0 = MTANH(a0,x0);
ne(idx) = ne0;
neg = meshgrid(ne,y);

end

function y=MTANH(a,x)  % MTANH function to generate the density profile for H-mode plasma
A=a(1);
B=a(2);
alpha=a(3);
x0=a(4);
w=a(5);
z=(x0-x)./w;
mtanh=((1+alpha*z).*exp(z)-(1-(0.04)*z).*exp(-z))./(exp(z)+exp(-z));
y=(A*mtanh+B)*1.e19;
end
function [Anx,Any]=setLaunchAntenna(Ax,waveGL,waveGW,AL,ddx,ddy)
    Anx = [];
    Any=[];
    waveguidthickness=2;
    setVaule=0;
    w=waveGW;
   for i = fix((Ax-w/1)/ddx)-waveguidthickness:fix((Ax-w/1)/ddx)
       for j = 1:fix(waveGL/ddy)
%            dz(i,j)=setVaule;
           Anx = [Anx i];
           Any = [Any j];
       end
   end
      for i = fix((Ax+w/1)/ddx):fix((Ax+w/1)/ddx)+waveguidthickness
       for j = 1:fix(waveGL/ddy)
%            dz(i,j)=setVaule;
           Anx = [Anx i];
           Any = [Any j];
       end
      end
   xm1 = Ax+w;
   ym1 = waveGL;
   xm2 = Ax-w;
   ym2 = waveGL;
   ap=80*pi/180;
   for j = fix(ym1/ddy):fix((ym1+AL)/ddy)
       y1 = j*ddy;
       x1 = (y1-ym1)/tan(ap)+xm1;
       ix = fix(x1/ddx):fix(x1/ddx)+waveguidthickness;
%        dz(ix,j)=setVaule;
       Anx = [Anx ix];
       Any = [Any j*ix./ix];
   end
   for j = fix(ym1/ddy):fix((ym1+AL)/ddy)
       y2 = j*ddy;
       x2 = xm2-(y2-ym2)/tan(ap);
       ix = fix(x2/ddx)-waveguidthickness:fix(x2/ddx);
%        dz(ix,j)=setVaule;
       Anx = [Anx ix];
       Any = [Any j*ix./ix];
   end

end
function [Anx,Any]=setReceiverAntenna(Ax,waveGL,waveGW,AL,ddx,ddy)
    waveguidthickness=2;
    setVaule=0;
    w=waveGW;
    Anx=[];
    Any=[];
   for i = fix((Ax-w/1)/ddx)-waveguidthickness:fix((Ax-w/1)/ddx)
       for j = 1:fix(waveGL/ddy)
           Anx = [Anx i];
           Any = [Any j];
       end
   end
      for i = fix((Ax+w/1)/ddx):fix((Ax+w/1)/ddx)+waveguidthickness
       for j = 1:fix(waveGL/ddy)
           Anx = [Anx i];
           Any = [Any j];
       end
      end
   xm1 = Ax+w;
   ym1 = waveGL;
   xm2 = Ax-w;
   ym2 = waveGL;
   ap=80*pi/180;
   for j = fix(ym1/ddy):fix((ym1+AL)/ddy)
       y1 = j*ddy;
       x1 = (y1-ym1)/tan(ap)+xm1;
       ix = fix(x1/ddx):fix(x1/ddx)+waveguidthickness;
       Anx = [Anx ix];
       Any = [Any j*ix./ix];
   end
   for j = fix(ym1/ddy):fix((ym1+AL)/ddy)
       y2 = j*ddy;
       x2 = xm2-(y2-ym2)/tan(ap);
       ix = fix(x2/ddx)-waveguidthickness:fix(x2/ddx);
       Anx = [Anx ix];
       Any = [Any j*ix./ix];
   end
end
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

