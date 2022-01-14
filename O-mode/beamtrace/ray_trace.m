

mode = 'O';
IT=0;
c = 3.e8;
frequency = 30.e9;
lamda = c/frequency;
k0 = 2*pi/lamda;
ddx = lamda/20;
ddy=1*ddx;
dt=ddx/c/2;
x = 0:ddx:0.4;
y = 0:ddx:0.5;
ne=Caune(x,y);
[nx,ny] = size(ne);
R=x;
Z=y;
% [R,Z] = meshgrid(x,y);
% figure
% plot3(R,Z,ne)
% figure
% contour(R,Z,ne')
w0 = -1*10^(-2);
w1 = 1*10^(-2);
theta = 90;
RZ0=[0.2,0.08];
ray_number=10;
dtau = -0.000001;
step=1000;

zR=pi*w0^2/lamda;
k0=2*pi*frequency/c; % wave-number in vacuum
theta=theta/180*pi; % the cline angle, conversion to radian
R0=RZ0(1); % the coordinate of the beam centre
Z0=RZ0(2)
d=linspace(-w1/2,w1/2,ray_number); % the distance of the intial points from the beam axis
x0=R0+d*sin(theta);
y0=Z0+d*cos(theta);

interindex = fix(length(d)/2);

k=k0+k0*((d.^2)./(2*zR^2)); 
kx0=-k*cos(theta); %the initial kx
ky0=k*sin(theta); %the initial ky
I0=-(d.^2)./(w0^2); % the initial value of I
P0=4*d.^2/(w0^4); %P=(ki)^2=(grad I)^2

standnum=170;
beta_x0=standnum*(x0-R0)/w0^4;
beta_y0=standnum*(y0-Z0)/w0^4;
[x1,y1]=ray_tracing(x0(interindex),y0(interindex),kx0(interindex),ky0(interindex),I0(interindex),P0(interindex),beta_x0(interindex),beta_y0(interindex),dtau,step,R,Z,ne,IT,frequency,2,mode);
% calculate the cutoff layer
R_judge=min(x1)+0.05;
% for the other ray
[x,y,kx,ky,P,beta_x,beta_y]=ray_tracing(x0,y0,kx0,ky0,I0,P0,beta_x0,beta_y0,dtau,step,R,Z,ne,IT,frequency,R_judge,mode);

hold on
plot(x,y,'linewidth',1)









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
function [x,y,kx,ky,P,beta_x,beta_y]=ray_tracing(x0,y0,kx0,ky0,I0,P0,beta_x0,beta_y0,dtau,cycle_length,R,Z,ne,IT,freq,R_judge,mode)
% ray tracing calculation based on complex einkonal method for X-mode wave
% using RK4 numerical integration method
c=3*10^8; %light velocity
L=length(x0);
k0=2*pi*freq/c; %wave-number in vacuum
x=[];
y=[];
kx=[];
ky=[];
P=[];
beta_x=[];
beta_y=[];

i=1;
Dx=0.001;
Dy=0.001;

while i<cycle_length
    x=[x x0];
    y=[y y0];
    kx=[kx kx0];
    ky=[ky ky0];
    P=[P P0];
    beta_x=[beta_x beta_x0];
    beta_y=[beta_y beta_y0];
    x1=x0;
    y1=y0;
    kx1=kx0;
    ky1=ky0;
    P1=P0;
    beta_x1=beta_x0;
    beta_y1=beta_y0;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK1_x,RK1_y,RK1_kx,RK1_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    x1=x0+RK1_x*dtau/2;
    y1=y0+RK1_y*dtau/2;
    kx1=kx0+RK1_kx*dtau/2;
    ky1=ky0+RK1_ky*dtau/2;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK2_x,RK2_y,RK2_kx,RK2_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    x1=x0+RK2_x*dtau/2;
    y1=y0+RK2_y*dtau/2;
    kx1=kx0+RK2_kx*dtau/2;
    ky1=ky0+RK2_ky*dtau/2;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK3_x,RK3_y,RK3_kx,RK3_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    x1=x0+RK3_x*dtau;
    y1=y0+RK3_y*dtau;
    kx1=kx0+RK3_kx*dtau;
    ky1=ky0+RK3_ky*dtau;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK4_x,RK4_y,RK4_kx,RK4_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    RK_x=(RK1_x+2*RK2_x+2*RK3_x+RK4_x)/6;
    RK_y=(RK1_y+2*RK2_y+2*RK3_y+RK4_y)/6;
    RK_kx=(RK1_kx+2*RK2_kx+2*RK3_kx+RK4_kx)/6;
    RK_ky=(RK1_ky+2*RK2_ky+2*RK3_ky+RK4_ky)/6;
    x1_temp=x0+RK_x*dtau;
    y1_temp=y0+RK_y*dtau;
    kx1_temp=kx0+RK_kx*dtau;
    ky1_temp=ky0+RK_ky*dtau;
    
    if any(x1_temp<R_judge)
        P1_temp=zeros(1,L);
        beta_x0=zeros(1,L);
        beta_y0=zeros(1,L);
    else
        [k_ix_temp,k_iy_temp]=grad_I(x0,y0,I0,x1_temp,y1_temp,I0);
        P1_temp=k_ix_temp.^2+k_iy_temp.^2;
        [beta_x0,beta_y0]=grad_I(x0,y0,P0,x1_temp,y1_temp,P1_temp);
    end
    x0=x1_temp;
    y0=y1_temp;
    kx0=kx1_temp;
    ky0=ky1_temp;
    P0=P1_temp;
    i=i+1;
    if mod(i,200)==0
        disp(i)
    end
    if any([x0>=max(R),x0<=min(R),y0>=max(Z),y0<=min(Z)])
       break;   
    end
end

x=reshape([x x0],L,[])';
y=reshape([y y0],L,[])';
kx=reshape([kx kx0],L,[])';
ky=reshape([ky ky0],L,[])';
P=reshape([P P0],L,[])';
beta_x=reshape([beta_x beta_x0],L,[])';
beta_y=reshape([beta_y beta_y0],L,[])';
end

function [x,y,kx,ky,P,beta_x,beta_y]=ray_tracing_ver_old(x0,y0,kx0,ky0,I0,P0,beta_x0,beta_y0,dtau,cycle_length,R,Z,ne,IT,freq,R_judge)
%% ray tracing calculation based on complex einkonal method for X-mode wave
%% using RK4 numerical integration method
c=3*10^8; %light velocity
L=length(x0);
k0=2*pi*freq/c; %wave-number in vacuum
x=[];
y=[];
kx=[];
ky=[];
P=[];
beta_x=[];
beta_y=[];

i=1;
Dx=0.001;
Dy=0.001;
while i<cycle_length
    x=[x x0];
    y=[y y0];
    kx=[kx kx0];
    ky=[ky ky0];
    P=[P P0];
    beta_x=[beta_x beta_x0];
    beta_y=[beta_y beta_y0];
    x1=x0;
    y1=y0;
    kx1=kx0;
    ky1=ky0;
    P1=P0;
    beta_x1=beta_x0;
    beta_y1=beta_y0; 
    [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
    [RK1_x,RK1_y,RK1_kx,RK1_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     x1=x0+RK1_x*dtau/2;
     y1=y0+RK1_y*dtau/2;
     kx1=kx0+RK1_kx*dtau/2;
     ky1=ky0+RK1_ky*dtau/2;
     [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
     [RK2_x,RK2_y,RK2_kx,RK2_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     x1=x0+RK2_x*dtau/2;
     y1=y0+RK2_y*dtau/2;
     kx1=kx0+RK2_kx*dtau/2;
     ky1=ky0+RK2_ky*dtau/2;
     [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
     [RK3_x,RK3_y,RK3_kx,RK3_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     x1=x0+RK3_x*dtau;
     y1=y0+RK3_y*dtau;
     kx1=kx0+RK3_kx*dtau;
     ky1=ky0+RK3_ky*dtau;
     [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
     [RK4_x,RK4_y,RK4_kx,RK4_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     RK_x=(RK1_x+2*RK2_x+2*RK3_x+RK4_x)/6;
     RK_y=(RK1_y+2*RK2_y+2*RK3_y+RK4_y)/6;
     RK_kx=(RK1_kx+2*RK2_kx+2*RK3_kx+RK4_kx)/6;
     RK_ky=(RK1_ky+2*RK2_ky+2*RK3_ky+RK4_ky)/6;
     x1_temp=x0+RK_x*dtau;
     y1_temp=y0+RK_y*dtau;
     kx1_temp=kx0+RK_kx*dtau;
     ky1_temp=ky0+RK_ky*dtau;
     
     if x1_temp<R_judge
        P1_temp=zeros(1,L);
        beta_x0=zeros(1,L);
        beta_y0=zeros(1,L);
     else
     [k_ix_temp,k_iy_temp]=grad_I(x0,y0,I0,x1_temp,y1_temp,I0);
     P1_temp=k_ix_temp.^2+k_iy_temp.^2;
     [beta_x0,beta_y0]=grad_I(x0,y0,P0,x1_temp,y1_temp,P1_temp);
     end
     x0=x1_temp;
     y0=y1_temp;
     kx0=kx1_temp;
     ky0=ky1_temp;
     P0=P1_temp;
     i=i+1;
end
x=reshape([x x0],L,[])';
y=reshape([y y0],L,[])';
kx=reshape([kx kx0],L,[])';
ky=reshape([ky ky0],L,[])';
P=reshape([P P0],L,[])';
beta_x=reshape([beta_x beta_x0],L,[])';
beta_y=reshape([beta_y beta_y0],L,[])';

end



function [dI_dx,dI_dy]=grad_I(x0,y0,I0,x1,y1,I1)
  % the first point i=1 should be on the beam axis. Beam axis ray obey the
  % classical ray tracing, i.e. no need to consider the complex einkonal
   L=length(x0);
   dI_dx=zeros(1,L);
   dI_dy=zeros(1,L);
   for i=2:L;
        r1=[x0(i) y0(i)];
        r2=[x1(i) y1(i)];
        r3=[x1(i-1) y1(i-1)];
        [dI_dx(i),dI_dy(i)]=pd_3(r1,I0(i),r2,I1(i),r3,I1(i-1));
   end
end
    
function [df_dx,df_dy]=pd_3(x1,f1,x2,f2,x3,f3)
        %calculate the partial derivatives from three points
        df1=f2-f1;
        dx_1=x2(1)-x1(1); 
        dy_1=x2(2)-x1(2);
        df2=f3-f2;
        dx_2=x3(1)-x2(1);
        dy_2=x3(2)-x2(2);
        df_dx=(df2*dy_1-df1*dy_2)/(dx_2*dy_1-dx_1*dy_2);
        df_dy=(df2*dx_1-df1*dx_2)/(dy_2*dx_1-dy_1*dx_2);
end
        
function [RK_x,RK_y,RK_kx,RK_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode)
switch mode
    
    case 'X'
        f = outputs(1,:);
        df_dx = outputs(2,:);
        df_dy = outputs(3,:);
        dg_dx = outputs(4,:);
        dg_dy = outputs(5,:);
        RK_x=-2*kx1.*f;
        RK_y=-2*ky1.*f;
        RK_kx=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dx+k0^2*dg_dx-beta_x1.*f;
        RK_ky=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dy+k0^2*dg_dy-beta_y1.*f;
        
    case 'O'
        dX_dx = outputs(1,:);
        dX_dy = outputs(2,:);
        RK_x=-2*kx1;
        RK_y=-2*ky1;
        RK_kx=k0^2*dX_dx-beta_x1;
        RK_ky=k0^2*dX_dy-beta_y1;
        
end
end
   
function [RK_x,RK_y,RK_kx,RK_ky]=RK4_Xmode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy) 
        RK_x=-2*kx1.*f;
        RK_y=-2*ky1.*f;
        RK_kx=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dx+k0^2*dg_dx-beta_x1.*f;
        RK_ky=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dy+k0^2*dg_dy-beta_y1.*f;
end      
        
function [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x0,y0,Dx,Dy,R,Z,ne,IT,freq)
         e=1.602*10^(-19);
         epsilon0=8.854*10^(-12);
         me=9.1096*10^(-31);
         c=3.0*10^8;
         L=length(x0);
         x=[x0-Dx x0 x0 x0 x0+Dx];
         y=[y0 y0+Dy y0 y0-Dy y0];
         Bt=(IT*1.7)./(4086*x); % EAST toroidal magnetic field 
         index1=find(R>(min(x)-0.1)&R<(max(x)+0.1));
         index2=find(Z>(min(y)-0.1)&Z<(max(y)+0.1));
         n=interp2(R(index1),Z(index2),ne(index1,index2)',x,y,'spline'); 
         Wpe=sqrt(e^2*n./(epsilon0*me));
         Wce=e*Bt./me;
         X=(Wpe.^2)./((2*pi*freq).^2);
         Y=Wce./(2*pi*freq);
         f1=1-(X./(1-Y.^2));
         g1=((X.^2).*Y.^2)./(1-Y.^2).^2;
         f1=reshape(f1,L,[]);
         g1=reshape(g1,L,[]);
         f=f1(:,3)';
         df_dx=(((f1(:,5)-f1(:,3))./Dx+(f1(:,3)-f1(:,1))./Dx)*1/2)';
         dg_dx=(((g1(:,5)-g1(:,3))./Dx+(g1(:,3)-g1(:,1))./Dx)*1/2)';
         df_dy=(((f1(:,2)-f1(:,3))./Dy+(f1(:,3)-f1(:,4))./Dy)*1/2)';
         dg_dy=(((g1(:,2)-g1(:,3))./Dy+(g1(:,3)-g1(:,4))./Dy)*1/2)';
end
    function outputs=XO_par(x0,y0,Dx,Dy,R,Z,ne,IT,freq,mode)
        e=1.602*10^(-19);
        epsilon0=8.854*10^(-12);
        me=9.1096*10^(-31);
        c=3.0*10^8;
        L=length(x0);
        xx=[x0-Dx x0 x0 x0 x0+Dx];
        yy=[y0 y0+Dy y0 y0-Dy y0];
        Bt=(IT*1.7)./(4086*xx); % EAST toroidal magnetic field
        index1=find(R>(min(xx)-0.1)&R<(max(xx)+0.1));
        index2=find(Z>(min(yy)-0.1)&Z<(max(yy)+0.1));
        while length(index1)==1
            thres = 0.1*2;
            index1=find(R>(min(xx)-thres)&R<(max(xx)+thres));
        end
        while length(index2)==1
            thres = 0.1*2;
            index2=find(Z>(min(yy)-thres)&Z<(max(yy)+thres));
        end
        
        n=interp2(R(index1),Z(index2),ne(index1,index2)',xx,yy,'spline');
        Wpe=sqrt(e^2*n./(epsilon0*me));
        Wce=e*Bt./me;
        X=(Wpe.^2)./((2*pi*freq).^2);
        Y=Wce./(2*pi*freq);
        
        switch mode
            case 'X'
                f1=1-(X./(1-Y.^2));
                g1=((X.^2).*Y.^2)./(1-Y.^2).^2;
                f1=reshape(f1,L,[]);
                g1=reshape(g1,L,[]);
                f=f1(:,3)';
                df_dx=(((f1(:,5)-f1(:,3))./Dx+(f1(:,3)-f1(:,1))./Dx)*1/2)';
                dg_dx=(((g1(:,5)-g1(:,3))./Dx+(g1(:,3)-g1(:,1))./Dx)*1/2)';
                df_dy=(((f1(:,2)-f1(:,3))./Dy+(f1(:,3)-f1(:,4))./Dy)*1/2)';
                dg_dy=(((g1(:,2)-g1(:,3))./Dy+(g1(:,3)-g1(:,4))./Dy)*1/2)';
            case 'O'
                X = reshape(X,L,[]);
                dX_dx = (((X(:,5)-X(:,3))./Dx+(X(:,3)-X(:,1))./Dx)*1/2)';
                dX_dy = (((X(:,2)-X(:,3))./Dy+(X(:,3)-X(:,4))./Dy)*1/2)';
        end
        if mode == 'X'
            outputs=[f;df_dx;df_dy;dg_dx;dg_dy];
        elseif mode == 'O'
            outputs=[dX_dx;dX_dy];
        end
    end  
    function [f,g]=X_par_1(x0,y0,R,Z,ne,IT,freq)
          e=1.602*10^(-19);
          epsilon0=8.854*10^(-12);
          me=9.1096*10^(-31);
          c=3.0*10^8;
          Bt=(IT*1.7)./(4086*x0); % EAST toroidal magnetic field 
          n=interp2(R,Z,ne',x0,y0,'spline'); 
          Wpe=sqrt(e^2*n./(epsilon0*me));
          Wce=e*Bt./me;
          X=(Wpe.^2)./((2*pi*freq).^2);
          Y=Wce./(2*pi*freq);
          f=1-(X./(1-Y.^2));
          g=((X.^2).*Y.^2)./(1-Y.^2).^2;
    end