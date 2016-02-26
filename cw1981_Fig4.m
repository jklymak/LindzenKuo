% CW1981: Fig. 4
clear
% model parameters
I=256; J=128;
omega=1.45e-4; f=1.06e-4; % tidal frequency and Coriolis parameter
Nsq=0.01^2; % buoyance
sig=sqrt((omega^2-f^2)/(Nsq-omega^2)); % characteristic slope for no thermocline tilting

x=linspace(0e3,80e3,I); dx=x(2)-x(1);
 
% geometry/topography
H=225-75*tanh(0.1e-3*(x-30e3)); %H=300+x*0; 
                                %H=225-75*tanh(0.04e-3*(x-40e3)); %H=300+x*0; 
if 0
    H = H*0.+225;
    ind = find((x>30e3)&(x<50e3));
    H(ind)=linspace(225,225-75,length(ind));
    H(ind(end):end)=225-75;    
end 
Hx=gradient(H,x); Hxx=gradient(Hx,x);


P1=zeros(J,J); P2=P1;
for j=2:J-1
    P1(j,j-1)=-1; P1(j,j+1)=1;
    P2(j,j-1)=1; P2(j,j)=-2; P2(j,j+1)=1;
end
P1(1,2)=1; P1(J,J-1)=-1;
P2(1,1)=-2;  P2(1,2)=1; P2(J,J-1)=1; P2(J,J)=-2;

% On the incident boundary to find alpha_1 and beta_1
dz=H(1)/(J+1); z=(-H(1)+dz:dz:-dz)/H(1); dz=dz/H(1);
for j=1:J
    Z(j,j)=z(j);
end

A=-(1/dx*eye(J)+Hx(1)/H(1)/2/dz*Z*P1);
B=1/dx*eye(J);

for m=1:J               
    for n=1:J
        kn=n*pi/H(1)*sig;
        E1(m,n)=-sqrt(2/H(1))*sin(n*pi*z(m))*exp(-1i*kn*x(1));
        E1d(m,n)=E1(m,n)*1i*kn;
    end
end

% normal modes of the incident wave
amp=1;
psin=-amp*sin(pi*z');
k0=pi/H(1)*sig; psinkn=-amp*sin(pi*z')*k0;

alp(:,:,1)=-inv(A+E1d*inv(E1))*B;
bta(:,1) =  inv(A+E1d*inv(E1))*(1i*psinkn+E1d*inv(E1)*psin)*exp(1i*k0*x(1));

% find all alpha and beta
for i=2:I-1
    dz=H(i)/(J+1); z=(-H(i)+dz:dz:-dz)/H(i); dz=dz/H(i);
    for j=1:J
        Z(j,j)=z(j);
    end
    
    G2=-2*Hx(i)/H(i)*Z;
    G3=-(Hxx(i)*H(i)-2*Hx(i)^2)/H(i)^2*Z;
    G4=(Hx(i)^2*Z*Z-sig^2*eye(J))/H(i)^2;
    A=1/dx^2*eye(J)-1/4/dx/dz*G2*P1; 
    B=-2/dx^2*eye(J)+1/2/dz*G3*P1+1/dz^2*G4*P2;
    C=1/dx^2*eye(J)+1/4/dx/dz*G2*P1;
    
    temp1=inv(A*alp(:,:,i-1)+B);
    alp(:,:,i)=-temp1*C;
    bta(:,i)=-temp1*A*bta(:,i-1);
end

% E2 on the transmitted boundary
dz=H(I)/(J+1); z=(-H(I)+dz:dz:-dz)/H(I); dz=dz/H(I);
for j=1:J
    Z(j,j)=z(j);
end
for m=1:J
    for n=1:J
        kn=n*pi/H(I)*sig; 
        E2(m,n)=-sqrt(2/H(I))*sin(n*pi*z(m))*exp(1i*kn*x(I));
        E2d(m,n)=E2(m,n)*1i*kn;
    end
end

phi(:,I)=inv(eye(J)-Hx(I)/H(I)/2/dz*dx*Z*P1-alp(:,:,I-1)-dx*E2d*inv(E2))*bta(:,I-1);
for i=I-1:-1:1
    phi(:,i)=alp(:,:,i)*phi(:,i+1)+bta(:,i);
end

% coordinate transform for plotting
phi1=zeros(J+2,I); phi1(2:J+1,:)=phi;
[xxx,zzz]=meshgrid(x,1:J+2);
u=xxx*0; w=u;
for i=1:length(x)
    dzzz=H(i)/(J+1);
    zzz(:,i)=-H(i):dzzz:0;
    u(:,i)=-gradient(phi1(:,i),zzz(:,i));
end
w(:,1)=(phi1(:,2)-phi1(:,1))/dx-Hx(1)./H(1)*zzz(:,1).*u(:,1);
w(:,I)=(phi1(:,I)-phi1(:,I-1))/dx-Hx(I)./H(I)*zzz(:,I).*u(:,I);
for i=2:I-1
    w(:,i)=(phi1(:,i+1)-phi1(:,i-1))/dx/2-Hx(i)./H(i)*zzz(:,i).*u(:,i);
end


save -V7 cw1981_Fig4.mat xxx zzz u w 
figure
subplot(2,1,1)
% pcolor(xxx/1e3,zzz,abs(phi1)); shading flat; colorbar
contour(xxx/1e3,zzz,abs(phi1),0:0.1:1.2); colorbar
hold on; fill([x/1e3 x(end)/1e3],[-H -500],[0.7 0.7 0.7])
axis([x(1)/1e3 x(end)/1e3 -300 0]); caxis([0 1])
xlabel('x (km)','FontSize',12); ylabel('z (m)','FontSize',12)
title('\psi (m^2 s^-^1)','FontSize',15)
set(gca,'FontSize',12)
subplot(2,1,2)
pcolor(xxx/1e3,zzz,real(u)); shading flat; colorbar
hold on; fill([x/1e3 x(end)/1e3],[-H -500],[0.7 0.7 0.7])
axis([x(1)/1e3 x(end)/1e3 -300 0]); caxis([-0.02 0.02])
xlabel('x (km)','FontSize',12); ylabel('z (m)','FontSize',12)
title('u (m s^-^1)','FontSize',15)
set(gca,'FontSize',12)