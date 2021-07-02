close all; clear all;clc;
load('Porous5.mat'); % Load Wall Model
%% Initialization Parameters
lx=300;
ly=300;
uin=0.0;
vin=-0.003; % Average Velocity of inlet
Tin=0.0;
d=1000; % Density
Rein=100; % Scaled Reynolds Number
Poro=0.3; % Porosity of wall
Cs=1/sqrt(3); % Dimensionless sound of speed
Dair=1.4e-5; % Mass Diffusivity of Air
U=1.7; % Macroscopic velocity
c=1e-2; % Inlet Velocity
cD=9e-4; % Dimensionless mass diffusivity
Ma=Dair/(0.01*cD); % Mach Number before open
Pe=1./(D*Ma); % Peclet Number
MaC=U/(Cs*c); % Mach Number after open
Re=Rein/(Cs*MaC); % Mesoscopic Reynolds Number

tm=1000000;
dt=1;
dx=1;
dy=dx;

%% Tracker initialization
num=tm/1000;
fcenter=zeros(1,num);% Fluid Concentration
scenter=zeros(1,num);% Solid Concentration
outletC=zeros(1,num);% Outlet Concentration
Gcenter=zeros(1,num);% Center Concentration

%% Number of nodes
nx = lx/dx+1;
ny = ly/dy+1;

%% Initialization of equation
F=d*ones(ny,nx); % Momentum Equation
G=th; % Scalar Concentration(Mass Transfer)
Gref=zeros(1,tm+1); % Reference Concentration
Open= false; % Ventilation systems activation
u=zeros(ny,nx);
v=zeros(ny,nx);

%% LBM Setting 
Tf= 1/(Cs.^2 *Re)+0.5; % Relaxation Time(Momentum)
Tg= 1./(Cs.^2 *Pe)+0.5; % Relaxation Time(Mass Transfer)
Wm=1/Tf;
Ws=1./Tg;
wk=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9]; % Weight Factor
Guvk=3*cat(3,u,v,-u,-v,u+v,-u+v,-u-v,u-v,u*0);
Fuvk=3*cat(3,u,v,-u,-v,u+v,-u+v,-u-v,u-v,u*0);
fk=zeros([size(F) 9]);
feq=zeros([size(F) 9]);
gk=zeros([size(F) 9]);
geq=zeros([size(F) 9]);
fy=zeros([size(F) 9]); % Dimensionless Force term

for i=1:9
	fk(:,:,i)=wk(i)*F;
	gk(:,:,i)=wk(i)*G;
end
% Unit Vector(Discrete velocity)
ck=[1 0 -1 0 1 -1 -1 1 0;0 1 0 -1 1 1 -1 -1 0];

%% Computation settings
n=0;m=1;
tic;
% For Streaming Step
x1=1.5:nx-0.5;
y1=1.5:ny-0.5;
[x2,y2]=meshgrid(1:10:301,301:-10:1);
[x3,y3]=meshgrid(0:dx:lx,0:dy:ly);

while n<tm
	%% Center Concentration determination
    if Gref(1,m)>=1
        Open=true;
    else
        Open=false;
    end
    if Open == true
        Pe=Re*Sc;
        Tg=1./(Cs.^2 *Pe)+0.5;
        Ws=1./Tg;
       
        for i=1:9
            Fuvk(:,:,i)=3*(ck(1,i)*u+ck(2,i)*v)+...
                3/2*(3*(ck(1,i)*u+ck(2,i)*v).^2-(u.^2+v.^2));
            feq(:,:,i)=wk(i)*F.*(1+Fuvk(:,:,i));
            % Collision Step
            fk(:,:,i)=fk(:,:,i)*(1-Wm)+Wm*feq(:,:,i)+wk(i)*fy(:,:,i);
            % Streaming Step
            fk(round(-0.5*ck(2,i)+y1),round(0.5*ck(1,i)+x1),i)=fk(round(0.5*ck(2,i)+y1),round(-0.5*ck(1,i)+x1),i);
       end
       %% Flow B.C.
       % North inlet B.C.
       dn=(fk(1,136:165,9)+fk(1,136:165,1)+fk(1,136:165,3)+...
           2*(fk(1,136:165,2)+fk(1,136:165,6)+fk(1,136:165,5)))./(1+v(1,136:165));
       fk(1,136:165,4)=fk(1,136:165,2)-(2/3)*dn.*v(1,136:165);
       fk(1,136:165,7)=fk(1,136:165,5)+0.5*(fk(1,136:165,1)-fk(1,136:165,3))...
           -dn.*v(1,136:165)/6-0.5*dn.*u(1,136:165);
       fk(1,136:165,8)=fk(1,136:165,6)+0.5*(fk(1,136:165,3)-fk(1,136:165,1))...
           -dn.*v(1,136:165)/6+0.5*dn.*u(1,136:165);
       
       % Inlet Bounce Back
       fk(1:31,136,1)=fk(1:31,136,3);
       fk(1:31,136,5)=fk(1:31,136,7);
       fk(1:31,136,8)=fk(1:31,136,6);
       fk(1:31,165,3)=fk(1:31,165,1);
       fk(1:31,165,7)=fk(1:31,165,5);
       fk(1:31,165,6)=fk(1:31,165,8);
       % North Bounce Back
       fk(31,31:135,4)=fk(31,31:135,2);
       fk(31,31:135,7)=fk(31,31:135,5);
       fk(31,31:135,8)=fk(31,31:135,6);    
       fk(31,166:end-31,4)=fk(31,166:end-31,2);
       fk(31,166:end-31,7)=fk(31,166:end-31,5);
       fk(31,166:end-31,8)=fk(31,166:end-31,6);
       
       % West Bounce Back
       fk(31:end-31,31,1)=fk(31:end-31,31,3);
       fk(31:end-31,31,5)=fk(31:end-31,31,7);
       fk(31:end-31,31,8)=fk(31:end-31,31,6);
       
       % East Bounce Back
       fk(31:end-31,end-31,3)=fk(31:end-31,end-31,1);
       fk(31:end-31,end-31,7)=fk(31:end-31,end-31,5);
       fk(31:end-31,end-31,6)=fk(31:end-31,end-31,8);
       
       % South Bounce Back
       fk(end-31,31:60,2)=fk(end-31,31:60,4);
       fk(end-31,31:60,5)=fk(end-31,31:60,7);
       fk(end-31,31:60,6)=fk(end-31,31:60,8);
       fk(end-31,91:210,2)=fk(end-31,91:210,4);
       fk(end-31,91:210,5)=fk(end-31,91:210,7); 
       fk(end-31,91:210,6)=fk(end-31,91:210,8); 
       fk(end-31,241:end-30,2)=fk(end-31,241:end-30,4);
       fk(end-31,241:end-30,5)=fk(end-31,241:end-30,7); 
       fk(end-31,241:end-30,6)=fk(end-31,241:end-30,8);
       % South Outlet Bounce Back
       fk(270:end,61,1)=fk(270:end,61,3);
       fk(270:end,61,5)=fk(270:end,61,7);
       fk(270:end,61,8)=fk(270:end,61,6);
       fk(270:end,90,3)=fk(270:end,90,1);
       fk(270:end,90,7)=fk(270:end,90,5);        
       fk(270:end,90,6)=fk(270:end,90,8);        
       fk(270:end,211,1)=fk(270:end,211,3);        
       fk(270:end,211,5)=fk(270:end,211,7);        
       fk(270:end,211,8)=fk(270:end,211,6);        
       fk(270:end,240,3)=fk(270:end,240,1);        
       fk(270:end,240,7)=fk(270:end,240,5);        
       fk(270:end,240,6)=fk(270:end,240,8);        
       % South Outlet        
       fk(end,61:90,2)=2*fk(end-1,61:90,2)-fk(end-2,61:90,2);        
       fk(end,61:90,5)=2*fk(end-1,61:90,5)-fk(end-2,61:90,5);        
       fk(end,61:90,6)=2*fk(end-1,61:90,6)-fk(end-2,61:90,6);        
       fk(end,211:240,2)=2*fk(end-1,211:240,2)-fk(end-2,211:240,2);        
       fk(end,211:240,5)=2*fk(end-1,211:240,5)-fk(end-2,211:240,5);        
       fk(end,211:240,6)=2*fk(end-1,211:240,6)-fk(end-2,211:240,6);
       
       F=sum(fk,3);
       u =zeros(size(u));
       v =zeros(size(v));

        for i =1:9
            u=u+ck(1,i)*fk(:,:,i);		
            v=v+ck(2,i)*fk(:,:,i);        
        end
        u=u./F;
        v=v./F;
        for i=1:9
            fy(:,:,i)=3*Poro*ck(2,i)*v.*(fk(:,:,i)./F);
        end
        % Wall velocity
        u(1:30,1:135)=0;u(1:30,166:end)=0;
        u(31:end,1:30)=0;u(31:end,end-30:end)=0;
        u(end-30:end,31:60)=0;u(end-30:end,91:210)=0;    
        u(end-30:end,241:end-31)=0;	
        v(1:30,1:135)=0;v(1:30,166:end)=0;    
        v(31:end,1:30)=0;v(31:end,end-30:end)=0;    
        v(end-30:end,31:60)=0;v(end-30:end,91:210)=0;    
        v(end-30:end,241:end-31)=0;
        % Inlet and Outlet velocity    
        v(1,136:165)=vin;    
        v(end,61:90)=v(end-1,61:90);v(end,211:240)=v(end-1,211:240);    
    end
    
    for i=1:9
        Guvk(:,:,i)=3*(ck(1,i)*u+ck(2,i)*v)+...
            3/2*(3*(ck(1,i)*u+ck(2,i)*v).^2-(u.^2+v.^2));
        geq(:,:,i)=wk(i)*G.*(1+Guvk(:,:,i));
        % Collision Step 
        gk(:,:,i)=gk(:,:,i).*(1-Ws)+Ws.*geq(:,:,i);
        % Streaming Step
        gk(round(-0.5*ck(2,i)+y1),round(0.5*ck(1,i)+x1),i)=gk(round(0.5*ck(2,i)+y1),round(-0.5*ck(1,i)+x1),i);
    end
    
    %% Boundary Conditions
    if Open == true
        % North Inlet B.C. (Adiabatic)    
        gk(end,136:165,2)=(wk(4)+wk(2))*G(end,136:165)-gk(end,136:165,4);   
        gk(end,136:165,5)=(wk(7)+wk(5))*G(end,136:165)-gk(end,136:165,7);
        gk(end,136:165,6)=(wk(8)+wk(6))*G(end,136:165)-gk(end,136:165,8);
    else
        for i=1:9
            gk(end,136:165,i)=gk(end-1,136:165,i);
        end
    end
    % Adiabatic B.C.
    for i=1:9
        % North
        gk(end,1:135,i)=gk(end-1,1:135,i);
        gk(end,166:end,i)=gk(end-1,166:end,i);
        % West
        gk(:,1,i)=gk(:,2,i);
        % East
        gk(:,end,i)=gk(:,end-1,i);
        % South
        gk(1,1:60,i)=gk(2,1:60,i);
        gk(1,91:210,i)=gk(2,91:210,i);
        gk(1,241:end,i)=gk(2,241:end,i);
    end

    % South B.C. (Outlet)
    if Open == true
        for i=1:9
            gk(2,61:90,i)=gk(1,61:90,i);
            gk(2,211:240,i)=gk(1,211:240,i);
        end
    else
        for i=1:9
            gk(1,61:90,i)=gk(2,61:90,i);
            gk(1,211:240,i)=gk(2,211:240,i);
        end
    end
    
	G=sum(gk,3);
	
	G(1,:)=G(2,:);G(:,1)=G(:,2);
    G(:,end)=G(:,end-1);G(end,:)=G(end-1,:);
    if Open== true
       G(end,136:165)=Tin; 
       G(2,61:90)=G(1,61:90);
       G(2,211:240)=G(1,211:240);
    else
       G(end,136:165)=G(end-1,136:165);
       G(1,61:90)=G(2,61:90);
       G(1,211:240)=G(2,211:240);
    end
    
	%% Tracking system per 1000 timestep
    if mod(n,1000)==0
        fcenter(1,ceil(m/1000))=sum(sum(G.*psi))/sum(sum(psi));
        scenter(1,ceil(m/1000))=sum(sum(G.*(1-psi)))/sum(sum(1-psi));
        outletC(1,ceil(m/1000))=mean(mean(G(1,61:90)+G(1,211:240)));
        Gcenter(1,ceil(m/1000))=G(150,150);
        save('MeanC.mat','fcenter','scenter','outletC','Gcenter');
    end
    
    %% Output file
    if n<1500
        if mod(n,50)==0
            save('CMiscPe.mat','Pe');
            save(['Ck' num2str(n) '.mat'],'G','geq','gk','Ws');
            save(['Uv' num2str(n) '.mat'],'u','v','F','feq','fk','Wm');
        elseif n==1484
            save(['Ck' num2str(n) '.mat'],'G','geq','gk','Ws');
            save(['Uv' num2str(n) '.mat'],'u','v','F','feq','fk','Wm');
        end
    else
       if mod(n,5000)==1484
           save(['LBM' num2str(n) '.mat'],'n','fy','Pe','Re');
           save(['Ck' num2str(n) '.mat'],'G','geq','gk','Ws');
           save(['Uv' num2str(n) '.mat'],'u','v','F','feq','fk','Wm');
       end 
    end
    
	n=n+dt;
    m=m+1;
    Gref(1,m)=G(150,150);
    
    %% Check Plot
    %{
    if mod(n,10)==1
        quiver(x2,y2,u(1:10:301,1:10:301),v(1:10:301,1:10:301),3);
        title(sprintf('Number of time step is %d .',n));
        axis equal;
        pause (0.1)
    end
    %}
    if mod(n,100)==1       
            contourf(x3,y3,G,30);      
            title(sprintf('TimeStep: %d.',n));       
            axis equal;       
            xlim([0 300]);       
            ylim([0 300]);
            caxis([0 4]);
            pause(0.1);    
    end   
    
end
save('LBMReference.mat','Gref','tm','n','Pe','Re');
save('FinalC.mat','G','geq','gk','Ws');
save('FinalUV.mat','u','v','F','feq','fk','fy','Wm');
%% Plot Setting
Ct=toc;
x=0:dx:lx;
y=ly:-dy:0;
[xx,yy]=meshgrid(x,y);
s=contourf(x3,y3,G,30);
colorbar('eastoutside');
caxis([0 4]);
axis equal;
xlim([0 300]);
ylim([0 300]);
xlabel('X');
ylabel('Y');
title(sprintf('Timestep: %d\nTotal time: %f',n,Ct));