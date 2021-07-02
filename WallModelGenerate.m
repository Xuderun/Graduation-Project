clear all;close all;
%% Wall Model Data
lx=301;ly=301;
lx1=30;lx2=60;lx3=90;lx4=135;
lx5=165;lx6=210;lx7=240;lx8=270;
ly1=30;ly2=270;

%% Initialization of Governing Parameters
sDist=ones(ly,lx); % Wall Model Function
psi=zeros(ly,lx); % Local Porosity
th=zeros(ly,lx); % Dimensioless Concentration
r=sqrt(8); % Growth Radius
D=ones(ly,lx); % Dimensionless D*
Sc=ones(ly,lx); % Schimdt number

%% Initialization of Random System
shuflim=150; 
p1=0;p2=0;
l=0;m=0;

%% Wall Model Area
Area1=(ly-ly2)*(lx4+lx-lx5); % Left Wall
Area2=(ly2-ly1+1)*(lx1+lx-lx8); % Middle
Area3=(ly1)*(lx2+lx6-lx3+lx-lx7); % Right Wall
TRange=Area1+Area2+Area3; % Total Area for particle growth 
p1Range=TRange*0.3; % Phase 2 Material Volume Fraction 
p2Range=TRange*0.1; % Phase 3 Material Volume Fraction
tic;

%% Elimination Part (Inner Room and ventilation port)
for i=1:ly
    for j=1:lx
        if i>=ly2 && i<=ly && j>=lx4 && j<lx5
            sDist(i,j)=0;
        elseif i>=ly1 && i<ly2 && j>=lx1 && j<lx8
            sDist(i,j)=0;
        elseif i>0 && i<ly1
            if (j>=lx2 && j<lx3) || (j>=lx6 && j<lx7)
                sDist(i,j)=0;
            end
        end
    end
end
%% Wall Model Generation
while p1<p1Range || p2<p2Range
    % X,Y Random Position
    x=randperm(lx,shuflim);
    y=randperm(ly,shuflim);
    px=randperm(lx,shuflim);
    py=randperm(ly,shuflim);
    %% Phase 2 Generation
    for i=1:shuflim
        touch= false; % Overlapped boolean
       for j=-ceil(r)+1:ceil(r)
           if j<0
                %% X,Y Boundary
                if x(1,i) <3
                    x1Low =1;
                    x1Up=x(1,i)-j;
                elseif x(1,i)>lx-2 && x(1,i)<=lx
                    x1Low=x(1,i)+j;
                    x1Up=lx;
                else
                    x1Low=x(1,i)+j;
                    x1Up=x(1,i)-j;
                end
                if y(1,i) <3
                    y1Low =1;
                    y1Up=y(1,i)-j;
                elseif y(1,i)>lx-2 && y(1,i)<=lx
                    y1Low=y(1,i)+j;
                    y1Up=ly;
                else
                    y1Low=y(1,i)+j;
                    y1Up=y(1,i)-j;
                end
                %% Determine overlapped
                if y1Low>ly2 && y1Up <=ly
                    if (x1Low>0 && x1Up <=lx4) || (x1Low>lx5 && x1Up <=lx)
                        l=sum(sum(sDist(y1Low:y1Up,x1Low:x1Up)==2));
                        m=sum(sum(sDist(y1Low:y1Up,x1Low:x1Up)==3));
                    end
                elseif y1Low>ly1 && y1Up <=ly2
                    if (x1Low>0 && x1Up<=lx1) || (x1Low>lx8 && x1Up<=lx)
                        l=sum(sum(sDist(y1Low:y1Up,x1Low:x1Up)==2));
                        m=sum(sum(sDist(y1Low:y1Up,x1Low:x1Up)==3));
                    end
                elseif y1Low>0 && y1Up<=ly1
                    if (x1Low>0 && x1Up<=lx2) || (x1Low>lx3 && x1Up<lx6) ...
                            || (x1Low>lx7 && x1Up<=lx)
                        l=sum(sum(sDist(y1Low:y1Up,x1Low:x1Up)==2));
                        m=sum(sum(sDist(y1Low:y1Up,x1Low:x1Up)==3));
                    end
                end

				%% Result of overlapped
                if l==0 && m==0
                    touch=false;
                else
                    touch=true;
                    break
                end
            end
          for k=-ceil(r):ceil(r)
              x1=x(1,i)+j;
              y1=y(1,i)+k;
              p1Radius=sqrt((x1-x(1,i)).^2+(y1-y(1,i)).^2);
              % Grow Material
              if touch == false && p1Radius<=r
                  if y1>0 && y1<ly1
                      if ((x1>0 && x1<lx2)||(x1>=lx3 && x1<lx6)||...
                              (x1>lx7 && x1<=lx))
                          sDist(y1,x1)=2;
                      end
                  elseif y1>=ly1 && y1<ly2
                      if ((x1>0 && x1<lx1)||(x1>=lx8 && x1<=lx))
                          sDist(y1,x1)=2;
                      end
                  elseif y1>=ly2 && y1<=ly
                      if ((x1>0 && x1<lx4) || (x1>=lx5&& x1<=lx))
                          sDist(y1,x1)=2;
                      end
                  end
              end
			  %% Sum of the Phase 2 Material Volume
              p1=sum(sum(sDist==2));
              if p1>=p1Range
                  break
              end
          end
       end
    end
    %% P2 Generate
     for i=1:shuflim
        touch2=false; % Overlapped boolean
        for j=-ceil(r)+1:ceil(r)
            if j<0
                %% X,Y Boundary
                if px(1,i) <3
                    x2Low =1;
                    x2Up=px(1,i)-j;
                elseif px(1,i)>lx-2 && px(1,i)<=lx
                    x2Low=px(1,i)+j;
                    x2Up=lx;
                else
                    x2Low=px(1,i)+j;
                    x2Up=px(1,i)-j;
                end
                if py(1,i) <3
                    y2Low =1;
                    y2Up=py(1,i)-j;
                elseif py(1,i)>lx-2 && py(1,i)<=lx
                    y2Low=y(1,i)+j;
                    y2Up=ly;
                else
                    y2Low=py(1,i)+j;
                    y2Up=py(1,i)-j;
                end
                %% Determine Overlapped
                if y2Low>ly2 && y2Up <=ly
                    if (x2Low>0 && x2Up <=lx4) || (x2Low>lx5 && x2Up <=lx)
                        l=sum(sum(sDist(y2Low:y2Up,x2Low:x2Up)==2));
                        m=sum(sum(sDist(y2Low:y2Up,x2Low:x2Up)==3));
                    end
                elseif y2Low>ly1 && y2Up <=ly2
                    if (x2Low>0 && x2Up<=lx1) || (x2Low>lx8 && x2Up<=lx)
                        l=sum(sum(sDist(y2Low:y2Up,x2Low:x2Up)==2));
                        m=sum(sum(sDist(y2Low:y2Up,x2Low:x2Up)==3));
                    end
                elseif y2Low>0 && y2Up<=ly1
                    if (x2Low>0 && x2Up<=lx2) || (x2Low>lx3 && x2Up<=lx6) ...
                            || (x2Low>lx7 && x2Up<=lx)
                        l=sum(sum(sDist(y2Low:y2Up,x2Low:x2Up)==2));
                        m=sum(sum(sDist(y2Low:y2Up,x2Low:x2Up)==3));
                    end
                end
				
				%% Result of overlapped
                if l==0 && m==0
                    touch2=false;
                else
                    touch2=true;
                    break
                end
            end
            for k=-ceil(r):ceil(r)
                x2=px(1,i)+j;
                y2=py(1,i)+k;
                p2Radius=sqrt((x2-px(1,i)).^2+(y2-py(1,i)).^2);
                % Phase 3 Material Growth
                if p2Radius<=r && touch2==false
                    if y2>0 && y2<ly1
                        if ((x2>0 && x2<lx2)||(x2>=lx3 && x2<lx6)||...
                                (x2>lx7 && x2<=lx))
                            sDist(y2,x2)=3;
                        end
                    elseif y2>=ly1 && y2<ly2
                        if ((x2>0 && x2<lx1)||(x2>=lx8 && x2<=lx))
                            sDist(y2,x2)=3;
                        end
                    elseif y2>=ly2 && y2<=ly
                        if ((x2>0 && x2<lx4) || (x2>=lx5&& x2<=lx))
                            sDist(y2,x2)=3;
                        end
                    end
                end

				%% Sum of Phase 3 Material Volume
                p2=sum(sum(sDist==3));
                if p2>=p2Range
                    break
                end
            end
            if p2>=p2Range 
                break
            end
        end
    end
	%% Sum of Phase 1 Material Volume
    p0=sum(sum(sDist==1));
end
%% Local Porosity Generate(determinination of solid wall)
for i=1:ly
    for j=1:lx
		% If is in the enclosure ,psi=1
        if sDist(i,j)==0
            psi(i,j)=1;
        end
    end
end
%% C* Concentration, D* and Schimdt number holds in the wall
for i=1:ly
    for j=1:lx
        if sDist(i,j)==1
            Sc(i,j)=10;
        elseif sDist(i,j)==2
            Sc(i,j)=1000;
            th(i,j)=100;
        elseif sDist(i,j)==3
            Sc(i,j)=100;
        end
        D(i,j)=1/Sc(i,j);
    end
end

%% Plot and Save settings
save('Porous5.mat','sDist','psi','th','D','Sc');
s=pcolor(sDist);
colormap([1 1 1;0 0 0;1 0 0;0 0 1;]);
axis equal;
xlim([1 301]);ylim([1 301]);
set(s,'EdgeColor','none');