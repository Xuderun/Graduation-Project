close all; clear all;
%{
%% Plot C*
[x,y]=meshgrid(0:1:300,0:1:300);
n=6573;

while n<1000000
    load(['Ck' num2str(n) '.mat']);
    contourf(x,y,G,30);
    colorbar('northoutside');
    caxis([0 4]);
    axis equal;
    xlim([0 300]);
    ylim([0 300]);
    saveas(gcf,['C' num2str(n) '.jpg']);
    n=n+5000;
end
%}
%{
%% Plot UV
load('FinalUV.mat');
x=0:1:300;
y=300:-1:0;
fill([0 135 135 30 30 60 60 0],[300 300 270 270 30 30 0 0],'k'...
    ,[90 210 210 90],[30 30 0 0],'k',[166 300 300 240 240 270 270 166]...
    ,[300 300 0 0 30 30 270 270],'k');
hold on;
[xx,yy]=meshgrid(x,y);
streamslice(xx,yy,u,v,3);
axis equal;
xlim([1 300]);
ylim([1 300]);
saveas(gcf,'UV.jpg');
%}
%{
%% Plot Wall
load('Porous3.mat');
s=pcolor(sDist);
colormap([1 1 1;0 0 0;1 0 0;0 0 1;]);
axis equal;
xlim([1 301]);ylim([1 301]);
set(s,'EdgeColor','none');
saveas(gcf,'Wall.jpg');
%save('Porous3.mat','D','Sc','th','sDist','psi');
%}