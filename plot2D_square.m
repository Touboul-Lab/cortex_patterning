%% Read in the mesh and solution data
close all
L = 10;
p1 = 0.48;
p2 = 1 - p1;
times = 1;
run=sprintf('run%d',0);

addpath('ffmatlib');
%Load the mesh
[p,b,t,nv,nbe,nt,labels]=ffreadmesh(strcat('4speciesKS', num2str(times),run,'.msh'));
%Load the finite element space connectivity
vh=ffreaddata(strcat('4speciesKS_vh', num2str(times),run,'.txt'));
%Load scalar data
u=ffreaddata(strcat('solution_Ent', num2str(times),run,'.txt'));
v=ffreaddata(strcat('solution_Neo', num2str(times),run,'.txt'));
%% Plots
figure(1);
ffpdeplot(p,b,t,'VhSeq',vh,'XYData',u,'Mesh','off','Boundary','off', ...
    'XLim',[0 L],'YLim',[0 L]);
xticks([0 20 40]);
yticks([0 20 40]);
colormap pink;
set(gca,'fontsize',20);
caxis([0 1.35]);
%colorbar('off');
hold on;
plot([p1*L p1*L],[0 L],'LineWidth',1.5,'Color','r');
hold on;
plot([p2*L p2*L],[0 L],'LineWidth',1.5,'Color','r');

figure(2);
ffpdeplot(p,b,t,'VhSeq',vh,'XYData',v,'Mesh','off','Boundary','off', ...
    'XLim',[0 L],'YLim',[0 L]);
xticks([0 20 40]);
yticks([0 20 40]);
colormap pink;
set(gca,'fontsize',20);
caxis([0 1.35]);
%colorbar('off');
hold on;
plot([p1*L p1*L],[0 L],'LineWidth',1.5,'Color','r');
hold on;
plot([p2*L p2*L],[0 L],'LineWidth',1.5,'Color','r');

