clear all;
close all;
addpath('ffmatlib');
[p,b,t,nv,nbe,nt,labels,regions]=ffreadmesh('4speciesKS_final.mesh');
vh=ffreaddata('4speciesKS_vh.txt');
u=ffreaddata('solution_CEnt_final.txt');
v=ffreaddata('solution_CNeo_final.txt');
% Plot the mesh
%figure(1);
%ffpdeplot3D(p,b,t,'Mesh','on',...
%    'XYZStyle','monochrome','ColorMap',parula,'ColorBar','on','BoundingBox','on');
%axis([0,10,0,5]);
% x-cordinates -> p(1,:), y-cordinates -> p(2,:), z-cordinates -> p(3,:)

%% Plot for Entorhinal cortex
figure(1);
set(1,'InvertHardcopy', 'off');
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,...
    'XYZStyle','interp','Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
axis([-3,10,0,10,0,10]);
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,...
    'XYZStyle','interp','Slice',[10 0 0],[10 0 10],[10 10 0],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,...
    'XYZStyle','interp','Slice',[-3 0 0],[-3 0 10],[-3 10 0],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,...
    'XYZStyle','interp','Slice',[10 0 0],[-3 0 10],[10 0 10],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',u,...
    'XYZStyle','interp','Slice',[-3 0 0],[-3 0 10],[10 0 10],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
colormap pink;
caxis([0.2 1.2]);
colorbar off;
%hold on;
%ffpdeplot3D(p,b,t,'BDLabels',[30,31],'XYZStyle','monochrome');
%custom_map = [
   % linspace(0.8,0,100)' linspace(0,0.8,100)' linspace(0,0,100)'];
%colormap(custom_map);
view([53.3,19.6]);
%xlabel('x-axis');
%ylabel('y-axis');
%zlabel('z-axis');
saveas(1,'3D_print_slice_ent.svg', 'svg');
%% Plot for neocortex
figure(2);
set(2,'InvertHardcopy', 'off');
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',v,...
    'XYZStyle','interp','Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
axis([-3,10,0,10,0,10]);
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',v,...
    'XYZStyle','interp','Slice',[10 0 0],[10 0 10],[10 10 0],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',v,...
    'XYZStyle','interp','Slice',[-3 0 0],[-3 0 10],[-3 10 0],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',v,...
    'XYZStyle','interp','Slice',[10 0 0],[-3 0 10],[10 0 10],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
hold on;
ffpdeplot3D(p,b,t,'VhSeq',vh,'XYZData',v,...
    'XYZStyle','interp','Slice',[-3 0 0],[-3 0 10],[10 0 10],'Boundary','on','ColorMap',parula,'Mesh','off','ColorBar','on','BoundingBox','off');
colormap pink;
caxis([0.2 1.2]);
%colorbar off;
%hold on;
%ffpdeplot3D(p,b,t,'BDLabels',[30,31],'XYZStyle','monochrome');
%custom_map = [
   % linspace(0.8,0,100)' linspace(0,0.8,100)' linspace(0,0,100)'];
%colormap(custom_map);
view([53.3,19.6]);
%xlabel('x-axis');
%ylabel('y-axis');
%zlabel('z-axis');
saveas(2,'3D_print_slice_neo_colorbar.svg', 'svg');




