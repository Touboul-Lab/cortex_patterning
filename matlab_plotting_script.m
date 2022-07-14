%% Read in the mesh and solution data
close all
%p1_vect= [0.15,0.25,0.333,0.45,0.4];
p1_vect = 0.15;
for i=0
    times = [2000];
    L = 40;
    run=sprintf('run%d',i);
    
    p1=p1_vect;
    p2=1-p1_vect;
    
    addpath('ffmatlib');
    %Load the mesh
    [p,b,t,nv,nbe,nt,labels]=ffreadmesh(strcat('4speciesKS', num2str(times),run,'.msh'));
    %Load the finite element space connectivity
    vh=ffreaddata(strcat('4speciesKS_vh', num2str(times),run,'.txt'));
    %Load scalar data
    u=ffreaddata(strcat('solution_Ent', num2str(times),run,'.txt'));
    v=ffreaddata(strcat('solution_Neo', num2str(times),run,'.txt'));
    %% Plots
    figure(i+1);
    [p_sorted, p_order] = sort(p(1,:));
    L=max(p_sorted);
    u_sorted = u(p_order,:);
    v_sorted = v(p_order,:);
    
    colors = v_sorted' > u_sorted';
    y_lim = 0:0.001:1.4; % fix range
    shade_map = repmat(colors, length(y_lim), 1);
    
    imagesc([0 L], [0 max(y_lim)],shade_map);
    shading interp;
    custom_map = [1 1 1
        235/255 235/255 235/255];
    colormap(custom_map);
    set(gca,'YDir','normal');
    
    hold on;
    plot(p_sorted,u_sorted,'LineWidth',4);
     hold on;
    plot(p_sorted,v_sorted,'LineWidth',4);
    hold on;
    plot([p1*L p1*L],[0 1.4],'LineWidth',4,'Color','r');
    hold on;
    plot([p2*L p2*L],[0 1.4],'LineWidth',4,'Color','r');
    ylim([0 1.4]);
    
    set(gca,'FontSize',24);
    
    %if i == 4 || i == 0 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6
    saveas(i+1,run,'svg');
    %end
    %title(sprintf('distance=%.2f',(p2-p1)*L))
end
