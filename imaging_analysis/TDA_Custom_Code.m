%cd ~/Dropbox/Postdoc_Denis/cortex_gradient_patterns/TDA/cell_detection/
tic;
fileName = {
    'control1.png'
    'control2.png'
    'control3.png'
    'COUP1.png'
    'COUP2.png'
    'COUP3.png'
    'PcdhKD1.png'
    'PcdhKD2.png'
    'PcdhKD3.png'};

%fileName={
%    'control1.png'};

its = 100;
Nfiles=length(fileName);
N_TDA=its;
linS = {'-','--',':'};
MaxSlope=zeros(1,Nfiles);
ConnectedComponents=zeros(Nfiles,N_TDA);
ConnectedComponentsArea=zeros(Nfiles,N_TDA);
CCs_areas = [];
Red=zeros(Nfiles,N_TDA);
max_pos = zeros(1,Nfiles);
max_pos_mm = zeros(1,Nfiles);
num_control = 3;
N0 = zeros(3,1);
N02 = zeros(Nfiles,1);

Control_cells = [];
COUP_cells = [];
Pcdh_cells = [];

storeArea = [];

%% compute average size of initial CCs in the control phenotypes
for FileIndex = 1:num_control
    file = fileName{FileIndex};
    %control1.png COUP1.png PcdhKD1.png
    im=imread(file);
    threshold=180;
    if FileIndex==4
        threshold=140;
    end
    smooth_radius=4;
    
    % plot images
    %figure(1);
    %subplot(3,3,FileIndex), imagesc(im);
    %set(gca,'XTick',[], 'YTick', [])
    
    im=im(:,:,1);
    
    H = fspecial('gaussian',smooth_radius,smooth_radius);
    imsmooth = imfilter(im,H,'replicate');
    Imagebin=imsmooth>threshold;
    Imagebin=imclose(imopen(Imagebin,strel('disk',4)),strel('disk',4));
    CC = bwconncomp(Imagebin);
    
    A = regionprops(CC,'Area');
    %figure(1);
    %subplot(1,3,FileIndex), histogram(vertcat(A.Area),5+ceil(sqrt(length(vertcat(A.Area)))));
    storeArea = [storeArea vertcat(A.Area)'];
end
a_bar = median(storeArea); % median or average pixels per CC in WT cortices
r = sqrt(a_bar/pi); % average cell radius assuming circular cells
%% thickening of CCs
for FileIndex= 1:9
    file = fileName{FileIndex};
    %control1.png COUP1.png PcdhKD1.png
    im=imread(file);
    threshold=180;
    if FileIndex==4
        threshold=140;
    end
    smooth_radius=4;
    
    % plot images
    %figure(1);
    %subplot(3,3,FileIndex), imagesc(im);
    %set(gca,'XTick',[], 'YTick', [])
    
    im=im(:,:,1);
    
    H = fspecial('gaussian',smooth_radius,smooth_radius);
    imsmooth = imfilter(im,H,'replicate');
    Imagebin=imsmooth>threshold;
    Imagebin=imclose(imopen(Imagebin,strel('disk',4)),strel('disk',4));
    
    %     plot thresholded images
    %     figure(2);
    %     subplot(3,3,FileIndex), imagesc(Imagebin);
    %     set(gca,'XTick',[], 'YTick', [])
    %     colormap(gray);
    
    % DD=strel('disk',1);
    % DD2=strel('disk',1);
    % Imagebin=imclose(imopen(Imagebin,DD2),DD);
    for radius = 0:its
        if radius == 0
            Dilated_Image = Imagebin;
            CC = bwconncomp(Dilated_Image);
            S = regionprops(CC,'Centroid');
            ImNaN=double(Dilated_Image).*double(im);
            ImNaN(Dilated_Image==0)=NaN;
            A = regionprops(CC,'Area');
            STD(FileIndex ,radius+1)=nanstd(ImNaN(:));%/sum(vertcat(A.Area));
            ImNaN=double(1-Dilated_Image).*double(im);
            ImNaN(Dilated_Image==1)=NaN;
            %     hold on
            %     for i=1:length(S)
            %         Center=S(i).Centroid;
            %         plot(Center(1),Center(2),'or');
            %     end
            % need to normalize number of CCs by the average cell size
            %figure(44);
            CCs_areas = vertcat(A.Area); % number of cells per CC
            N02(FileIndex) = sum(CCs_areas);
            if FileIndex < 4
                Control_cells = [Control_cells CCs_areas'/a_bar];
                N0(1) = N0(1) + sum(CCs_areas);
            elseif FileIndex < 7
                COUP_cells = [COUP_cells CCs_areas'/a_bar];
                N0(2) = N0(2) + sum(CCs_areas);
            else
                Pcdh_cells = [Pcdh_cells CCs_areas'/a_bar];
                N0(3) = N0(3) + sum(CCs_areas);
            end
            %subplot(3,3,FileIndex), histogram(CCs_adjusted,3+ceil(sqrt(length(CCs_adjusted))));
            %title('Cells per CC');
            ConnectedComponents(FileIndex,radius+1) = length(S);
            ConnectedComponentsArea(FileIndex,radius+1) = sum(CCs_areas);
            Red(FileIndex,radius+1)=sum(Dilated_Image(:));
        else
            DD=strel('disk',radius);
            Dilated_Image=imdilate(Imagebin,DD);
            %         imagesc(Dilated_Image);
            %         colormap(gray);
            CC = bwconncomp(Dilated_Image);
            S = regionprops(CC,'Centroid');
            
            ImNaN=double(Dilated_Image).*double(im);
            ImNaN(Dilated_Image==0)=NaN;
            A = regionprops(CC,'Area');
            STD(FileIndex ,radius+1)=nanstd(ImNaN(:));%/sum(vertcat(A.Area));
            ImNaN=double(1-Dilated_Image).*double(im);
            ImNaN(Dilated_Image==1)=NaN;
            %     hold on
            %     for i=1:length(S)
            %         Center=S(i).Centroid;
            %         plot(Center(1),Center(2),'or');
            %     end
            CCs_areas = vertcat(A.Area);
            ConnectedComponents(FileIndex, radius+1) = length(S);
            ConnectedComponentsArea(FileIndex,radius+1) = sum(CCs_areas);
            Red(FileIndex,radius+1) = sum(Dilated_Image(:));
        end
    end
    
%     figure(19);
%     hold on
%     if FileIndex < 4
%         plot(0:its,ConnectedComponents(FileIndex,:),'b','LineWidth',2,'linestyle',linS{mod(FileIndex,3)+1});
%     elseif FileIndex < 7
%         plot(0:its,ConnectedComponents(FileIndex,:),'r','LineWidth',2,'linestyle',linS{mod(FileIndex,3)+1});
%     else
%         plot(0:its,ConnectedComponents(FileIndex,:),'g','LineWidth',2,'linestyle',linS{mod(FileIndex,3)+1});
%     end
%     title('Connected components');
%     
%     figure(24);
%     mmp = 13;
%     str = sprintf('Slope of connected components (moving mean param. = %d)',mmp);
%     title(str);
%     hold on;
%     slopes_vec = movmean(diff(ConnectedComponents(FileIndex,:)),mmp);
%     %slopes_vec = diff(movmean(ConnectedComponents(FileIndex,:),mmp));
%     [ww max_pos_mm(FileIndex)] = max(abs(slopes_vec));
%     if FileIndex < 4
%         plot(0:its-1,slopes_vec,'b','LineWidth',2,'linestyle',linS{mod(FileIndex,3)+1});
%     elseif FileIndex < 7
%         plot(0:its-1,slopes_vec,'r','LineWidth',2,'linestyle',linS{mod(FileIndex,3)+1});
%     else
%         plot(0:its-1,slopes_vec,'g','LineWidth',2,'linestyle',linS{mod(FileIndex,3)+1});
%     end
end
%% This section analyzes the number of cells per CC in the initial images (and their distributions)
figure(33);
title('Cells per CC');
subplot(3,1,1), histogram(Control_cells,5+ceil(sqrt(length(Control_cells))),'Normalization','probability'), title('Control: Cells per CC');
subplot(3,1,2), histogram(COUP_cells,5+ceil(sqrt(length(COUP_cells))),'Normalization','probability'), title('COUP: Cells per CC');
subplot(3,1,3), histogram(Pcdh_cells,5+ceil(sqrt(length(Pcdh_cells))),'Normalization','probability'), title('Pcdh KD: Cells per CC');

disp('Table 1: Initial cells per connected component');
[h1, p1] = ttest2(Control_cells,COUP_cells);
[h2, p2] = ttest2(Control_cells,Pcdh_cells);
[h3, p3] = ttest2(COUP_cells,Pcdh_cells);
tab_data = [h1 h2 h3; p1 p2 p3];
T_cells = array2table(tab_data,'VariableNames',{'ControlCOUP','ControlPcdhKD','COUPPcdhKD'})

% max_pos_mm = max_pos_mm - 1; % -1 since indices start at 1 but that is radius 1, not zero
% figure(24);
% legend(fileName,'Location','Southeast');

%% Weighted average calculation of typical distance between cells
N0 = N0/a_bar; % total number of cells per already merged CC -> weight for r
merge_data = -(diff(ConnectedComponents')'); % number of merges at each radius
dists = (0:its) + round(r);
control_merges = [N0(1) sum(merge_data(1:3,:))];
COUP_merges = [N0(2) sum(merge_data(4:6,:))];
Pcdh_merges = [N0(3) sum(merge_data(7:9,:))];

% figure();
% subplot(1,3,1), bar(control_merges), title('Control: cell merges');
% subplot(1,3,2), bar(COUP_merges), title('COUP: cell merges');
% subplot(1,3,3), bar(Pcdh_merges), title('Pcdh: cell merges');
% for i = 1:3
%     subplot(3,3,i), bar(merge_data(i,:)), title('Control: cell merges');
% end
% for i = 4:6
%     subplot(3,3,i), bar(merge_data(i,:)), title('COUP: cell merges');
% end
% for i = 7:9
%     subplot(3,3,i), bar(merge_data(i,:)), title('Pcdh: cell merges');
% end
avDistControl = sum(dists.*control_merges)/sum(control_merges);
M1 = sum(control_merges>0);
MM1 = (M1-1)/M1;
stdDistControl = sqrt( sum(control_merges.*((dists-avDistControl).^2))/( MM1*sum(control_merges) ) );

avDistCOUP = sum(dists.*COUP_merges)/sum(COUP_merges);
M2 = sum(COUP_merges>0);
MM2 = (M2-1)/M2;
stdDistCOUP = sqrt( sum(COUP_merges.*((dists-avDistCOUP).^2))/( MM2*sum(COUP_merges) ) );

avDistPcdh = sum(dists.*Pcdh_merges)/sum(Pcdh_merges);
M3 = sum(Pcdh_merges>0);
MM3 = (M3-1)/M3;
stdDistPcdh = sqrt( sum(Pcdh_merges.*((dists-avDistPcdh).^2))/( MM3*sum(Pcdh_merges) ) );

figure();
x = 1:3;
data = [avDistControl avDistCOUP avDistPcdh]';
%errhigh = ([std(max_pos_mm(1:3)) std(max_pos_mm(4:6)) std(max_pos_mm(7:9))])/sqrt(3);
%errlow  = [std(max_pos_mm(1:3)) std(max_pos_mm(4:6)) std(max_pos_mm(7:9))]/sqrt(3);

bar(x,data);
xticklabels({'Control','COUP','Pcdh KD'})
ylabel('Typical distance');
set(gca,'FontSize',20);

%% Compute standard deviation CC by CC
r_control = 8;%floor(avDistControl);
r_COUP = 8;%floor(avDistCOUP);
r_Pcdh = 8;%floor(avDistPcdh);

% dilate the images to the appropriate merge radius and compute stdev of
% each CC
control_std = [];
COUP_std = [];
Pcdh_std = [];

for i = 1:length(fileName)
    file = fileName{i};
    %control1.png COUP1.png PcdhKD1.png
    im=imread(file);
    threshold=180;
    if FileIndex==4
        threshold=140;
    end
    smooth_radius=4;
    im=im(:,:,1);
    H = fspecial('gaussian',smooth_radius,smooth_radius);
    imsmooth = imfilter(im,H,'replicate');
    Imagebin=imsmooth>threshold;
    Imagebin=imclose(imopen(Imagebin,strel('disk',4)),strel('disk',4));
    if i < 4
       radius = r_control; 
    elseif i < 7
        radius = r_COUP;
    else
        radius = r_Pcdh;
    end
    DD=strel('disk',radius);
    Dilated_Image=imdilate(Imagebin,DD);
    CC = bwconncomp(Dilated_Image);
    for j = 1:length(CC.PixelIdxList)
        if i < 4
            control_std = [control_std std(double(im(CC.PixelIdxList{j})))];
        elseif i < 7
            COUP_std = [COUP_std std(double(im(CC.PixelIdxList{j})))];
        else
            Pcdh_std = [Pcdh_std std(double(im(CC.PixelIdxList{j})))];
        end
    end
end
figure();
x = 1:3;
data = [mean(control_std); mean(COUP_std); mean(Pcdh_std)]';
errlow  = ([1.96*std(control_std)/sqrt(length(control_std)); ...
    1.96*std(COUP_std)/sqrt(length(control_std));...
    1.96*std(Pcdh_std)/sqrt(length(control_std))]);

bar(x,data);
xticklabels({'Control','COUP','Pcdh KD'})
ylabel('Heterogeneity');
set(gca,'FontSize',20);
hold on

er = errorbar(x,data,errlow, 'LineWidth', 2);
er.Color = [0 0 0];
er.LineStyle = 'none';

disp('Table 2: Stdev at typical scale');
[h1, p1] = ttest2(control_std,COUP_std);
[h2, p2] = ttest2(control_std,Pcdh_std);
[h3, p3] = ttest2(COUP_std,Pcdh_std);
tab_data = [h1 h2 h3; p1 p2 p3];
disp(array2table(tab_data,'VariableNames',{'ControlCOUP','ControlPcdhKD','COUPPcdhKD'}))
%%
N02 = N02/a_bar;
frac = N02./(ConnectedComponents(:,1)+N02);
%% fraction of connected components merged before thickening (i.e. in original image)
figure(); 
bar([mean(frac(1:3)); mean(frac(4:6)); mean(frac(7:9))]);
hold on;
errorbar([1 2 3],[mean(frac(1:3)); mean(frac(4:6)); mean(frac(7:9))],[std(frac(1:3)); std(frac(4:6));std(frac(7:9))]/sqrt(3));
title('fraction pre-merged');

disp('Table 3: Differences in CCs initially merged');
[h1,p1] = ttest2(frac(1:3),frac(4:6));
[h2,p2] = ttest2(frac(1:3),frac(7:9));
[h3,p3] = ttest2(frac(7:9),frac(4:6));
tab_data = [h1 h2 h3; p1 p2 p3];
disp(array2table(tab_data,'VariableNames',{'ControlCOUP','ControlPcdhKD','COUPPcdhKD'}))

%% fraction of connected components merged after thickening
r0 = 8;
control_stat = sum(merge_data(:,r0:end),2)./(ConnectedComponents(:,1)+N02(:));
disp('Table 4: Differences in CC merging');
[h1, p1] = ttest2(control_stat(1:3),control_stat(4:6));
[h2, p2] = ttest2(control_stat(1:3),control_stat(7:9));
[h3, p3] = ttest2(control_stat(7:9),control_stat(4:6));
tab_data = [h1 h2 h3; p1 p2 p3];
disp(array2table(tab_data,'VariableNames',{'ControlCOUP','ControlPcdhKD','COUPPcdhKD'}))

figure(); 
bar([mean(control_stat(1:3)); mean(control_stat(4:6));mean(control_stat(7:9))]);
hold on;
errorbar([1 2 3],[mean(control_stat(1:3)); mean(control_stat(4:6));mean(control_stat(7:9))],[std(control_stat(1:3)); std(control_stat(4:6));std(control_stat(7:9))]/sqrt(3));
title('fraction merged after r0');
toc;
