%%
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

for j = 1:length(fileName)
    imagebin = imread(fileName{j});
    
    figure(1);
    bw_image = rgb2gray(imagebin);
    subplot(3,3,j), imshow(bw_image','DisplayRange',[min(min(bw_image)) max(max(bw_image))]);
    title(fileName(j));
    set(gca,'XTick',[], 'YTick', []);
    
    imagebin = imagebin(:,:,1);
    dists = sum(imagebin>180,1);
    num_bins = 100;%ceil(sqrt(length(dists)));
    M = floor(length(dists)/num_bins);
    temp = [];
    for i = 1:num_bins
        temp(i) = mean(dists( (1+(i-1)*M):1+i*M ));
    end
    %figure(j+9);
    figure(2);
    subplot(3,3,j), plot(dists,'LineWidth',1.5,'Color','k');
    title(fileName(j));
    hold on;
    height = mean(dists);
    plot([1 length(dists)],[height height],'LineWidth',3);
    axis tight;
    set(gca,'XTick',[], 'YTick', []);
    
    %figure(j+18);
    figure(3);
    subplot(3,3,j), plot(temp,'LineWidth',1.5,'Color','k');
    title(fileName(j));
    axis tight;
    set(gca,'XTick',[], 'YTick', []);
    %%
    %temp = find(dists>0);
    %xmax = ceil(length(dists)/bin_size)*bin_size;
    %[h p] = kstest(temp,'CDF',[(bin_size:xmax)' (cdf('Uniform',bin_size:xmax,bin_size,xmax))'])
    
    dists_binned = temp;
    xmax = max(cumsum(dists_binned));
    temp = cumsum(dists_binned);
    temp = temp/xmax;
    
    figure(4);
    subplot(3,3,j), plot(temp, 'LineWidth',2);
    delta = 1/(length(temp)-1);
    temp2 = 0:delta:1;%
    %unifcdf(0:(length(dists)-1),1,length(dists)+2);
    figure(4);
    hold on;
    subplot(3,3,j), plot(temp2, 'LineWidth',2);
    if j ==1
        legend('Empirical CDF','Uniform CDF','Location','NorthWest');
    end
    [h, p] = kstest(temp,'CDF',[(temp2)' (temp2)']);
    if h == 0
        str = sprintf('Accept H0 (uniform), p-value = %1.3f',p);
        title(str);
    else
        str = sprintf('Reject H0 (not uniform), p-value = %1.3f',p);
        title(str);
    end
    
    %     if j == 1 || j == 5
    %         figure(30+j);
    %         plot(temp, 'LineWidth',2);
    %         xticks([0 50 100]);
    %         yticks([0 0.5 1]);
    %         delta = 1/(length(temp)-1);
    %         temp2 = 0:delta:1;%
    %         hold on;
    %         plot(temp2, 'LineWidth',2);
    %         xticks([0 50 100]);
    %         yticks([0 0.5 1]);
    %         if j ==1
    %             legend('Empirical CDF','Uniform CDF','Location','NorthWest');
    %         end
    %         [h, p] = kstest(temp,'CDF',[(temp2)' (temp2)']);
    %         if h == 0
    %             str = sprintf('Accept H0 (uniform), p-value = %f',p);
    %             title(str);
    %         else
    %             str = sprintf('Reject H0 (not uniform), p-value = %f',p);
    %             title(str);
    %         end
    %     end
%     if j < 4
%         figure(42);
%         plot(temp, 'LineWidth',2);
%         xticks([0 50 100]);
%         yticks([0 0.5 1]);
%         delta = 1/(length(temp)-1);
%         temp2 = 0:delta:1;%
%         hold on;
%         plot(temp2, 'LineWidth',2);
%         xticks([0 50 100]);
%         yticks([0 0.5 1]);
%         if j ==1
%             legend('Empirical CDF','Uniform CDF','Location','NorthWest');
%         end
%         [h, p] = kstest(temp,'CDF',[(temp2)' (temp2)']);
%         if h == 0
%             str = sprintf('Accept H0 (uniform), p-value = %f',p);
%             title(str);
%         else
%             str = sprintf('Reject H0 (not uniform), p-value = %f',p);
%             title(str);
%         end
%     end
end
toc;


