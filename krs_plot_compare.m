function [] = krs_plot_compare(Y,X,varargin)


%% calc distance and kmeans clustering. 
dfunc = 'euclidean'; 
% distY = pdist(zscore(Y), dfunc); 
% distX = pdist(zscore(X), dfunc); 
% [idx,C,SUMD,D] = kmeans(Y,nK,'Distance','sqeuclidean','Replicates',5);
% [mdsY,e] = cmdscale(distY);
% maxerrY = max(abs(pdist(Y) - pdist(mdsY(:,1:2))));

n = size(Y,1);
sv = size(varargin,2) +1;

%% create a figure
figure;

% for s = 1:sv;

% plot 1
%loc1 = 1+(s-1)*3;
subplot(1,2,1);
[cdatY,dY] = krs_plot_mds(Y);
% hold on
% for k = 1:nK
% %for k = 1:n
% colors = {'r.' 'b.' 'y.' 'g.' 'k.' 'm.' 'c.'};     
% colornum = k - (ceil((k/length(colors)))-1)*length(colors);
% plot(mdsY(idx==k,1),mdsY(idx==k,2),colors{colornum},'MarkerSize',12)
% %plot(mdsY(k,1),mdsY(k,2),colors{colornum},'MarkerSize',12)
% end
% xlabel('maxerrY')
title 'data Y: MDS Y + K-means Y'
hold off

% plot 2 
    
%loc2 = 2+(s-1)*3;
subplot(1,2,2);
[cdatX,dX] = krs_plot_mds(X,'cdat',cdatY);
title 'data X: MDS X + K-means Y'


% subplot(1,3,3)
% scatter(dX,dY);
% c = corr(dX',dY','Type','Spearman');
% title(['Correlation: ' num2str(c)])
% xlabel('Distance X')
% ylabel('Distance Y')
        

% 
% if s == 1
%    Xdat = X;
% else
%    Xdat = varargin{s-1};
% end
% disterr = fsdist_mindist(Xdat,Y);   
% distX = pdist(zscore(Xdat), dfunc); 
% [mdsX,e] = cmdscale(distX);
% maxerrX = max(abs(distX) - pdist(mdsX(:,1:2)));
% % % nomrmalize X ? 
% mX = max(mdsX); 
% mX = repmat(mX,size(X,1),1);
% mdsX = mdsX./mX;

% hold on
% for k = 1:nK
% %for k = 1:n
% colors = {'r.' 'b.' 'y.' 'g.' 'k.' 'm.' 'c.'};     
% colornum = k - (ceil((k/length(colors)))-1)*length(colors);
% plot(mdsX(idx==k,1),mdsX(idx==k,2),colors{colornum},'MarkerSize',12)
% %plot(mdsX(k,1),mdsX(k,2),colors{colornum},'MarkerSize',12)
% end
% title (['MDS X + K-means Y. Distance Error: ' num2str(disterr)])
% xlabel(maxerrX)
% hold off

% % plot 3 
% %loc3 = 3+(s-1)*3;
% subplot(1,3,3);
% [A,B,r,U,V,STATS] = canoncorr(X,Y);
% plot(U(:,1),V(:,1),'.')
% title(['CanCorr: ' num2str(r(1))])
%xlabel('0.0025*Disp+0.020*HP-0.000025*Wgt')
%ylabel('-0.17*Accel-0.092*MPG')

% disp(['max correlation: ' num2str(max(max(corr(Y,X)))
% disp(['max correlation: ' num2str(max(max(abs(corr(Y,X)))))])
% disp(['max CanCorr: ' num2str(r(1))])
% end

% also compare this... NEEDS TO BE FIXED FOR prediction data.
[eCorr cCorr compdat] = nwa_compare_conn(nwa_reshape(dY,'vec2mat'),nwa_reshape(dX,'vec2mat'),'plot');

end

