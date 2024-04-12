function [cdat,dist] = krs_plot_mds(dat,varargin)

% defaults
% ------------------
dofilter = 0;
idx = [];
dfunc = 'euclidean';
cdat = [];
nK = 2;
fig=0;

%%  get input
%-------------------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'dfunc', dfunc = varargin{i+1};
            case 'cdat', cdat = varargin{i+1};
            case 'idx', idx = varargin{i+1};
            case 'fig', fig = 1;
        end
    end
end

% clustering ..
if isempty(idx)
    % [idx,C,SUMD,D] = kmeans(dat,nK,'Distance',dfunc,'Replicates',5);
    % 'eucledian' not working? see help kmeans...
    if size(dat,1)==1;
        dat = dat';
    end
    [idx,C,SUMD,D] = kmeans(dat,nK,'Distance','sqeuclidean','Replicates',5);
end


%% calc distance or take predictions from a model
% 20200705 - Fix to take squaare distance matrix...
s = size(dat);

if s(1)==s(2)
    dist = dat;
    mD = mean(dat);
    
elseif isvector(dat);
    if size(dat,2)==1;
        dist = dat;
    else
        dist = dat';
    end
    try
        distM = squareform(dist);
    catch
        dist = pdist(zscore(dat), dfunc);
        distM = squareform(dist);
    end
    mD = mean(distM);
    
else
    %if s(1)>1 & s(2)>1 - this is not relevant
    dist = pdist(zscore(dat), dfunc);
    distM = squareform(dist);
    mD = mean(distM);
end


%% filter for outliers..
% if dofilter == 1;
% filterval = 6;
% zmD = zscore(mD);
% filter = abs(zmD)<filterval;
% dat = dat(filter,:);
% idx = idx(filter);
% dist = pdist(zscore(dat), dfunc);
% distM = squareform(dist);
% mD = mean(distM);
% cdat = cdat(filter);
% disp(['filtered out: ' num2str(sum(filter==0)) ' subjects'])
% end

%% color data
if isempty(cdat)
    ns = 7;
    ls = linspace(min(mD),max(mD),ns);
    % recode
    cdat = [];
    for j = 1:(ns-1)
        cdat(mD>=ls(j)&mD<=ls(j+1))=j;
    end
    cdat(idx==2) = cdat(idx==2)*-1;
end


%% MDS
[mdsY,e] = cmdscale(dist);
%[mdsY,stress] = mdscale(distY,2);

n = size(dat,1);
%sv = size(varargin,2) +1;
nK = length(unique(idx));

%% create a figure
if fig==1
    figure;
    hold on
end
% colormap hot
nwa_colors %- gives an error, check why?!
colormap(cmap_redblue);
scatter(mdsY(:,1),mdsY(:,2),50,cdat,'filled')

end



