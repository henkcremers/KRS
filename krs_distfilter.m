function [fdat,filter,nf] = krs_distfilter(X,varargin)

dfunc = 'euclidean';
filter = [];
doplot = 0;
fthresh = 3;
dat{1} = X;

%%  get input
%-------------------
% if length(varargin)>0
% if isnumeric(varargin{1})
%     dat{2} = varargin{1};
% end
datcount = 1;
for i = 1:length(varargin)
    arg = varargin{i};
    if isnumeric(arg)
        datcount = datcount +1;
        dat{datcount}=arg;
        
    elseif ischar(arg)
        switch arg
            case 'filter',dofilter = 1;
            case 'dfunc', dfunc = varargin{i+1};
            case 'fthresh', fthresh = varargin{i+1};
            case 'plot', doplot = 1;
        end
    end
end


%% filter each data modality
nf = [];
for j = 1:length(dat)
    
    d = dat{j};
    if ismatrix(d);
        dist = pdist(zscore(d), dfunc);
    elseif isvector(d);
        dist = d;
    end
    
    %% filter outliers..
    distM = squareform(dist);
    mD = mean(distM);
    %smD = mD./std(mD);
    smD = zscore(mD);
    f = (smD<=fthresh);
    filter(:,j) = f';
    nf(j) = sum(f==0);
    
%     if doplot == 1;
%         fsdist_mds(d,'cdat',f)
%     end
    
    % dist = pdist(zscore(dat), dfunc);
    % distM = squareform(dist);
    % mD = mean(distM);
    % disp(['filtered out: ' num2str(sum(filter==0)) ' subjects'])
end
filter = sum(filter,2)==length(dat);

for j = 1:length(dat)
    fdat{j} = dat{j}(filter,:);
end

if length(dat)==1;
    fdat = fdat{1};
end

% plot data
if doplot == 1
    
    figure
    for j = 1:length(dat)
        
        hold on
        p1 = j;
        subplot(1,datcount,p1)
        krs_plot_mds(dat{j});
        hold on
        
%         p2 = datcount+j;
%         subplot(2,datcount,p2)
%         krs_plot_mds(fdat{j});
%         hold on
        
    end
    
    figure
    for j = 1:length(dat)

        subplot(1,datcount,j)
        krs_plot_mds(fdat{j});
        hold on
        
    end
    
end

return



