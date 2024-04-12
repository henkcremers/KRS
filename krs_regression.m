function KRS = krs_regression(X,Y,varargin)
% cross-validated regression on distance kernel.
% ==============================================
% USE:
% IN:
%
%
% OUT:
% ==============================================

% select part of data 
% ns = 300;
% X = X(1:ns,:);
% Y = Y(1:ns,:);


%% Defaults
% ===============================================
niter = 10;
kfold = 10;
nperm = 10;

% regression function.. just ridge for now, update later
% ------------------------------------------------------
regfunc = @(X,Y)ridge(Y,X,2,0);

% beta = ridge(distY_total,distX,k,0);

%%  set up the data..
% UPDATE THIS FUNCTION
% [dX,dY,dX_tot,dY_tot] = krs_kernel(X,Y);
X = zscore(X);
Y = zscore(Y);

%% Run the regression
% ======================================================
count = 0;
progressbar_new(['running the regression'])
for i = 1:niter;
    
    CVP = cvpartition(size(Y,1),'kfold',kfold);
    for k = 1:kfold;
        
        clc; disp(['k: ' num2str(k)])
        
        kloc = CVP.test(k);
        
        % split the data
        Xtr  = X(~kloc,:);
        Xte  = X(kloc,:);
        Ytr  = Y(~kloc,:);
        Yte  = Y(kloc,:);
        
        %%  set up the data..
        % When Y is multivariate - change to distance 
        
        if size(Y,2)>1
            [dXtr,dYtr,dXtr_tot,dYtr_tot] = krs_kernel(Xtr,Ytr);
            [dXte,dYte,dXte_tot,dYte_tot] = krs_kernel(Xte,Yte);
            disp('switching to distance regression')
        else   
            dYtr_tot = Ytr;
            dYte_tot = Yte;
            dXtr = Xtr;
            dXte = Xte;
            dXtr_tot = Xtr*ones(size(Xtr,2),1);
            dXte_tot = Xte*ones(size(Xte,2),1);
            disp('standard regression')
        end
        

        
        % partition
        
        %     CVP = cvpartition(length(dY_tot),'kfold',kfold);
        %     for k = 1:kfold;
        %
        count = count +1;
        progressbar_new(count/(niter*kfold))
        %
        %         kloc = CVP.test(k);
        %
        %         % split the data
        %         dXtr  = dX(~kloc,:);
        %         dXte  = dX(kloc,:);
        %         dXtr_tot  = dX_tot(~kloc);
        %         dXte_tot  = dX_tot(kloc);
        %         dYtr_tot  = dY_tot(~kloc);
        %         dYte_tot  = dY_tot(kloc);
        
        rawtr(i,k) = corr(dYtr_tot,dXtr_tot,'type','Spearman');
        rawte(i,k) = corr(dYte_tot,dXte_tot,'type','Spearman');
        
        
        % =================================================================
        % Train the model
        % =================================================================
        
        
        % (1) no correction
        % beta = ridge(dYtr_tot,[dXtr],0.5);
        
        % beta = regfunc(dXtr,dYtr_tot);
        % beta = ridge(dYtr_tot,[dXtr dXtr_cm],0.5);
        
        
        % (2) correct for dependency - Lukas & Steven.
        % ----------------------
%         dX_conf   = depmat(size(Xtr,1));
% % % %         tic
% % % %         dXtr_corr = Xcorrect(dXtr,dX_conf);
% % % %         toc
%         dXtr_corr = [dXtr dX_conf];
        dXtr_corr = dXtr;
        
        
        %beta = ridge(dYtr_tot,dXtr_corr,0.3,0);
        % beta = ridge(dYtr_tot,dXtr_corr,0.3,0);
        
        [beta,FitInfo] = lasso(dXtr_corr,dYtr_tot,'lambda',0.0);
%         model.type = 'lasso';
%         model.B = B;
%         model.FitInfo = FitInfo;
%         model.Yfit = distX*B;
%         model.dfunc = dfunc;
        
        
        % todo: save model parameters
        
        
        % prediction
        % ------------
        % dY_pred = beta(1)+dXte*beta(2:end)';
        
        
        % dXte_con = dXte; % not correction
        % dXte_corr = depcorrect(dXte,size(Xte,1));
        dXte_corr = dXte;
        % dY_pred = dXte_corr*beta(2:size(dXte_corr,2)+1);
        dY_pred = dXte_corr*beta;
        
        % sanity check..
        % dY_fit_tr = [dXtr dXtr_cm]*beta; %this is the version with all
        % variables in the model
        % dY_fit_tr = dXtr_corr*beta(2:end);
        dY_fit_tr = dXtr_corr*beta;
        fit(i,k) = corr(dYtr_tot,dY_fit_tr,'type','Spearman');
        
        % test evaluation features
        % --------------------------
        c(i,k) = corr(dYte_tot,dY_pred,'type','Spearman');
        
        
        %         % permutation
        %         % --------------------------
        %         for n = 1:nperm
        %             r = randperm(length(beta));
        %             beta_perm = beta(r);
        %             dY_perm = dXte*beta_perm;
        %             perm(i,k,n) = corr(dYte_tot,dY_perm,'type','Spearman');
        %
        %         end
        
        % mse =
        
    end
    
    
end

KRS.stats.c  = c;
KRS.stats.fit  = fit;
%KRS.stats.perm = perm;

KRS.stats.rawtr = rawtr;
KRS.stats.rawte = rawte;


% create dep matrix
    function dX_conf = depmat(n)
        % create a binary matrix n*(n-1)/2 by n.
        for j = 1:n
            emp = zeros(n,n);
            conn = 1:n;
            conn(j) = [];
            emp(conn,j) = 1;
            empf = nwa_reshape(emp,'flipdiag');
            emp = emp+empf;
            dX_conf(:,j) = nwa_reshape(emp,'mat2vec');
        end
    end

% correct for dependency


%% plot the results
figure;
subplot(2,1,1)
hist(hist(KRS.stats.fit(:)))
B = mean(mean(KRS.stats.fit));
title(['Fit: ' num2str(B)])
subplot(2,1,2)
hist(hist(KRS.stats.c(:)))
CV = mean(mean(KRS.stats.c));
title(['Predict: ' num2str(CV)])
end

% switch method
%
%     case 'ridge'
%         % ridge regression
%         k = 0.01; %:0.05:0.9;
%         model.type = 'ridge';
%         B = ridge(distY_total,distX,k,0);
%         for j = 1:size(B,2);
%             distY_pred(:,j) = B(1,j)+distX*B(2:end,j);
%         end
%         err = fsdist_error(distY_pred,distY_total);
%         model.dfunc = dfunc;
%         model.B = B;
%         model.trErr = err;
%         model.Yfit = distY_pred;
%         model.Ytrain = distY_total;
%
%     case 'rvr'
%         [w,alpha,beta,ll] = prt_rvr(distX*distX',distY_total);
%         predictions=distX*distX'*w(1:end-1)+w(end);
%         model.type = 'rvr';
%         % model.B = B;
%         model.Yfit = predictions;
%         model.dfunc = dfunc;
%
%     case 'lasso'
%         [B,FitInfo] = lasso(distX,distY_total,'lambda',1);
%         model.type = 'lasso';
%         model.B = B;
%         model.FitInfo = FitInfo;
%         model.Yfit = distX*B;
%         model.dfunc = dfunc;
%
%    case 'svr'
%         Mdl = fitrsvm(distX,distY_total);
%         model.type = 'svr';
%         model.B = Mdl.Beta;
%         %model.FitInfo = FitInfo;
%         model.Yfit = predict(Mdl,X);
%
%         model.dfunc = dfunc;
%
%     case 'tree'
%         model.type = 'tree';
%         model.dfunc = dfunc;
%         tree = fitrtree(distX,distY_total,'MaxNumSplits',100);
%         %tree = fitrtree(distX,distY_total);
%         model.tree = tree;
%         Yfit = predict(tree,distX);
%         model.Yfit = Yfit;
%         model.Ytrain = distY_total;
%         err = fsdist_error(Yfit,distY_total);
%         model.trErr = err;
%
% end
% end

