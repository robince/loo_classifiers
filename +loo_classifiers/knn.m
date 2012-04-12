function [conmtx,info] = knn(data,k,metric)
%LOO_CLASSIFIERS.KNN k-Nearest Neighbours with LOOCV.
%   [CONMTX,INFO] = LOO_CLASSIFIERS.KNN(DATA,K,METRIC) performs leave-one-out
%   cross validation using a k-nearest neighbours classifier.
%   Assigment is by majority rule with nearest point tie break.
%
%   DATA should be a (Nftr, Ncls, Ntrl) 3D array. It is assumed that there 
%   are the same number of trials for each class, and class priors are 
%   uniform.
%
%   K is the integer number of neighbours. METRIC is the distance metric 
%   to use (see CREATENS documentation).
%
%   CONMTX is the (Ncls, Ncls) confusion matrix where the i,j entry gives
%   the percentage of class i trials which were classsified as class j.
%   mean(diag(CONMTX)) gives the average correct performance.
%
%   INFO is the information in the confusion matrix.
%
%nsmethod = 'kdtree';
nsmethod = 'exhaustive';
%includeties = false;
includeties = true;

[Nftr, Ncls, Ntrl] = size(data);
% flatten trials and classes
datflt = data(:,:)';
% build search object
% kdtree is faster, but exhaustive allows more metrics
% including user defined function

NS = createns(datflt,'NSMethod',nsmethod,'Distance',metric);
cls = repmat(0:(Ncls-1),[1 Ntrl]);
prdstm = zeros(Ntrl,Ncls);
conmtx = zeros(Ncls,Ncls);
prctrl = 100 / Ntrl;

% do all nearest neighbours at once
% (only one pdist2 calculation in exhaustive search case)
% k+1 neighbours because we will obtain current point in results
[idx, dst] = knnsearch(NS, datflt,'k',k+1,'IncludeTies',includeties);
if ~includeties
   idx = num2cell(idx,2);
   dst = num2cell(dst,2);
end
curidx = 1;
for ti=1:Ntrl
    for ci=1:Ncls
        % remove current leave out
        % knncls indexes classes from 0 (for bincount)
        if length(idx{curidx}==1)
            % with includeties, sometimes only the point itself is returned
            % not sure why this is
            [idx{curidx} dst{curidx}] = knnsearch(NS, datflt(curidx,:),'k',k+1,'IncludeTies',false);
        end
        knncls = cls(idx{curidx}(idx{curidx}~=curidx));
        % count votes
        votcnt = bincount(knncls,Ncls);
        [cnt, prestm] = max(votcnt);
        if cnt < k/2
            % might have a tie
            tieidx = find(votcnt==cnt);
            Ntie = length(tieidx);
            if Ntie>1
                % have a tie
                tiemin = zeros(1,length(tieidx));
                for tie=1:Ntie
                   tiemin(tie) = min(dst{curidx}(knncls == (tieidx(tie)-1)));
                end
                [mindst, mintie] = min(tiemin);
                % if distances are also tied chose at random
                % use == since usually will be integers
                dsttieidx = find(tiemin==mindst);
                if length(dsttieidx)>1
                    mintie = dsttieidx( randi(length(dsttieidx),1) );
                end
                prestm = tieidx(mintie);
            end % tie
        end % tie check

        prdstm(ti,ci) = prestm;
        conmtx(ci, prestm) = conmtx(ci, prestm) + prctrl;
        curidx = curidx + 1;
    end % classes
end % trials

opts.method = 'dr';
opts.bias   = 'pt';
opts.btsp   = 0;
opts.nt     = Ntrl;
info = information(reshape(prdstm,[1 Ntrl Ncls]),opts,'I');

