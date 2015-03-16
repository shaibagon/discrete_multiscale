function [P dc w] = buildEnergyPyramid(Dc, sC, W, NS)
%
% Construct multi-scale variable coarsening in an adaptive manner
%
% Usage:
%    [P dc w] = buildEnergyPyramid(Dc, sC, W, NS)
%
% Inputs:
%   Dc - data term 
%   sC - smoothness
%   w  - problem graph
%   NS - max number of scales
%
% Outputs:
%   P  - interp. matrix
%   dc - multi scale unary term
%   w  - multi-scale w
%        (note that sC is unchanged)
%


% some parameters - these settigs were used throughout our experiments
beta = 0.2; % controls the ratio of coarse to fine variables
maxD = 3; % max number of coarse representatives for each fine node
K = 10; % number of random inits
rLim = 10; % number of iterations for random inits

[n nl] = size(Dc);

dc = cell(NS+1,1);
P = cell(NS,1);
w = cell(NS+1,1);

dc{1} = Dc;
w{1} = W;

warning('off','ICM:itr');

weps = 1e-6;
wsig = -reallog(weps)/max(sC(:));

% work up the scales
for si=1:NS
    
    fprintf(1, 'Working fine->coarse %d (%d) with %dx%d \n',...
        si, NS, n, nl);

    
    %
    % 1. estimate correlations under given energy between variables
    l = zeros(n,K);    
    [ii jj] = find(w{si});
    dij = zeros(numel(ii),K);
    
    
    % try greedy init
    [iwt, ig] = min(dc{si},[],2);
    [~, ord] = sort(iwt, 'descend');  
        
    for sti=1:K % K random inits
                
        l(ord,sti) = ICM(dc{si}(ord,:), sC,...
            w{si}(ord,ord),...
            ig, rLim);     
        
        
        dij(:,sti) = sC( ([l(ii,sti) l(jj,sti)]-1)*[1; nl] + 1 );
    
        % random init for next iteration
        ord = randperm(n);
        ig = randi(nl, n, 1);
    end

    % variable correlations:
    cij = exp( -wsig*mean(dij,2) );
    cij(cij<=weps)=0;
    
    %
    % 2. from correlations to interpolation P (AMG algorithm)
    [~, P{si}] = fine2coarse(sparse(ii, jj, cij, n, n), beta,  -std(l,1,2), maxD);
    
    if size(P{si},2) <= 2 || ...
            size(P{si},2) + 3 >= size(P{si},1) % less than 3 var. reduction
        % stop growing
        NS = si-1;
        P = P(1:NS);
        dc = dc(1:(NS+1));
        w = w(1:(NS+1));
        break;
    end
  
    %
    % 3. Given P construct coarse energy
    w{si+1} = P{si}'*w{si}*P{si};
    diagw = spdiags(w{si+1},0);
    
    
    % remove diagonal
    w{si+1} = w{si+1} - spdiags( diagw, 0, size(w{si+1},1), size(w{si+1},2));

    dc{si+1} = P{si}'*dc{si}; % assuming daigv==0 + full(diagw(:))*(diagv(:)');        
    [n nl] = size(dc{si+1});
end

warning('on','ICM:itr');

%-------------------------------------------------------------------------%
function [c P] = fine2coarse(W, beta, weights, maxD)
%
% Performs one step of coarsening
%

% normalized attraction part
nW = full(sum(W,2)); 
nW(nW==0) = 1;
nW = bsxfun(@rdivide, W, nW); %  spmtimesd(W, 1./nW, []);


[~, ord] = sort(weights, 'descend');
c = ChooseCoarseGreedy_mex(nW, ord, beta);

% compute the interp matrix
ci=find(c);
P = W(:,ci);

% enforce max degree
f = find(~c);
deg = sum(P(~c,:)>0,2);
if any(deg>maxD)
    P(f(deg>maxD),:) = AdjustPdeg( P(f(deg>maxD),:)', maxD)' ;
end

P = bsxfun(@rdivede, P, sum(P,2)); % spmtimesd(P, 1./full(sum(P,2)), []);
% make sure coarse points are directly connected to their fine counterparts
[jj ii pji] = find(P'); 
sel = ~c(ii); % select only fine points
mycat = @(x,y) vertcat(x(:),y(:));
P = sparse( mycat(jj(sel), 1:sum(c)),...
    mycat(ii(sel), ci),...
    mycat(pji(sel), ones(1,sum(c))), size(P,2), size(P,1))';

% numerical stability - remove close to zero elements
P = spfun(@(x) x.*(x>1e-2), P);
P = bsxfun(@rdivide, P, sum(P,2)); % spmtimesd(P, 1./full(sum(P,2)), []);

%-------------------------------------------------------------------------%
