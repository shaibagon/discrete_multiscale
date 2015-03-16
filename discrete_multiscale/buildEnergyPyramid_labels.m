function [P dc sc] = buildEnergyPyramid_labels(Dc, sC, NS)
%
% Construct multi-scale LABEL coarsening 
%
% Usage:
%    [P dc sc] = buildEnergyPyramid_labels(Dc, sC, NS)
%
% Inputs:
%   Dc - data term 
%   sC - smoothness
%   NS - max number of scales
%
% Outputs:
%   P  - interp. matrix
%   dc - multi scale unary term
%   sc - multi-scale smoothness
%        (note that W is unchanged)

%
% typical values:
beta = .5;
thr = 0;


n = size(sC,1);

P = cell(1,NS);

dc = cell(1,NS+1);
sc = cell(1,NS+1);

maxD = 2;


sc{1} = sC;
dc{1} = Dc;
assert(all(nonzeros(sc{1})>0));

fine = 1:n;


for si=1:NS
    
    [tmp_c P{si}] = fine2coarse(V2Aff(sc{si}), beta, maxD);    
    
    if size(P{si},2) <= 2 || ...
            size(P{si},2) == size(P{si},1)  
        % stop growing
        NS = si-1;
        P = P(1:NS);        
        
        sc = sc(1:(NS+1));
        return;
    end
    
    
    fine = fine(tmp_c);
    
    sc{si+1} = P{si}'*sc{si}*P{si};       
    dc{si+1} = dc{si}*P{si};
    
end





%-------------------------------------------------------------------------%
function [c P] = fine2coarse(V, beta, maxD)
%
% Performs one step of coarsening
%

n = size(V,1);
V(1:(n+1):end) = 0;

% normalized 
nV = full(sum(V,2)); nV(nV==0) = 1;
nV = bsxfun(@rdivide, V, nV);

% % % select coarse nodes
c = ChooseCoarseGreedy_mex(sparse(nV), [2:2:n, 1:2:n], beta);


% compute the interp matrix
ci=find(c);
P = sparse(V(:,ci));


% enforce max degree
f = find(~c);
deg = sum(P(~c,:)>0,2);
if any(deg>maxD)
    P(f(deg>maxD),:) = AdjustPdeg( P(f(deg>maxD),:)', maxD)' ;
end



% make sure coarse points are directly connected to their fine counterparts
[jj ii pji] = find(P'); 
sel = ~c(ii); % select only fine points
mycat = @(x,y) vertcat(x(:),y(:));
P = sparse( mycat(jj(sel), 1:sum(c)),...
    mycat(ii(sel), ci),...
    mycat(pji(sel), ones(1,sum(c))), size(P,2), size(P,1))';

P = bsxfun(@rdivide, P, sum(P,1));

assert(all(P(:)>=0));
assert(all( abs(sum(P,1)-1) < 1e-10 ) );

%-------------------------------------------------------------------------%
function A = V2Aff(V)
%
% converts distances matrix into affinity matrix
%
V = V./max(abs(V(:)));

A = 1./(V+eps);

