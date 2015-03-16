function [l itr] = ICM(u, v, w, il, LIMIT)
%
% ICM for pair-wise energy
%
% Minimize energies over n variables with l labels defined by:
%   E(l) = \sum_i u(i, li)  + \sum_ij w_ij v(li,lj)
%
% usage:
%   l = ICM(u, v, w, il)
%
% Inputs:
%   u  - unary term (n)x(l)
%   v  - label-to-label cost (l)x(l)
%   w  - sparse adjacency matrix (n)x(n)
%   il - initial guess (n)x1 \in {1,...,l}^n
%   LIMIT - max number of iterations
%
% Note: if V is Potts (i.e. V(li,li)=0, for all li, and all other entries
% for li!=lj V(li,lj) are equal), then ICM_potts should be faster.
%

n = size(w,1);
L = size(v,1);
if size(u,1)~= L
    u = u';
end
assert(numel(il)==n && all(size(u)==[L n]));

if nargin < 5 || ~isnumeric(LIMIT)
    LIMIT = 1000*ceil(reallog(n));
end

% ignore diagonal
w = w - spdiags( spdiags(w,0), 0, n, n);

itr = 1;

ord = 1:n;

while itr<LIMIT
    l(ord) = gen_eng_icm_iter_mex(u(:,ord), v, w(ord,ord), il(ord));
    
    if isequal(l,il) && all(l>0)
        break;
    end
    il = l;
    itr = itr+1;
    
    iord = ord(end:-1:1);
    % revese order
    l(iord) = gen_eng_icm_iter_mex(...
        u(:,iord), v, w(iord,iord), il(iord));
    
    if isequal(l,il) && all(l>0)
        break;
    end
    il = l;
    itr = itr+1;
    ord = randperm(n);
end



if itr >= LIMIT
    warning('ICM:itr','reached LIMIT');
end
