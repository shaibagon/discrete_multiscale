function l = ICM_Potts(u, w, il, LIMIT)
%
% ICM for pair-wise potts model
%
% Minimize energies over n variables with l labels defined by:
%   E(l) = \sum_i u(i, li)  + \sum_ij w_ij [li != lj]
% where [.] is the indicator function
%
% usage:
%   l = ICM_Potss(u, w, il)
%
% Inputs:
%   u  - unary term (n)x(l)
%   w  - sparse adjacency matrix (n)x(n)
%   il - initial guess (n)x1 \in {1,...,l}^n
%   LIMIT - max number of iterations
%



n = size(w,1);

if nargin < 4
    LIMIT = 1000*ceil(reallog(n));
end

% ignore diagonal
w = w - spdiags( spdiags(w,0), 0, n, n);

itr = 1;
ord = 1:n;
while itr<LIMIT
    l(ord) = potts_icm_iter_mex(u(:,ord), w(ord,ord), il(ord));
    
    if isequal(l,il) && all(l>0)
        break;
    end
    il = l;
    itr = itr+1;
    
    iord = ord(end:-1:1);
    % revese order
    l(iord) = potts_icm_iter_mex(...
        u(:,iord), w(iord,iord), il(iord));
    
    if isequal(l,il) && all(l>0)
        break;
    end
    il = l;
    itr = itr+1;
    ord = randperm(n);
end


if itr == LIMIT
    warning('ICM_Potts:itr','reached LIMIT');
end
