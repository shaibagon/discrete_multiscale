function u=l2u(l,nl)
%
% convert with fixed number of labels
%
n = numel(l);
u=zeros(n,nl);

u([1 n]*([1:n;l(:)']-1)+1) = 1;
