function [de se te] = Energy(dc, v, w, l)
%
% Computes energy
%
nl=size(dc,1);
u=l2u(l,nl);

de = dc.*(u');
de = sum(de(:));

se = v.*(u'*(tril(w,0)*u));
se = sum(se(:)); % w is symmetric - we want to avoid double counting

te=se+de;

if nargout<=1
    de = te;
end
