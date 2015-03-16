function l = single_scale_swap(u, v, w, il, varargin)
%
% Using ab-swap as a "black-box" single-scale optimization 
%
gch=GraphCut('open', u', v, w);
gch=GraphCut('set',gch, il-1);
[gch l]=GraphCut('swap',gch);
l=double(l)+1;
gch=GraphCut('close',gch);
