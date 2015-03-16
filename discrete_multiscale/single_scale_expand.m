function l = single_scale_expand(u, v, w, il, varargin)
%
% Using a-expand as a "black-box" single scale optimization
%
gch=GraphCut('open', u', v, w);
gch=GraphCut('set',gch, il-1);
[gch l]=GraphCut('expand',gch);
l=double(l)+1;
gch=GraphCut('close',gch);
