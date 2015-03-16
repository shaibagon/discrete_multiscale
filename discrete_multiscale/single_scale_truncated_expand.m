function l = single_scale_truncated_expand(u, v, w, il, varargin)
%
% using a-expand (with truncation) as a "black-box" single scale
% optimization
%
gch=GraphCut('open', u', v, w);
gch=GraphCut('set',gch, il-1);
gch = GraphCut('truncate',gch, true); % allow truncation of non-submodular terms
[gch l]=GraphCut('expand',gch);
l=double(l)+1;
gch=GraphCut('close',gch);
