function l = discreteMultiscaleOptimization_labels(P, dc, sc, W, SingleScaleOptFn)
%
% Performing discrete multiscale optimization of a given energy pyramid
% with coarser LABELS
%
% Usage:
%   l = discreteMultiscaleOptimization_labels(P, dc, sc, W, SingleScaleOptimizationFunction);
%
% Inputs:
%   P   - cell array with a cascade of interpolation matrices 
%         (obtained using buildEnergyPyramid_labels.m)
%   dc  - cell array with a cascade of data costs for the different scales
%         (obtained using buildEnergyPyramid_labels.m)
%   sc  - cell array with a cascade of smoothness terms for the different scales
%         (obtained using buildEnergyPyramid_labels.m)  
%   W   - weighted graph 
%         
%   SingleScaleOptimizationFunction - a function handle for the single
%         scale primal optimization method
%         The function should have the following signiture:
%         function l = SingleScaleOpt(d, v, w, il, opt_params)
%         where
%           d - data cost at the relevant scale
%           v - smoothness cost at relevant scale (always sC)
%           w - wighted graph for relevant scale
%           il - initial guess
%           opt_params - strust with additional parameters
%         output of SingleScaleOpt is a valid discrete labeling for the
%         relevant scale.
%         This package contains several examples of single-scale
%         optimization algorithms:
%         @ICM
%         @single_scale_swap                
%         @single_scale_expand  
%         @single_scale_truncated_expand
%         (the use of large move making single-scale methods requires the
%         installation of GraphCut wrapper:
%         www.wisdom.weizmann.ac.il/~bagon/matlab.html#graphcut)
%
% Output:
%   l   - cell array with the discrete labels of each scale in the pyramid.
%         finest scale solution in l{1}
%


NS = numel(P);

iu = dc{NS+1};

% work down the scales
for si=(NS+1):-1:1
    [n nl] = size(iu);
    
    fprintf(1, 'Working coarse->fine %d (%d) with %dx%d \n',...
        si, NS+1, n, nl);

    % try greedy init
    [iwt, ig] = min(iu,[],2);
    [~, ord] = sort(iwt); % , 'descend');
    
    % note that nodes are ordered fine then coarse
    [l{si}(ord)] = SingleScaleOptFn(dc{si}(ord,:), sc{si}, W(ord,ord),...
        ig(ord));
    
    
    
    if si>1       
        iu = -(l2u(l{si}, nl))*P{si-1}';
    end
    
    
end


