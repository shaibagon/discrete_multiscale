function discrete_multiscale_demo()
%
% Demo energies from middlebury
%

addpath('../discrete_multiscale');
addpath('../icm');

demo_energies = {'demo_denoise_penguin', 'demo_stereo_ven'};

%%-----------------------------------------------------------------------%%
% First demo: Coarsening variables
fprintf(1, 'Discrete Multiscale Energy minimization demo:\n');
fprintf(1, 'Coarsening variables\n');


% read energy - 
for di=1:numel(demo_energies)
    fprintf(1, '\n\tEnergy: %s\n', demo_energies{di});
    
    % load energy
    load([demo_energies{di},'.mat']);
    
    % construct energy pyramid
    % (number of scales is chosen very large so the function can build the
    % pyramid as coarse as it can)
    [P dc w] = buildEnergyPyramid(Dc, sC, W, 100);
    
    % optimize
    licm = discreteMultiscaleOptimization(P, ...
        dc, sC, w, ...
        @ICM);
    lswap = discreteMultiscaleOptimization(P, ...
        dc, sC, w, ...
        @single_scale_swap);
    lexpand = discreteMultiscaleOptimization(P, ...
        dc, sC, w, ...
        @single_scale_truncated_expand);
end

%%-----------------------------------------------------------------------%%
% Second demo: Coarsening labels

for di=1:numel(demo_energies)
    fprintf(1, '\n\tEnergy: %s\n', demo_energies{di});
    
    % load energy
    load([demo_energies{di},'.mat']);
    
    % construct energy pyramid
    % (number of scales is chosen very large so the function can build the
    % pyramid as coarse as it can)
    [P dc sc] = buildEnergyPyramid_labels(Dc, sC, 100);
    
    lswap = discreteMultiscaleOptimization_labels(P, dc, sc, W, @single_scale_swap);
end
