function installDiscreteMultiscaleOpt()
%
% compile all associated cpp code
% update path
%

flags = '-O -largeArrayDims';

subd = {'discrete_multiscale', 'icm'};

for si=1:numel(subd)
    
    fprintf('Entering %s...\n', subd{si});
    
    cd(subd{si});
    
    cfls = dir('*.cpp');
    
    for ci=1:numel(cfls)
        
        cmd = sprintf('mex %s %s', flags, cfls(ci).name);
        
        fprintf(1, '\tCompiling %s (%s)\n', cfls(ci).name, cmd);
        
        eval(cmd);
    end
    fprintf(1,'Exiting %s.\n', subd{si});
    cd('..');
end
fprintf(1,'Done mex compilation\n');

% setting the path
fprintf(1, '\nUpdating path\n');
cwd = pwd();
addpath( fullfile(cwd, 'icm') ); % add apth to ICM
addpath( fullfile(cwd, 'discrete_multiscale') );

fprintf(1, '\nInstallation complete.\n');
fprintf(1, 'To get started run:\n\tdemo/discrete_multiscale_demo();\n\n');

