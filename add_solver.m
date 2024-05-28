% if strcmp(getenv('COMPUTERNAME'),'OY1902032') == 1 || strcmp(getenv('COMPUTERNAME'),'JALAVA1') == 1
% soft_store = 'C:\Users\thanguye20\OneDrive - Oulun yliopisto\Research\OptSoftware\';
% soft_store = 'F:\OneDrive - Oulun yliopisto\Research\OptSoftware\';
soft_store = 'C:\Users\thanguye20\OneDrive - Oulun yliopisto\Research\OptSoftware\';

pwd
currentFolder = pwd; % save current folder
% add CVX
% cd(strcat(soft_store,'cvx-w64\cvx')); cvx_setup

% add YALMIP
cd(strcat(soft_store,'YALMIP-master')); addpath(genpath(pwd))

% add mosek
cd(strcat(soft_store,'mosek\9.2\toolbox')); addpath(genpath(pwd))
mosek_env = [getenv('PATH') strcat(';',soft_store,'mosek\9.2\tools\platform\win64x86\bin')];
setenv('PATH', mosek_env);

% add Manifold Opt
% addpath(pwd);
% cd(strcat(soft_store,'manopt')); 
% addpath(genpath(pwd));
% cd(currentFolder)

% back to current folder
cd(currentFolder)
% end