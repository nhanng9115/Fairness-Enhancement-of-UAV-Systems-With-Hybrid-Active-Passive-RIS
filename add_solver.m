soft_store = 'your folder containing YALMIP and MOSEK'; % <--- edit here

pwd
currentFolder = pwd; % save current folder


% add YALMIP
cd(strcat(soft_store,'YALMIP-master')); addpath(genpath(pwd))

% add mosek
cd(strcat(soft_store,'mosek\9.2\toolbox')); addpath(genpath(pwd))
mosek_env = [getenv('PATH') strcat(';',soft_store,'mosek\9.2\tools\platform\win64x86\bin')];
setenv('PATH', mosek_env);

% back to current folder
cd(currentFolder)
