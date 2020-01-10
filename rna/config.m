function config()
% configures functions and libraries that driver.m uses

% confirm that necessary functions and libraries exist
path = pwd;
if ~exist(strcat(path,'/Libraries'), 'dir')
    error('Couldn''t locate Libraries folder; please place in working directory');
elseif ~exist(strcat(path,'/Functions'), 'dir')
    error('Couldn''t locate Functions folder; please place in working directory');
end

% add folders to matlab search path
addpath(genpath(strcat(path, '/Libraries')));
addpath(genpath(strcat(path, '/Functions')));
    