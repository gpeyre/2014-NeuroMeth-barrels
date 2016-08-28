clear all;
close all;

% Go into data directory
data_dir = 'C:\Users\Isabelle\Dropbox\master-lorraine\data';
cd(data_dir);

% Add directories to the path
addpath(genpath('C:\Users\Isabelle\Dropbox\master-lorraine\code\lib'));
addpath(genpath('C:\Users\Isabelle\Dropbox\master-lorraine\code\Registration'));
addpath(genpath('C:\Users\Isabelle\Dropbox\master-lorraine\code\GUI\Barrel_Cortex_registration'));

%% Run registration GUI !
load_gui