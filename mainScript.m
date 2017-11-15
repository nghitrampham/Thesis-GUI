% main script used to run spectral clustering algorithm 

disp('Starting spectral clustering algorithm');

currentFile = mfilename();
getPath = mfilename('fullpath');

getPath = getPath(1:end-size(currentFile,2));

% add all folders and subfolders to current directory 
path(genpath([getPath 'ZPclustering']), path)
path(genpath([getPath 'data']), path)
% path(genpath([getPath 'func']), path)

% load data 

% load mat4SVM.mat % data is confidential 

 SpectralGUI



