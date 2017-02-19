clc; clear; close all;
fileID = fopen('graph_brca.txt');
C = textscan(fileID,'%f %f %f');
i = C{1,1};
j = C{1,2};
v = C{1,3};
network = sparse(i,j,v,529,529);
save('graph_brca.mat', 'network');

%%
clc; clear; close all;

load('graph_brca.mat');