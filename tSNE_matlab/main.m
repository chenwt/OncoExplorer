% show.m
clc; clear; close all;
figure('Color', [1 1 1]);
Iter = 0;
for can = {'pancan', 'brca', 'gbm', 'ov'}
%can = 'ov';
Iter = Iter + 1;
can = can{1};
load(sprintf(['bond_',can,'.mat']));

input = tdfread(sprintf(['index2name_',can,'.txt']));
gene = input.gene;
notch1 = input.notch1;
p53 = input.p53;
pi3k = input.pi3k;
inpathway = [notch1,p53,pi3k];
N = size(gene,1);


train_X = corr;
train_labels = [];
% load('mnist_train.mat');
% ind = randperm(size(train_X,1));
% num = 500;
% train_X = train_X(ind(1:num),:);
% train_labels = train_labels(ind(1:num));
%


% Set parameters
no_dims = 2;
initial_dims = 30;
perplexity = 10;
mappedX=tsne(train_X, [], no_dims, initial_dims, perplexity);

%
ydata = mappedX;
%gscatter(mappedX(:,1), mappedX(:,2));
%

iter = 0;
for pathway = {'notch1', 'p53', 'pi3k'}
    pathway = pathway{1};
iter = iter + 1;
    subplot(3,4,iter*4-3+Iter-1);

%title(pathway{1});
%gscatter(mappedX(:,1), mappedX(:,2));
hold on;
for i = 1:N
    g = gene(i,:);
    %text(ydata(i,1), ydata(i,2), g, 'Color',[0 0 0]);
    plot(ydata(i,1), ydata(i,2),'b.','MarkerFaceColor','b');
end

for i = 1:N
    g = gene(i,:);
    %text(ydata(i,1), ydata(i,2), g, 'Color',[0 0 0]);
    if inpathway(i,iter) == 1
       plot(ydata(i,1), ydata(i,2),'ro','MarkerFaceColor','r'); 
       %text(ydata(i,1), ydata(i,2), g, 'Color',[0 0 0]);
    end
    %if ~isempty(strfind(gene, 'tp53')) || ~isempty(strfind(gene, 'pten'))
    %if strfind(gene, 'tp53')
    %    gene
    %    text(ydata(id,1), ydata(id,2), gene, 'Color',[0 0 1]);
    %end
end
    %tt = ;
    %sprintf(['index2name_',can,'.txt'])
    title(sprintf([can, '-', pathway]));
end
end
%EOF.
