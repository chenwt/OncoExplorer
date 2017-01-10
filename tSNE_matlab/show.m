% show.m
clc; clear; close all;

can = 'brca';
input = tdfread(sprintf(['linked_',can,'_n.txt']));
Source = input.Source;
Target = input.Target;
Weight = input.Weight;
Distance = 1./Weight;



input = tdfread(sprintf(['lut_',can,'.txt']));
Id = input.Id;
Gene = input.Gene;

%txt1 = '\leftarrow sin(\pi) = 0';
%text(x1,y1,txt1)

%min(Distance)
maxDist = max(Distance)
maxDist = maxDist*1;
%hist(Distance,500)
N = size(Id,1);
% D is the squared distance matrix.
D = maxDist* ( ones(N,N) - diag(ones(N,1)) );

num_nonzero = size(Distance,1);

for i = 1:num_nonzero
    src = Source(i,1);
    dst = Target(i,1);
    d = Distance(i,1);
    D(src,dst) = d;
    D(dst,src) = d;
end

% Diagnal elements should be zeros.


%
%load('mnist_train.mat');
%ind = randperm(size(train_X,1));
%num = 500;
%train_X = train_X(ind(1:num),:);
%train_labels = train_labels(ind(1:num));

% Set parameters
no_dims = 2;
initial_dims = 30;
perplexity = 30;
%ydata=tsne(train_X, [], no_dims, initial_dims, perplexity);
%X = train_X;
%ydata = tsne(X, labels, no_dims, initial_dims, perplexity)
labels = [];
% First check whether we already have an initial solution
if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
    perplexity = initial_dims;
else
    initial_solution = false;
end

% % Normalize input data
% X = X - min(X(:));
% X = X / max(X(:));
% X = bsxfun(@minus, X, mean(X, 1));
% 
% % Perform preprocessing using PCA
% if ~initial_solution
%     disp('Preprocessing data using PCA...');
%     if size(X, 2) < size(X, 1)
%         C = X' * X;
%     else
%         C = (1 / size(X, 1)) * (X * X');
%     end
%     [M, lambda] = eig(C);
%     [lambda, ind] = sort(diag(lambda), 'descend');
%     M = M(:,ind(1:initial_dims));
%     lambda = lambda(1:initial_dims);
%     if ~(size(X, 2) < size(X, 1))
%         M = bsxfun(@times, X' * M, (1 ./ sqrt(size(X, 1) .* lambda))');
%     end
%     X = bsxfun(@minus, X, mean(X, 1)) * M;
%     clear M lambda ind
% end

% Compute pairwise distance matrix
% %%
% X = [-2;-1;0;1;3]
% X .^ 2
% sum_X = sum(X .^ 2, 2);
% D = bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (X * X')));

%

%size(X): Nxn
%size(D): NxN

% Compute joint probabilities
P = d2p(D, perplexity, 1e-5);                                           % compute affinities using fixed perplexity
clear D

% Run t-SNE
if initial_solution
    ydata = tsne_p(P, labels, ydata);
else
    ydata = tsne_p(P, labels, no_dims);
end







%
%gscatter(mappedX(:,1), mappedX(:,2),train_labels);
%figure;
%
fig = figure();
gscatter(ydata(:,1), ydata(:,2));
%
hold on;
for i = 1:N
    id = Id(i);
    gene = Gene(i,:);
    text(ydata(id,1), ydata(id,2), gene, 'Color',[0 0 0]);
    
    if strcmp(gene, 'pten')
        gene
        text(ydata(id,1), ydata(id,2), gene, 'Color',[0 0 1]);
    end
end

%
set(fig,'Position',[-20, 800, 900, 900]);
%1
% Id = input.Id;
% Gene = input.Gene;

%EOF.
