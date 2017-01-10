% show.m
clc; clear; close all;

can = 'brca';
input = tdfread(sprintf(['linked_',can,'_n.txt']));
Source = input.Source;
Target = input.Target;
Weight = input.Weight;
Distance = 1./Weight/100;

%hist(Distance,1000)

%
perplexity = 14;

input = tdfread(sprintf(['lut_',can,'.txt']));
Id = input.Id;
Gene = input.Gene;

%min(Distance)
maxDist = max(Distance)
maxDist = maxDist;
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

% Set parameters
no_dims = 2;
initial_dims = 30;

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

fig = figure();
hold on;

gscatter(ydata(:,1), ydata(:,2), [], 'r',[],20);

set(fig,'Position',[-20, 800, 900, 900]);


if true
for i = 1:N
    id = Id(i);
    gene = Gene(i,:);
    text(ydata(id,1), ydata(id,2), gene, 'Color',[0 0 0]);
    if ~isempty(strfind(gene, 'tp53')) || ~isempty(strfind(gene, 'pten'))
    %if strfind(gene, 'tp53')
        gene
        text(ydata(id,1), ydata(id,2), gene, 'Color',[0 0 1]);
    end
end
end


%EOF.
