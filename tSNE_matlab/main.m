% show.m
clc; clear; close all;

load('mnist_train.mat');
ind = randperm(size(train_X,1));
num = 500;
train_X = train_X(ind(1:num),:);
train_labels = train_labels(ind(1:num));

% Set parameters
no_dims = 2;
initial_dims = 30;
perplexity = 20;
mappedX=tsne(train_X, [], no_dims, initial_dims, perplexity);

%
%gscatter(mappedX(:,1), mappedX(:,2),train_labels);
%figure;
gscatter(mappedX(:,1), mappedX(:,2));

%EOF.
