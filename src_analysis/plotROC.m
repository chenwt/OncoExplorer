% plotROC.m

% show_roc.m

% Comparison of parameters.
clc; clear all; close all;

n1 = 'eps=1e-4';
n2 = 'eps=1e-5';
n3 = 'eps=1e-6';
n4 = 'eps=1e-7';
n5 = 'alph=0.05';
roc1 = load(n1);
% % roc2 = load(n2);
% roc3 = load(n3);
% roc4 = load(n4);
% roc5 = load(n5);
figure('color',[1 1 1]);
hold on;
plot(roc1(:,2),roc1(:,3),'LineWidth',1.5,'Marker','o', 'MarkerFaceColor','k',...
    'MarkerEdgeColor','k','Color','k');
plot(roc2(:,2),roc2(:,3),'LineWidth',1.5,'Marker','o', 'MarkerFaceColor','r',...
    'MarkerEdgeColor','k','Color','k');
% plot(roc3(:,2),roc3(:,3),'LineWidth',1.5,'Marker','o', 'MarkerFaceColor','g',...
%     'MarkerEdgeColor','k','Color','k');
% plot(roc4(:,2),roc4(:,3),'LineWidth',1.5,'Marker','o', 'MarkerFaceColor','b',...
%     'MarkerEdgeColor','k','Color','k');
% plot(roc5(:,2),roc5(:,3),'LineWidth',1.5,'Marker','s', 'MarkerFaceColor','r',...
%     'MarkerEdgeColor','k','Color','k');
xlabel('precision');
ylabel('recall');
legend(n1,n2,n3,n4,n5);
%title('')
% 
xlim([0 1])
ylim([0 1])
box on