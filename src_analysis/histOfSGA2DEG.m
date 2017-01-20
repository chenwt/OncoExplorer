% histogram of SGA->DEG edge weight.

figure('color', [1 1 1]);
k = 0;
for can = {'pancan','brca','gbm','ov'}
    k = k + 1;
    subplot(2,2,k);
    filename = strcat('bond_',can,'.txt');
    filename = filename{1};
    bond = load(filename);
    hist(bond,200);
    mean(bond == 1) %+ mean(bond == 2)
    xlabel('sga2deg bond');
    ylabel('frequency');
    title(can);
end

