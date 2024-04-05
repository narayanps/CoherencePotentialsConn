function [true_diff,pval ]= perm_test(data_1,data_2)
if size(data_1,1) > 1
    data_1=data_1';
end

if size(data_2,1) > 1
    data_2=data_2';
end
ind = ~isnan(data_1);
data_1=data_1(ind);
ind = ~isnan(data_2);
data_2=data_2(ind);
all_data = cat(1,data_1', data_2');
n1 = size(data_1,2);
n2 = size(data_2,2);
labels = cat(1,ones(n1,1), 2*ones(n2,1));
true_diff = (mean(all_data(labels==1)) - mean(all_data(labels==2)));
%[~,~,true_diff] = kstest2(all_data(labels==1), all_data(labels==2));
iter=500;
perm_vals = zeros(iter,1);

for j=1:iter
    shuffle_labels = labels(randperm(n1+n2));
    perm_vals(j) = (mean(all_data(shuffle_labels==1)) - mean(all_data(shuffle_labels==2)));
%    [~,~,perm_vals(j)] = kstest2(all_data(shuffle_labels==1), all_data(shuffle_labels==2));
end

perm_mean = mean((perm_vals));
perm_std = std((perm_vals));
zdist = (true_diff - perm_mean) / perm_std;
pval = (1-normcdf(abs(zdist)));

pval = sum(abs(perm_vals) > abs(true_diff))/iter;
% 
%histogram(perm_vals,40)
%hold on
%plot([1 1]*true_diff, get(gca, 'ylim'), 'r--', 'LineWidth', 3)
%legend({'shuffled', 'true'})