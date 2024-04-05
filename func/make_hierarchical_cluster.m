function [ Z, I, cophenetic_CC] = make_hierarchical_cluster(corr_vec, method, plot_dend, no_of_leaves)
%AUTHOR : Narayan Subramaniyam

%This function performs agglomerative hierarchical clustering.
%INPUTS

% 1) corr_vec = lower / upper triangular entries of correlation matrix
% excluding diagonals in vector form. Size : mX1

% 2) method = method to compute linkage : 'average'(DEFAULT), 'weighted', 'complete',
%etc

% 3) c = cut-off value to make clusters
% 4) crit = criterion : 'distance' or 'inconsistency'

% 5) plot_dend = 1(yes) / 0 (no)
% 6) no_of_leaves = no_of_leaves to show in dendrogram

%OUTPUTS

% no_of_clusters = No of clusters formed
% clusters = struct providing info on how events are clustered
% I = Inconsistency matrix
% cophenetic_CC  = cophenetic correlation coefficient



% some mandatory checks to make sure inputs are correct
[~,n] = size(corr_vec);
if n>1
   error('Input must be a row vector!')
end

if nargin < 2
    
    method = 'average';
end

if nargin <3
    plot_dend = 0;
end

if nargin <4
    if plot_dend==1
        no_of_leaves = 0;
    end
end

%compue distance matrix
Y = 1-abs(corr_vec);

%get linkage information
Z = linkage(Y',method);

%plot dendrogram
if plot_dend==1
   dendrogram(Z,no_of_leaves)
end

%Verification of cluster tree : compute cophenetic correlation coefficient
%and inconsistency matrix
cophenetic_CC = cophenet(Z,Y);
I = inconsistent(Z);

% %cluster based on given cut-off and criterion
% T = cluster(Z,'cutoff',c, 'criterion', crit);
% no_of_clust = length(unique(T));
% 
% % form a cell with clustered events
% for i=1:1:no_of_clust
%     clusters{1,i}.events = find(T==i);
% end