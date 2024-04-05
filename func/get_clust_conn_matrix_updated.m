function [A_iei, A_iei_min, A_iei_max, A_min_cf, A_all_cf, A_all_cf_norm, A_min_cf_norm, A_min_cf_norm_clust, A_perc_cf, a] = get_clust_conn_matrix_updated(clust, chan_labels,Fs, min_ne, total_eve_per_chan)

%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%
count=0;
N=length(chan_labels);
clust_count = zeros(N,N);
for i=1:N
    for j=i+1:N
        for k=1:length(clust)
            if ~isempty(cell2mat(strfind(clust{k,1}.event_id, chan_labels{i,1})))
                count=count+1;
            end
            if ~isempty(cell2mat(strfind(clust{k,1}.event_id, chan_labels{j,1})))
                count=count+1;
            end
        end
        clust_count(i,j) = count;
        count=0;
    end
end


count=0;
ds=[];
d=0;
for k=1:length(clust)
    ds=[];
    d=0;
    if length(clust{k,1}.events) >= min_ne
        tmp=squeeze(split(clust{k,1}.event_id,'_'));
        clust{k,1}.chan = (tmp(:,1));
        chans=(tmp(:,1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [unique_chans,~,idx1] = unique(chans);
        eve_counts = histc(idx1, unique(idx1));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        count=count+1;
        for m=1:1:length(chans)
            for n=m+1:length(chans)
                if strcmp(chans{m}, chans{n}) ~= 1
                    d=d+1;
                    ds.chan_pairs(d,:) = [find(strcmp(chan_labels, chans{m})), find(strcmp(chan_labels, chans{n}))];
                    ds.conn_strength(d,1) = (double(abs(clust{k,1}.peak_index(m) - clust{k,1}.peak_index(n))))./Fs;
                end
                
            end
        end
        
    end
    if ~isempty(ds)
        ds_sorted = sort(ds.chan_pairs, 2);
        [unique_chans_pairs,~,idx] = unique(ds_sorted,'rows');
        for i=1:size(unique_chans_pairs,1)
            id =find(idx==i);
            a{k}.connectStrength(i)=mean(ds.conn_strength(id));
            a{k}.minconnectStrength(i)=min(ds.conn_strength(id));
            a{k}.maxconnectStrength(i)=max(ds.conn_strength(id));
            a{k}.chanPairs(i,:)=unique_chans_pairs(i,:);
           % norm_=(total_eve_per_chan(unique_chans_pairs(i,1))*total_eve_per_chan(unique_chans_pairs(i,2)));
            %a{k}.connectFreq(i)=length(id)/norm_;
            a{k}.minconnectFreq(i) = min(eve_counts(find(strcmp(unique_chans,chan_labels{unique_chans_pairs(i,1)}))),...
                eve_counts(find(strcmp(unique_chans,chan_labels{unique_chans_pairs(i,2)}))));
            a{k}.allconnectFreq(i) = eve_counts(find(strcmp(unique_chans,chan_labels{unique_chans_pairs(i,1)}))) * ...
                eve_counts(find(strcmp(unique_chans,chan_labels{unique_chans_pairs(i,2)})));
            
            %a{k}.minconnectFreqNorm(i) = a{k}.minconnectFreq(i) / norm_;
            %a{k}.allconnectFreqNorm(i) = a{k}.allconnectFreq(i) / norm_;


        end
    else
        a{k}=[];
    end
end
% if ~isempty(ds)
%     [unique_chans_pairs,~,idx] = unique(ds.chan_pairs,'rows');
%     a=[];
%     for i=1:1:max(idx)
%         id =find(idx==i);
%         a.connectStrength(i)=mean(ds.conn_strength(id));
%         a.connectFreq(i)=length(id)/(total_eve_per_chan(unique_chans_pairs(i,1))*total_eve_per_chan(unique_chans_pairs(i,2)));
%         a.chanPairs(i,:)=unique_chans_pairs(i,:);
%     end
% else
%     a=[];
% end
num_clust = length(a);
A_iei_clust = nan(num_clust,N,N);
A_iei_min_clust = nan(num_clust,N,N);
A_iei_max_clust = nan(num_clust,N,N);

A_min_cf_clust = zeros(num_clust,N,N);
A_all_cf_clust = zeros(num_clust,N,N);



for k=1:length(a)
    if ~isempty(a{k})
        for i=1:size(a{k}.chanPairs,1)
            A_iei_clust(k,a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.connectStrength(i);
            A_iei_min_clust(k,a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.minconnectStrength(i);
            A_iei_max_clust(k,a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.maxconnectStrength(i);
            A_min_cf_clust(k, a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.minconnectFreq(i);
            A_all_cf_clust(k, a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.allconnectFreq(i);
           % A_norm_min_cf_clust(k, a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.minconnectFreqNorm(i);
           % A_norm_all_cf_clust(k, a{k}.chanPairs(i,1), a{k}.chanPairs(i,2))  = a{k}.allconnectFreqNorm(i);
        end
    end
end

A_min_cf = squeeze(sum(A_min_cf_clust,1));
A_all_cf = squeeze(sum(A_all_cf_clust,1));

A_iei = squeeze(nanmean(A_iei_clust,1));
A_iei_min = squeeze(nanmean(A_iei_min_clust,1));
A_iei_max = squeeze(nanmean(A_iei_max_clust,1));

A_perc_cf=zeros(14,14);
for i=1:14
    for j=i+1:14
        A_perc_cf(i,j) = length(find(A_all_cf_clust(:,i,j)>=1))./num_clust;
    end
end
A_all_cf_norm=zeros(14,14);
A_min_cf_norm=zeros(14,14);
A_min_cf_norm_clust=zeros(14,14);
for i=1:14
    for j=i+1:14
        A_all_cf_norm(i,j) = A_all_cf(i,j)./(1*total_eve_per_chan(i)*total_eve_per_chan(j));
        A_min_cf_norm_clust(i,j) = A_min_cf(i,j)./(clust_count(i,j));
        A_min_cf_norm(i,j) = A_min_cf(i,j)./(min(total_eve_per_chan(i),total_eve_per_chan(j)));
    end
end


end



