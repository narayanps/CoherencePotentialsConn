function [A_iei, A_cf, a] = get_clust_conn_matrix(clust, chan_labels,Fs, min_ne, total_eve_per_chan)

%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%



count=0;
ds=[];
d=0;
N=length(chan_labels);
for k=1:length(clust)
    if length(clust{k,1}.events) >= min_ne
        tmp=split(clust{k,1}.event_id,'_');
        clust{k,1}.chan = tmp(:,:,1);
        chans=squeeze(tmp(:,:,1));
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
    
end
if ~isempty(ds)
    [unique_chans_pairs,~,idx] = unique(ds.chan_pairs,'rows');
    a=[];
    for i=1:1:max(idx)
        id =find(idx==i);
        a.connectStrength(i)=mean(ds.conn_strength(id));
        a.connectFreq(i)=length(id)/(total_eve_per_chan(unique_chans_pairs(i,1))*total_eve_per_chan(unique_chans_pairs(i,2)));
        a.chanPairs(i,:)=unique_chans_pairs(i,:);
    end
else
    a=[];
end
A_iei = nan(N,N);
A_cf = nan(N,N);
if ~isempty(a)
for i=1:size(a.chanPairs,1)
        A_iei(a.chanPairs(i,1), a.chanPairs(i,2))  = a.connectStrength(i); 
        A_cf(a.chanPairs(i,1), a.chanPairs(i,2))  = a.connectFreq(i); 
end
end

end

