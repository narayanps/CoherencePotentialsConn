function [C, eve_id, peak_index, long_seg, long_seg_chop] =  align_and_correlate_v2(data, chan_features,window, Fs, flags)

%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%

if isempty(window)
    window = 1;
end

window_samples = ceil(window*Fs);
num_elecs=size(data,1);
pnts=size(data,2);
count=1;
eve_id = [];
peak_index=[];

for ii=1:size(chan_features,1)
    features = chan_features(ii,:);
    empty_cells = cellfun(@isempty,features);
    if sum(empty_cells) < num_elecs
        for j=1:num_elecs
            if ~isempty(features{j})
                num_eve = length(features{j}.start_index);
                if ~isempty(num_eve)
                    for k=1:num_eve
                        s=features{j}.start_index(k);
                        e=features{j}.end_index(k);
                        if (s > window_samples && e+window_samples <= pnts)
                            if flags(ii)==0
                                seg{count} = -1*data(j,s:e);
                                long_seg{count} = -1*data(j,s-window_samples:e+window_samples);
                            elseif flags(ii)==1
                                seg{count} = data(j,s:e);
                                long_seg{count} = data(j,s-window_samples:e+window_samples);
                            end
                            count=count+1;
                            eve_id = [eve_id features{j}.event_id(k)];
                            peak_index = [peak_index features{j}.peak_index(k)];
                            
                        end
                    end
                end
            end
        end
        
    else
        disp('None of the channels have events detected according to the threshold specified')
        C=[];
        eve_id=[];
        peak_index=[];
        long_seg=[];
        long_seg_chop=[];
        return;
    end
end
    
    
    
    ne=length(seg);
    count=1;
    if ne>1
        for i=1:ne-1
            seg_1=seg{i};
            m1=find(seg_1==min(seg_1));
            for j=i+1:ne
                seg_2=seg{j};
                m2=find(seg_2==min(seg_2));
                leading=min(m1(1), m2(1));
                falling = min(length(seg_1)-m1(1), length(seg_2)-m2(1));
                err=0;
                if (((m1(1)-leading+1+window_samples)<1) || (m1(1)+falling+window_samples)>length(long_seg{i}))
                    err=1;
                end
                if (((m2(1)-leading+1+window_samples)<1) || (m2(1)+falling+window_samples)>length(long_seg{j}))
                    err=1;
                end
                if err==0
                    
                    seg11=long_seg{i}((m1(1)-leading+1+window_samples):(m1(1)+falling+window_samples));
                    seg22=long_seg{j}((m2(1)-leading+1+window_samples):(m2(1)+falling+window_samples));
                    C{count,1}=eve_id{i};
                    C{count,2}=eve_id{j};
                    C{count,3}=corr(seg11',seg22');
                    count=count+1;
                    long_seg_chop{i,j}=[seg11; seg22];
                    
                end
                
            end
        end
    else
        C = [];
    end





