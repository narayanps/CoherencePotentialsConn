function [f] = get_features(data, abs_sd, flag, fact, offset, Fs, ch_name)

%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%
if isempty(offset)
    offset=0.5;
end
peak_thr = -fact*abs_sd;
if flag==1
    data = data;
elseif flag==0
    data = -1*data;
else
    error('flag must be 1 or 0 or left empty (defaults to 1)')
end
offset_samples = offset*Fs;
f=[];
I =find(data < peak_thr);
if ~isempty(I)
    B=[];
    for j=offset_samples:length(data)-offset_samples
        if (data(j)>0 && data(j+1)<0) ||(data(j)<0 && data(j+1)>0)
            B = [B j];
        end
    end
    if length(I) > 0
        ind = find(diff(I)>1);
        if ~isempty(ind)
            ind = [ind ind(end)+1];
        else
            ind = length(I);
        end
        trigs=B;
        ends=B;
        c=0;
        if length(ind)>0
            for k=1:length(ind)
                t=find(trigs<I(ind(k)));
                e=find(ends>=I(ind(k)));
                
                if ~( isempty(t) || isempty(e))
                    t=trigs(t(length(t)));
                    e=ends(e(1));
                    c=c+1;
                    start_index(c) = t;
                    end_index(c) = e;
                else
                    continue
                end
                
            end
        end
    end
    if c>0
        start_index=unique(start_index, 'stable');
        end_index = unique(end_index, 'stable');
        for j=1:length(start_index)
            seg = data(start_index(j):end_index(j));
            [~,id] = nanmin(seg);
            peak_index(j) = start_index(j)+id-1;
        end
        
        f.start_index=start_index;
        f.end_index=end_index;
        f.peak_index=peak_index;
        f.event_id = arrayfun(@(i) sprintf(strcat(ch_name,'_%d'), i), 1:length(start_index), 'UniformOutput', false);
        
    end
else
f=[];
end
end

