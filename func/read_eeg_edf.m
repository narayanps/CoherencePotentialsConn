function [hdr,data] = read_eeg_edf(fname)

%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%

addpath(strcat(pwd, '/external'))
[data,hdr] = lab_read_edf(fname);
end

