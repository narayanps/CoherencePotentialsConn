function [abs_sd] = get_sd(data)
%%%%%%%%AUTHOR : NARAYAN P SUBRAMANIYAM%%%%%%%%%%%%
%%%%%%%%%%%%%2023 \ SAPIEN LABS%%%%%%%%%%%%%%%%%%%%%
abs_sd = mean(std(data,[],2));
