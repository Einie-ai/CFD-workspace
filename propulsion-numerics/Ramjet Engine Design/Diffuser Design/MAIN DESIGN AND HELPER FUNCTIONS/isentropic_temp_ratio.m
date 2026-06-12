function [t_ratio]= isentropic_temp_ratio(M)
gamma=1.4;
t_ratio=(1+((gamma-1)/2)*M.^2).^(gamma/(gamma-1));
end