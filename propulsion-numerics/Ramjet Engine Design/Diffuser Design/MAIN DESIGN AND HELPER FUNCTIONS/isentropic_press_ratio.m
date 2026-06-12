function [p_ratio]= isentropic_press_ratio(M)
gamma=1.4;
p_ratio=(1+((gamma-1)/2)*M.^2).^(gamma/(gamma-1));
end