function [p_ratio] = norm_shock_press_stag(M)
gamma=1.4;
p_ratio=((gamma+1)*M.^2/(2+(gamma-1)*M.^2)).^(gamma/(gamma-1))*((gamma+1)/(2*gamma*M.^2-(gamma-1))).^(1/(gamma-1));
end