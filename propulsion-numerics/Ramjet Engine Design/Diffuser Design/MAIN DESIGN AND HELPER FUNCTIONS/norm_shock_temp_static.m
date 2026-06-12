function [t_ratio] = norm_shock_temp_static(M)
gamma=1.4;
t_ratio=(2*gamma*M.^2-(gamma-1)).*(2+(gamma-1)*M.^2)./((gamma+1)^2.*M.^2);
end
