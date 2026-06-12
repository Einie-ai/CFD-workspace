function [res] = area_mach_func(M)
gamma=1.4;
res=(1./M).*((2/(gamma+1)).*(1+((gamma-1)/2).*M.^2)).^((gamma+1)/(2*(gamma-1)));
end