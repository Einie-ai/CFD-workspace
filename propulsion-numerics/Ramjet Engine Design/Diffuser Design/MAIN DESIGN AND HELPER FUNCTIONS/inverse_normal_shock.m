function [Ms] = inverse_normal_shock(static_pressure_ratio)
gamma=1.4;
Ms=((static_pressure_ratio*(gamma+1)+(gamma-1))/(2.0*gamma))^(0.5);
end