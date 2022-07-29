function [T, Y] = MoonShipModel1(t_vector, y0, u1, u2)
%t_vector - time
%y0 - starting conditions
%u1, u2 - control

  [T, Y] = ode45(@ShipModel, t_vector, y0);
  function dxdt = ShipModel(t, x)
% x(1) = V
% x(2) = theta
% x(3) = h
% x(4) = phi
% x(5) = m

dxdt = zeros(5,1);

P0 = 440;
R0 = 1738400;
g1 = 9.80665;
g0 = 1.623;
Pud = 319;

%u1 = 0;
%u2 = 80;
g = g0 .* (R0./(R0 + x(3))) .^ 2;
R = R0 + x(3);
W = g1.*(P0+u2) ./ x(5);

dxdt(1) = W.*cos(u1 - x(2)) - g.*cos(x(2));
dxdt(2) = (1./x(1)).*(W.*sin(u1-x(2)) + g.*sin(x(2)));
dxdt(3) = x(1) .* cos(x(2));
dxdt(4) = (x(1) .* sin(x(2))) ./ R;
dxdt(5) = -(P0 + u2) ./ Pud;
%dxdt = [dx1; dx2; dx3; dx4; dx5];
    end
end
