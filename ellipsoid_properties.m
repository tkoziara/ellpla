function [mass, volume, E0, J0] =  ellipsoid_properties (density, a, b, c, R)

  volume = (4./3.) * pi * a * b * c;
  mass = density * volume;
  
  Inertia = [mass * (b*b + c*c) / 5.0, 0, 0;
             0, mass * (a*a + c*c) / 5.0, 0;
             0, 0, mass * (b*b + a*a) / 5.0];
  
  % note that Inertia = Trace(Euler)*Identity - Euler,
  % hence Euler = 0.5 * Trace (Inertia)*Identity - Inertia,
  % as Trace(Inertia) = 3*Trace(Euler) - Trace(Euler)
  
  E0 = 0.5*trace(Inertia)*eye(3,3) - Inertia;
  
  E0 = R*E0*R';
  
  J0 = Inertia;
end