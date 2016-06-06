function [con_point, depth, A, B] = ellipsoid_plane (x, U1, R1, point, normal)
  
  % IMPORTANT: I am assuming that x, point, normal are column vectors here
  
  % what I think you need to do is simply move the plane into the natural
  % coordinate system of our input ellipsoid; it is centered at x and
  % oriented according to the R1 rotation, so I think that we need to do
  % something like:
  
  point = R1'*(point - x); % so that plane's point, relative to the origin of
                          % the ellipsoid is suitably shifted and rotated
                          % and this needs to be followed by

  normal = R1' * normal;   % to redefine the normal in the natural coordinate
                          % system of the ellipsoid
                          
  % now we extract a, b, c from U1 as follows
  a = 1/sqrt(U1(1,1));
  b = 1/sqrt(U1(2,2));
  c = 1/sqrt(U1(3,3));
  C1 = [1/a, 0, 0; 0, 1/b, 0; 0, 0, 1/c];
   
  v = rand(3,1);        % let's generate a random vector first
  r = cross(normal, v); % now r'*normal = 0 as required
  r = r / norm(r);      % lets normalize it (make it's length = 1)
  s = cross(normal, r); % now s'*normal = 0 and s'*r = 0 as required#
  q = point;
   
  % now lets check for condition 7 to be fullfilled 
  
  % by doing this transformation we are already satisfying condition 4-6 as
  % rr and ss are orthogonal to each other and to the normal vector
  
  % now lets choose omega and fulfill condition 7
  
  % here I will define all the dot products needed
  D1 = dot((C1*r),(C1*s));
  D2 = dot((C1*r),(C1*r));
  D3 = dot((C1*s),(C1*s));
  
  if (abs(D2-D3) < 1E-10)
    omega = 0.25*pi;
  else
    omega = 0.5 * atan((2*D1)/(D2-D3));
  end
  
  rr = cos(omega)* r + sin(omega) * s;  % here we are transforming r and s to fullfill our condition 7
  ss = -sin(omega) * r + cos(omega) * s;
  %D7 = dot((C1*rr),(C1*ss));
  
  r = rr;
  s = ss;
  
  D1 = dot((C1*r),(C1*s));
  D2 = dot((C1*r),(C1*r));
  D3 = dot((C1*s),(C1*s));
  D4 = dot((C1*q),(C1*r));
  D5 = dot((C1*q),(C1*s));
  D6 = dot((C1*q),(C1*q));
  
  % now lets find the ellipse of the intersection
  d = D6 -((D4)^2/D2)-((D5)^2/D3); % this is just a value
  
  % semi-axes of the ellipse
  Aa = (1-d)/D2;
  Bb = (1-d)/D3;
  
  if Aa > 0. && Bb > 0.  
    A  = sqrt(Aa);
    B  = sqrt(Bb);
  
    % now the centre of the ellipse
    t0 = -D4/D2;
    u0 = -D5/D3;
    x0 = q + t0*r + u0*s;
  
    % now lets find the contact point
    alpha = sqrt(1/((x0(1)^2/a^2) + (x0(2)^2/b^2) + (x0(3)^2/c^2)));
    con_point = alpha * x0;
    
    % the penetration depth (ellipsoid coodinates)
    depth = dot(normal,  x0 - con_point);
    
    % input coodinates
    con_point = R1*con_point + x;
  else
    con_point = [0, 0, 0];
    depth = -1;
    A = 0;
    B = 0;
  end 
end