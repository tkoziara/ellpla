% clear all
clear;

% ellipsoid meshing subdivisions
ndiv = 20;

% animation FPS
fps = 120;

% initial ellipsoid
X0 = [0; 0; 0.1];
a=0.02;
b=0.02;
c=0.05;
%R=eye(3,3);
R=expmap([0,pi/3.,0]);

% initial matrix representation of it
U=[1/a^2, 0, 0; 0, 1/b^2, 0; 0, 0, 1/c^2];
K = R*U*R';

% material properties
density = 0.5E3;
young = 1E6;
poisson = 0.3;

% mass properties from geometry and density
[mass, volume, E0, J0] = ellipsoid_properties (density, a, b, c, R);

% time integration step and duration
step = 0.0001;
duration = 1.0;

% gravity
grav = [0; 0; -10];

% contact damper constant
damper = 10;

% contact plane
normal = [0; 0; 1];
point = [0; 0; 0];

% referential point whose time history we are going to track
P0 = R*[0; 0; c];
P0_position = zeros(3,duration/step);
P0_velocity = zeros(3,duration/step);
P0_time = zeros(1,duration/step);

Q0 = R*[0; 0; -c];
P0Q0_length = zeros(1,duration/step);

% inverse of Euler tensor
iE0 = inv(E0);

% effective value of Young modulus!!!
EE0 = young/(1-poisson^2);

x = X0;         % current position of the mass centre
v = [0; 0; 0];  % initial velocity

F = eye(3,3);   % initial deformation gradient
vF= zeros(3,3); % initial deformation gradient velocity

% current time
time = 0.0;

% animation time history
ani_x = zeros(3,duration*fps+1);
ani_K = zeros(3,3,duration*fps+1);
ani_time = 0;
ani_step = 1/fps;
ani_frame = 1;
ani_x(:,1) = x;
ani_K(:,:,1) = K;

% rune simulation first
tic
nstep = 1;
while (time < duration)
     
    % update confirguration: first half-step
    x = x + 0.5*step*v;
    F = F + 0.5*step*vF;
    
    % evalute new ellipsoid representation
    iF = inv(F);
    K1 = iF'*K*iF;
    [R, U] = eig(K1);
  
    % rigid motion force
    f = mass*grav;
  
    % pseudo-rigid motion force
    M0 = zeros(3,3);
  
    % detect contact
    [con_point, depth, A, B] = ellipsoid_plane (x, U, R, point, normal);
    
    if  depth > 0.0                     % negative gap => penetration
        X = inv(F)*(con_point - x) + X0; % referential image of the contact point
        vX = vF*(X-X0) + v;
        vn = dot(vX,normal);
        % radii of curvature
        r1 = A^2 / depth;
        r2 = B^2 / depth;
        % effective radius of curvature
        eff_rad = sqrt(r1*r2);
        % contact force value
        cfv = 4/3* EE0 * sqrt(eff_rad) * (depth^1.5) - damper*vn;
        % contat force
        cf = cfv * normal;
        % contact moment
        CM = cf*(X-X0)';
        % accumulate contact force and moment
        f = f + cf;
        M0 = M0 + CM;
    end
    
      
    % update velocity
    v = v + step*f/mass;
    vF = vF + step*(M0 - volume*piola(young,poisson,F))*iE0;
  
    % update confirguration: second half-step
    x = x + 0.5*step*v;
    F = F + 0.5*step*vF;
   
    % increase time
    time = time + step;
    
    if time > (ani_time + ani_step)
        ani_frame = ani_frame + 1;
        ani_time = time;
        ani_x(:,ani_frame) = x;
        ani_K(:,:,ani_frame) = K1;
    end
    
    % update P0 time histories
    P0_position(:,nstep) = F*(P0-X0)+x;
    P0_velocity(:,nstep) = vF*(P0-X0)+v;
    P0_time(nstep) = time;
    
    % update Q0 time histories
    P0Q0_length(nstep) = norm(F*(P0-Q0));
    
    % update step counter
    nstep = nstep + 1;
end
toc

% plot P0 time histories
figure(1);
plot (P0_time, P0_position(1,:), P0_time, P0_position(2,:), '--', P0_time, P0_position(3,:), '-.')
title ('P0 displacement')
xlabel ('Time [s]')
ylabel ('Displacement [m]')
legend ('x_1(P0)', 'x_2(P0)', 'x_3(P0)')

figure(2);
plot (P0_time, P0_velocity(1,:), P0_time, P0_velocity(2,:), '--', P0_time, P0_velocity(3,:), '-.')
title ('P0 velocity')
xlabel ('Time [s]')
ylabel ('Velocity [m/s]')
legend ('v_1(P0)', 'v_2(P0)', 'v_3(P0)')

figure(3);
plot (P0_time, P0Q0_length)
title ('|P0-Q0| length')
xlabel ('Time [s]')
ylabel ('Diameter [m]')

% set up 3d view (copied from internet - read on)
figure(4);
view(3);
axis equal
axis tight
camlight
q = 2*max([a, b, c]);
axis([-q, q, -q, q, 0, q]); % user specified view size

% initial ellipsoid
x = ani_x(:,1);
K1 = ani_K(:,:,1);

% render initial ellipsoid patch
rmax = max([a, b, c]);
[eX,eY,eZ] = ndgrid(linspace(X0(1)-rmax, X0(1)+rmax,ndiv), ...
                    linspace(X0(2)-rmax, X0(2)+rmax,ndiv), ...
                    linspace(X0(3)-rmax, X0(3)+rmax,ndiv));
eV = zeros(size(eX));
for i=1:size(eX,1)
  for j = 1:size(eX,2)
    for k = 1:size(eX,3)
      x1 = eX(i, j, k);
      y1 = eY(i, j, k);
      z1 = eZ(i, j, k);
      w1 = [x1; y1; z1] - x;
      eV(i, j, k) = w1'* K * w1 - 1;
    end
  end
end

% create graphical representation of the ellipsoid
p = patch(isosurface(eX,eY,eZ,eV,0));
set(p,'FaceColor','b','EdgeColor','b','FaceAlpha',1.0);
 
% play animation
last_frame = ani_frame;
ani_frame = 1;
while (ani_frame < last_frame)
    
   % read ellipsoid
   x = ani_x(:,ani_frame);
   K = ani_K(:,:,ani_frame);
    
    % redraw the ellipsoid
    rmax = sqrt(1/min(spdiags(U)));
    [eX,eY,eZ] = ndgrid(linspace(x(1)-rmax, x(1)+rmax,ndiv), ...
                    linspace(x(2)-rmax, x(2)+rmax,ndiv), ...
                    linspace(x(3)-rmax, x(3)+rmax,ndiv));                  
    eV = zeros(size(eX));
    for i=1:size(eX,1)
      for j = 1:size(eX,2)
        for k = 1:size(eX,3)
          x1 = eX(i, j, k);
          y1 = eY(i, j, k);
          z1 = eZ(i, j, k);
          w1 = [x1; y1; z1] - x;
          eV(i, j, k) = w1'* K * w1 - 1;
        end
      end
    end
  
   % delete previous patch to create animation effect
   delete (p);
   p = patch(isosurface(eX,eY,eZ,eV,0));
   set(p,'FaceColor','b','EdgeColor','b','FaceAlpha',1.0);
 
   % draw and pause to create an effect of animation
   drawnow;
   pause(ani_step);
   
   % go to next frame
   ani_frame = ani_frame + 1;
end
