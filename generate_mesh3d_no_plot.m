%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% Problem parameters
L = 10.0;         % length of computational domain (um)
L_s = 6.0;        % Length of swimmer (um)
N = 128;           % number of Cartesian grid points
h = L/N;          % Cartesian grid spacing
rho = 1.0e-12;        % fluid density (g/um^3)
mu = 0.01e-4;        % fluid viscosity (g/(um*s))

r0 = 0.2067;        % radius of unstressed rod (um)
w = -2*pi/2.1361;   % wave number
k = 2.0;            % decay of helix
npts = 101;         % number of Lagrangian mesh nodes
ds = L_s/(npts-1);  % Lagrangian meshwidth

a1 = 3.5e-3; % Bending modulus (g um^3/s^2)
a2 = 3.5e-3; % Bending modulus (g um^3/s^2)
a3 = 1.0e-4; % Twist modluls   (g um^5/s^2)
a1_hook = a1*0.01;
a2_hook = a2*0.01;
a3_hook = a3*1.0;


b1 = 8.0e-1; % Shear modulus   (g um/s^2)
b2 = 8.0e-1; % Shear modulus   (g um/s^2)
b3 = 8.0e-1; % Stretch modulus (g um/s^2)

kappa2 = 1.3057; % Intrinsic curvature (1/um)
kappa1 = 0.0;    % Intrinsic curvature (1/um)
tau1 = -2.1475;  % Right-handed intrinsic twist (1/um)
tau2 = 2.1475;   % Left-handed intrinisic twist (1/um)
kappa1_hook = 0.0; % Hook Intrinsic curvature (1/um)
kappa2_hook = 0.0; % Hook intrinsic curvature (1/um)
tau1_hook = 0.0; % hook right-handed intrinsic twist (1/um)
tau2_hook = 0.0; % hook left-handed intrinsic twist (1/um)

gamma = 1.0e-3; % Twist gradient coefficient (g^(1/2) um^(3/2)/s)

nhook = 2;  % Length of hook

%%
% Initial Computations
s = ds*(0:npts);
s0 = s(nhook); % Location of last point of hook
%(x,y,z) points of hook
x = zeros(1,npts); y = x; z = x;
z(1:nhook) = ds*(0:nhook-1);
r = zeros(1,npts); dr = r; ddr = r;
z0 = z(nhook);
for i = nhook:npts
    r(i) = r0*(1-exp(-k*(s(i) - s0).^2));
    dr(i) = r0*2.0*k*(s(i)-s0)*exp(-k*(s(i)-s0).^2);
    ddr(i) = r0*2.0*k*exp(-k*(s(i)-s0).^2)*(1.0-2.0*k*(s(i)-s0).^2);
end

% construct flagellum
for i = nhook+1:npts
    z(i) = z(i-1) + sqrt((1.0-dr(i).^2)./(1.0+(w*r(i)).^2))*ds*0.5+sqrt((1.0-dr(i-1).^2)./(1.0+(w*r(i-1)).^2))*ds*0.5;
    x(i) = r(i)*cos(w*(z(i)-z0));
    y(i) = r(i)*sin(w*(z(i)-z0));
end

% Compute orthonormal triads
D3 = zeros(npts,3); D2 = D3; D1 = D3;
for i = 1:nhook
    D1(i,1) = 1.0;
    D1(i,2) = 0.0;
    D1(i,3) = 0.0;
    
    D2(i,1) = 0.0;
    D2(i,2) = 1.0;
    D2(i,3) = 0.0;
    
    D3(i,1) = 0.0;
    D3(i,2) = 0.0;
    D3(i,3) = 1.0;
end

for i = nhook+1:npts
    % D3 is unit tangent vector
    dz = sqrt((1.0-dr(i)^2)/(1.0+(w*r(i))^2));
    dx = dr(i)*cos(w*(z(i)-z0)) - r(i)*w*sin(w*(z(i)-z0))*dz;
    dy = dr(i)*sin(w*(z(i)-z0)) + r(i)*w*cos(w*(z(i)-z0))*dz;
    den = sqrt(dx^2+dy^2+dz^2);
    D3(i,1) = dx/den;
    D3(i,2) = dy/den;
    D3(i,3) = dz/den;
    
    % D1 is unit normal vector
    ddz = -w^2*r(i)*dr(i)*sqrt(1-dr(i)^2)/sqrt(1.0+(w*r(i))^2)^3-dr(i)*ddr(i)/sqrt(1.0-dr(i)^2)/sqrt(1.0+(w*r(i))^2);
    ddx = ddr(i)*cos(w*(z(i)-z0)) - 2.0*dr(i)*sin(w*(z(i)-z0))*w*dz - r(i)*cos(w*(z(i)-z0))*(w*dz)^2 - r(i)*sin(w*(z(i)-z0))*w*ddz;
    ddy = ddr(i)*sin(w*(z(i)-z0)) + 2.0*dr(i)*cos(w*(z(i)-z0))*w*dz - r(i)*sin(w*(z(i)-z0))*(w*dz)^2 - r(i)*cos(w*(z(i)-z0))*w*ddz;
    
    TT = dx*ddx+dy*ddy+dz*ddz;
    
    DtDs(1) = (ddx*den^2-dx*TT)/den^3;
    DtDs(2) = (ddy*den^2-dy*TT)/den^3;
    DtDs(3) = (ddz*den^2-dz*TT)/den^3;
    
    denD3 = sqrt(sum(DtDs.^2));
    D1(i,1) = DtDs(1)/denD3;
    D1(i,2) = DtDs(2)/denD3;
    D1(i,3) = DtDs(3)/denD3;
    
    % D2 is unit binormal vector
    D2(i,1) = D3(i,2)*D1(i,3) - D3(i,3)*D1(i,2);
    D2(i,2) = D3(i,3)*D1(i,1) - D3(i,1)*D1(i,3);
    D2(i,3) = D3(i,1)*D1(i,2) - D3(i,2)*D1(i,1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
vertex_fid = fopen(['Higdon_helix_' num2str(npts) '.vertex'], 'w');
rod_fid = fopen(['Higdon_helix_' num2str(npts) '.rod'], 'w');
director_fid = fopen(['Higdon_helix_' num2str(npts) '.director'], 'w');
target_fid = fopen(['Higdon_helix_' num2str(npts) '.anchor'], 'w');

fprintf(vertex_fid, '%d\n', npts);
fprintf(rod_fid, '%d\n', npts-1);
fprintf(director_fid, '%d\n', npts);
fprintf(target_fid, '%d\n', 1);

fprintf(target_fid, '%d %1.16e\n', 0);

for r = 0:nhook-1
    fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', x(r+1), y(r+1), z(r+1));
    fprintf(rod_fid, ['%6d %6d %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e ' ...
                    '%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n'], r, r+1, ds, ...
          a1_hook, a2_hook, a3_hook, b1, b2, b3, kappa1_hook, kappa2_hook, tau1_hook, tau2_hook, gamma, 1.0);
    fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D1(r+1,1), D1(r+1,2), D1(r+1,3));
    fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D2(r+1,1), D2(r+1,2), D2(r+1,3));
    fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D3(r+1,1), D3(r+1,2), D3(r+1,3));
end
for r = nhook:npts-1
  fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', x(r+1), y(r+1), z(r+1));
  if (r ~= npts-1)
  fprintf(rod_fid, ['%6d %6d %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e ' ...
                    '%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n'], r, r+1, ds, ...
          a1, a2, a3, b1, b2, b3, kappa1, kappa2, tau1, tau2, gamma, 0.0);
  end

  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D1(r+1,1), D1(r+1,2), D1(r+1,3));
  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D2(r+1,1), D2(r+1,2), D2(r+1,3));
  fprintf(director_fid, '%1.16e %1.16e %1.16e\n', D3(r+1,1), D3(r+1,2), D3(r+1,3));
end %for

fclose(vertex_fid);
fclose(rod_fid);
fclose(director_fid);
fclose(target_fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
