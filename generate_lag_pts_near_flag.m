clear;
% Problem parameters
L = 10.0;           % length of computational domain (um)
L_s = 6.0;          % Length of swimmer (um)
N = 128;            % number of Cartesian grid points
h = L/N;            % Cartesian grid spacing
rho = 1.0e-12;      % fluid density (g/um^3)
mun = 1.0e-6;       % fluid viscosity (g/(um*s))
mup = 0.5*mun;      % Polymeric viscosity (g/(um*s))

xl = 0;             % lower limit for lag density
xu = 6;             % Upper limit for lag density
sq_w = 1;           % width of square for lag points
num_h = 0.25;        % spacing for z direction for lag points
num_w = 15;         % number of points in width direction
npts = 201;         % number of Lagrangian mesh nodes
ds = L_s/(npts-1);  % Lagrangian meshwidth

vertex = zeros(3,num_w*num_w*length(xl:num_h:xu));
% for i = 1:length(xl:num_h:xu)
%     ih = xl+(i-1)*num_h;
%     [x,y,z] = meshgrid(
xgrid = linspace(-sq_w/2,sq_w/2,num_w);
ygrid = xgrid;
zgrid = xl:num_h:xu;
[x,y,z] = meshgrid(xgrid,ygrid,zgrid);
vertex(1,:) = x(:);
vertex(2,:) = y(:);
vertex(3,:) = z(:);

%%
figure(1); clf; hold on;
for i = 1:length(z(1,1,:))
    plot3(x(:,:,i),y(:,:,i),z(:,:,i),'x');
end

vertex_fid = fopen(['Vertex_pts_', num2str(length(vertex(1,:))), '.vertex'], 'w');
fprintf(vertex_fid, '%d\n',length(vertex(1,:)));
for i = 1:length(vertex(1,:))
    fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', vertex(1,i), vertex(2,i), vertex(3,i));
end
fclose(vertex_fid);