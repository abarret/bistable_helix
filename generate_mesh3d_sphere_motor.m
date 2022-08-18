%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% Problem parameters
L = 10;
N = 128;
dx = L/N;
sphere_radius = 2*dx;
%%
addpath('distmesh');
% Initial Computations
[p,t] = distmeshnd(@(p) sqrt(sum(p.^2,2))-sphere_radius,@huniform,dx,[-1,-1,-1;1,1,1],[]);
npts = length(p)
for i=1:length(t)
    idx = (i-1)*4+1;
    for j=1:4
        edge(1,idx) = t(i,j)-1;
        if(j ~= 4)
            edge(2,idx) = t(i,j+1)-1;
        else
            edge(2,idx) = t(i,1)-1;
        end
        idx = idx+1;
    end
end
edge = sort(edge);
edge = unique(edge','rows')';
tot_edge = length(edge);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
vertex_fid = fopen(['Sphere_motor_' num2str(npts) '.vertex'], 'w');
target_fid = fopen(['Sphere_motor_' num2str(npts) '.target'], 'w');
spring_fid = fopen(['Sphere_motor_' num2str(npts) '.spring'], 'w');

fprintf(vertex_fid, '%d\n', npts);
fprintf(spring_fid, '%d\n', tot_edge);
fprintf(target_fid, '%d\n', npts);

for r = 0:npts-1
  fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', p(r+1,1), p(r+1,2), p(r+1,3));
  fprintf(target_fid, '%d %1.16e \n', r, 1e0);
end
for r = 1:tot_edge
  fprintf(spring_fid, '%d %d %1.16e %1.16e\n', edge(1,r), edge(2,r), 0.0, length(r));
end

fclose(vertex_fid);
fclose(spring_fid);
fclose(target_fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
