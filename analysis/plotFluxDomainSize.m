clear;
set(groot, 'defaulttextinterpreter', 'tex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
%set(groot, 'defaultLegendInterpreter', 'latex');
markersize = 150;
linewidth = 2;

color_list = {[0 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
%color_list = {[0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410]};

root_path = "/scratch/bistable_helix/data";
paths = [];
colors = {};

%De_strs = ["0.5"];
%alpha_strs = ["0.001", "0.01", "0.05", "0.1", "0.3"];
%x_vals = [0.001, 0.01, 0.05, 0.1, 0.3];
alpha_strs = ["0.01", "0.05", "0.1", "0.2", "0.3"];
De_strs = ["0.01", "0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"];
x_vals = [20.0 40.0 60.0 80.0];

paths = [root_path + "/De=0.0", root_path + "/L/2.0", root_path + "/L/3.0", root_path + "/L/4.0"];
colors = {color_list{1}, color_list{2}, color_list{3}};

flux = getAdvectiveFlux(root_path + "/De=0.0");
avg_flux_1 = mean(flux(2, flux(1,:) > 0.04));

figure(1); clf; hold on;
f.Renderer = 'painters';
figure(2); clf; hold on;
for i = 1:length(paths)
    base_path = paths(i);
    fprintf("getting pitch from " + base_path + "\n");
    flux = getAdvectiveFlux(base_path);
    
    t = flux(1,:);
    t_idxs = t > 0.04;
    avg_flux = mean(flux(2,t_idxs));
    max_flux = max(flux(2,t_idxs));
    min_flux = min(flux(2,t_idxs));

    figure(2); hold on;
    e = scatter(x_vals(rem(i, 8)) / 6.0, avg_flux / avg_flux_1, markersize);
    e.MarkerEdgeColor = colors{ceil(i / 8) };
    e.MarkerFaceColor = colors{ceil(i / 8) };
    e.LineWidth = 1.0;
    e = errorbar(x_vals(rem(i, 8)) / 6.0, avg_flux(1) / avg_flux_1, (avg_flux(1) - min_flux(1)) / avg_flux_1, (max_flux(1) - avg_flux(1)) / avg_flux_1);
    e.Color = colors{ceil(i / 8) };
    e.LineWidth = linewidth;
    e.CapSize = 10;
end
h = gca;
set(h, 'fontsize', 24);
xlabel('Domain size / L');
ylabel('Normalized Flow Rate');
%set(h, 'xscale', 'log');
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
