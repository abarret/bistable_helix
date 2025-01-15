clear;
%set(groot, 'defaulttextinterpreter', 'latex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
%set(groot, 'defaultLegendInterpreter', 'latex');
markersize = 300;
linewidth = 4;
plot_norm = 1;

color_list = {[0 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};
%color_list = {[0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410], [0 0.4470 0.7410]};
marker_shape = ["d", "o", "s", "^", "<", ">"];

root_path = "/scratch/bistable_helix/data";
paths = [];
colors = {};

%De_strs = ["0.5"];
%alpha_strs = ["0.001", "0.01", "0.05", "0.1", "0.3"];
%x_vals = [0.001, 0.01, 0.05, 0.1, 0.3];
alpha_strs = ["0.01", "0.05", "0.1", "0.2", "0.3"];
De_strs = ["0.01", "0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"];
x_vals = [0.01, 0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0];


i = 1;
for alpha_str = alpha_strs
    for De_str = De_strs
        paths = [paths, root_path + "/giesekus/alpha=" + alpha_str + "/De=" + De_str];
        
    end
    colors{end+1} = color_list{i};
    i = i+1;
end

flux = getAdvectiveFlux(root_path + "/De=0.0");
torque = getTorques(root_path + "/De=0.0/torque");

t = flux(1,:);
t_idxs = t > 0.04;
Nhigh_flux = mean(flux(2,t_idxs));
Nhigh_torque = abs(mean(torque(4,torque(1,:) > 0.06)));

flux = getAdvectiveFlux(root_path + "/De=0.0_scaled");
torque = getTorques(root_path + "/De=0.0_scaled/torque");

t = flux(1,:);
t_idxs = t > 0.04;
Nlow_flux = mean(flux(2,t_idxs));
Nlow_torque = abs(mean(torque(4,torque(1,:) > 0.06)));

figure(1); clf; hold on;
tiledlayout(2,1);
f.Renderer = 'painters';
for i = 1:length(paths)
    base_path = paths(i);
    fprintf("getting pitch from " + base_path + "\n");
    flux = getAdvectiveFlux(base_path);
    
    t = flux(1,:);
    t_idxs = t > 0.06;
    avg_flux = mean(flux(2,t_idxs));
    max_flux = max(flux(2,t_idxs));
    min_flux = min(flux(2,t_idxs));
    
    nexttile(1); hold on;
    if (plot_norm)
        e = scatter(x_vals(rem(i-1, 8) + 1), avg_flux / Nhigh_flux, markersize, marker_shape(ceil(i / 8)));
        e.MarkerEdgeColor = colors{ceil(i / 8) };
        e.MarkerFaceColor = colors{ceil(i / 8) };
        e.LineWidth = 1.0;
        e = errorbar(x_vals(rem(i-1, 8) + 1), avg_flux(1) / Nhigh_flux, (avg_flux(1) - min_flux(1)) / Nhigh_flux, (max_flux(1) - avg_flux(1)) / Nhigh_flux);
        e.Color = colors{ceil(i / 8) };
        e.LineWidth = linewidth;
        e.CapSize = 10;
    else
        e = scatter(x_vals(rem(i-1, 8) + 1), avg_flux, markersize, marker_shape(ceil(i / 8)));
        e.MarkerEdgeColor = colors{ceil(i / 8) };
        e.MarkerFaceColor = colors{ceil(i / 8) };
        e.LineWidth = 1.0;
        e = errorbar(x_vals(rem(i-1, 8) + 1), avg_flux(1), (avg_flux(1) - min_flux(1)), (max_flux(1) - avg_flux(1)));
        e.Color = colors{ceil(i / 8) };
        e.LineWidth = linewidth;
        e.CapSize = 10;
    end
    
    % Torque
    base_path = paths(i) + '/torque';
    fprintf("getting torque from " + base_path + "\n");
    torque = getTorques(base_path);
    
    t = torque(1,:);
    t_idxs = t > 0.06;
    avg_torque = abs(mean(torque(4,torque(1,:) > 0.04)));
    
    nexttile(2); hold on;
    if (plot_norm)
        e = scatter(x_vals(rem(i-1, 8) + 1), avg_torque / Nhigh_torque, markersize, marker_shape(ceil(i / 8)));
        e.MarkerEdgeColor = colors{ceil(i / 8) };
        e.MarkerFaceColor = colors{ceil(i / 8) };
        e.LineWidth = 1.0;
    else
        e = scatter(x_vals(rem(i-1, 8) + 1), avg_torque, markersize, marker_shape(ceil(i / 8)));
        e.MarkerEdgeColor = colors{ceil(i / 8) };
        e.MarkerFaceColor = colors{ceil(i / 8) };
        e.LineWidth = 1.0;
    end
end
nexttile(1); hold on;
h = gca;
set(h, 'fontsize', 24);
%xlabel('De');
ylabel('Flow rate / N^{high}');
%set(h, 'xscale', 'log');
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
if (plot_norm)
    yline(Nlow_flux(1) / Nhigh_flux, '-', 'N^{low}', 'linewidth', 6, 'fontsize', 24);
else
    yline(Nhigh_flux(1), '-', 'N^{high}', 'linewidth', 6, 'fontsize', 24);
    yline(Nlow_flux(1), '-', 'N^{low}', 'linewidth', 6, 'fontsize', 24);
end

nexttile(2); hold on;
h = gca;
set(h, 'fontsize', 24);
xlabel('De');
ylabel('Torque / N^{high}');
%set(h, 'xscale', 'log');
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
if (plot_norm)
    yline(Nlow_torque(1) / Nhigh_torque, '-', 'N^{low}', 'linewidth', 6, 'fontsize', 24);
else
    yline(Nhigh_torque(1), '-', 'N^{high}', 'linewidth', 6, 'fontsize', 24);
    yline(Nlow_torque(1), '-', 'N^{low}', 'linewidth', 6, 'fontsize', 24);
end

c1 = plot([NaN NaN], [NaN NaN], marker_shape(1), 'Color', colors{1}, 'linewidth', linewidth);
c1.MarkerSize = 15;
c1.MarkerFaceColor = colors{1};
c2 = plot([NaN NaN], [NaN NaN], marker_shape(2), 'Color', colors{2}, 'linewidth', linewidth);
c2.MarkerSize = 15;
c2.MarkerFaceColor = colors{2};
c3 = plot([NaN NaN], [NaN NaN], marker_shape(3), 'Color', colors{3}, 'linewidth', linewidth);
c3.MarkerSize = 15;
c3.MarkerFaceColor = colors{3};
c4 = plot([NaN NaN], [NaN NaN], marker_shape(4), 'Color', colors{4}, 'linewidth', linewidth);
c4.MarkerSize = 15;
c4.MarkerFaceColor = colors{4};
c5 = plot([NaN NaN], [NaN NaN], marker_shape(5), 'Color', colors{5}, 'linewidth', linewidth);
c5.MarkerSize = 15;
c5.MarkerFaceColor = colors{5};
l = legend([c1, c2, c3, c4, c5], '\alpha = 0.01', '\alpha = 0.05', '\alpha = 0.1', '\alpha = 0.2', '\alpha = 0.3');
l.Layout.Tile = 'east';