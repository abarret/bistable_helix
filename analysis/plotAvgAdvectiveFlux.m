clear;
%set(groot, 'defaulttextinterpreter', 'latex');
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

t = flux(1,:);
t_idxs = t > 0.04;
N1_flux = mean(flux(2,t_idxs));

flux = getAdvectiveFlux(root_path + "/De=0.0_scaled");

t = flux(1,:);
t_idxs = t > 0.04;
N2_flux = mean(flux(2,t_idxs));

figure(1); clf; hold on;
f.Renderer = 'painters';
for i = 1:length(paths)
    base_path = paths(i);
    fprintf("getting pitch from " + base_path + "\n");
    flux = getAdvectiveFlux(base_path);
    
    t = flux(1,:);
    t_idxs = t > 0.04;
    avg_flux = mean(flux(2,t_idxs));
    max_flux = max(flux(2,t_idxs));
    min_flux = min(flux(2,t_idxs));
    
    nexttile(1); hold on;
    e = scatter(x_vals(rem(i, 8)+1), avg_flux, markersize);
    e.MarkerEdgeColor = colors{ceil(i / 8) };
    e.MarkerFaceColor = colors{ceil(i / 8) };
    e.LineWidth = 1.0;
    e = errorbar(x_vals(rem(i, 8)+1), avg_flux(1), (avg_flux(1) - min_flux(1)), (max_flux(1) - avg_flux(1)));
    e.Color = colors{ceil(i / 8) };
    e.LineWidth = linewidth;
    e.CapSize = 10;
end
h = gca;
set(h, 'fontsize', 24);
%xlabel('De');
ylabel('Radius (\mum)');
%set(h, 'xscale', 'log');
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
yline(N1_flux(1), '-', 'N1', 'linewidth', 4, 'fontsize', 24);
yline(N2_flux(1), '-', 'N2', 'linewidth', 4, 'fontsize', 24);

c1 = plot([NaN NaN], [NaN NaN], 'Color', colors{1}, 'linewidth', linewidth);
c2 = plot([NaN NaN], [NaN NaN], 'Color', colors{2}, 'linewidth', linewidth);
c3 = plot([NaN NaN], [NaN NaN], 'Color', colors{3}, 'linewidth', linewidth);
c4 = plot([NaN NaN], [NaN NaN], 'Color', colors{4}, 'linewidth', linewidth);
c5 = plot([NaN NaN], [NaN NaN], 'Color', colors{5}, 'linewidth', linewidth);
l = legend([c1, c2, c3, c4, c5], 'alpha = 0.01', 'alpha = 0.05', 'alpha = 0.1', 'alpha = 0.2', 'alpha = 0.3');
l.Layout.Tile = 'east';