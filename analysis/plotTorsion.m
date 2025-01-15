clear;
%set(groot, 'defaulttextinterpreter', 'latex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
%set(groot, 'defaultLegendInterpreter', 'latex');
markersize = 300;
linewidth = 4;
plot_norm = 0;

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
%x_vals = [1.0 2.0 3.0 4.0];


i = 1;
for alpha_str = alpha_strs
    for De_str = De_strs
        paths = [paths, root_path + "/giesekus/alpha=" + alpha_str + "/De=" + De_str];
        
    end
    colors{end+1} = color_list{i};
    i = i+1;
end
%paths = [root_path + "/L/2.0", root_path + "/L/3.0", root_path + "/L/4.0"];
%colors = {color_list{1}, color_list{2}, color_list{3}};

prs = getPitchAndRadius(root_path + "/De=0.0/pitch", [50, 100]);
prs = abs(prs);
[N1_mid, ~, ~] = getAvgPitchAndRadius(prs(:,:,1));
[N1_end, ~, ~] = getAvgPitchAndRadius(prs(:,:,2));
N1_mid = N1_mid(2) / (N1_mid(1)^2 + N1_mid(2)^2);
N1_end = N1_end(2) / (N1_end(1)^2 + N1_end(2)^2);

prs = getPitchAndRadius(root_path + "/De=0.0_scaled/pitch", [50, 100]);
prs = abs(prs);
[N2_mid, ~, ~] = getAvgPitchAndRadius(prs(:,:,1));
[N2_end, ~, ~] = getAvgPitchAndRadius(prs(:,:,2));
N2_mid = N2_mid(2) / (N2_mid(1)^2 + N2_mid(2)^2);
N2_end = N2_end(2) / (N2_end(1)^2 + N2_end(2)^2);

figure(1); clf; tiledlayout(2,1);
f.Renderer = 'painters';
for i = 1:length(paths)
    base_path = paths(i) + "/pitch";
    fprintf("getting pitch from " + base_path + "\n");
    prs = getPitchAndRadius(base_path, [50, 100]);
    prs = abs(prs);
    [avg_mid, max_mid, min_mid] = getAvgPitchAndRadius(prs(:,:,1));
    [avg_end, max_end, min_end] = getAvgPitchAndRadius(prs(:,:,2));
    
    torsion_mid = avg_mid(2) / (avg_mid(1).^2 + avg_mid(2).^2);
    torsion_end = avg_mid(2) / (avg_mid(1).^2 + avg_mid(2).^2);
    
    nexttile(1); hold on;
    e = scatter(x_vals(rem(i-1, 8) + 1), torsion_mid, markersize, marker_shape(ceil(i / 8)));
    e.MarkerEdgeColor = colors{ceil(i / 8) };
    e.MarkerFaceColor = colors{ceil(i / 8) };
    e.LineWidth = 1.0;
    nexttile(2); hold on;
    e = scatter(x_vals(rem(i-1, 8) + 1), torsion_end, markersize, marker_shape(ceil(i / 8)));
    e.MarkerEdgeColor = colors{ceil(i / 8) };
    e.MarkerFaceColor = colors{ceil(i / 8) };
    e.LineWidth = 1.0;
end
nexttile(1);
h = gca;
set(h, 'fontsize', 24);
%xlabel('De');
ylabel('Radius');
title('Middle');
%set(h, 'xscale', 'log');
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
if (plot_norm)
    yline(N2_mid(1) / N1_mid(1), '-', 'N^{low}', 'linewidth', 6, 'fontsize', 24);
else
    yline(N1_mid(1), '-', 'N1', 'linewidth', 4, 'fontsize', 24);
    yline(N2_mid(1), '-', 'N2', 'linewidth', 4, 'fontsize', 24);
end

nexttile(2);
h = gca;
set(h, 'fontsize', 24);
%xlabel('De');
title('End');
%set(h, 'xscale', 'log');
h.XAxis.TickLength = [0 0];
h.YAxis.TickLength = [0 0];
if (plot_norm)
    yline(N2_end(1) / N1_end(1), '-', 'N^{low}', 'linewidth', 6, 'fontsize', 24);
else
    yline(N1_end(1), '-', 'N1', 'linewidth', 4, 'fontsize', 24);
    yline(N2_end(1), '-', 'N2', 'linewidth', 4, 'fontsize', 24);
end

nexttile(2);
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
