%set(groot, 'defaulttextinterpreter', 'latex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
%set(groot, 'defaultLegendInterpreter', 'latex');

color_list = {[0 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};

root_path = "/scratch/bistable_helix/data";
paths = [root_path + "/De=0.0"];
legends = ["N1"];
colors = {[0 0 0]};
plot_styles = {"-."};

De_strs = ["0.5"];
alpha_strs = ["0.0", "0.001", "0.01", "0.05", "0.1", "0.3"];

i = 1;
for alpha_str = alpha_strs
    for De_str = De_strs
        paths = [paths, root_path + "/giesekus/alpha=" + alpha_str + "/De=" + De_str];
        legends = [legends, "$\\alpha= " + alpha_str + ",\\De=" + De_str + "$"];
        colors{end+1} = color_list{i};
        plot_styles{end+1} = "-";
        i = i+1;
    end
end

paths = [paths, root_path + "/De=0.0_scaled"];
legends = [legends, "N2"];
colors{end+1} = [0 0 0];
plot_styles{end+1} = ":";

figure(1); clf; tiledlayout(2,2);
f.Renderer = 'painters';
for i = 1:length(paths)
    base_path = paths(i) + "/pitch";
    fprintf("getting pitch from " + base_path + "\n");
    prs = getPitchAndRadius(base_path, [50, 100]);
    prs(3,:) = smoothdata(prs(3,:), 'gaussian', 6);
    prs(2,:) = smoothdata(prs(2,:), 'gaussian', 6);
    
    nexttile(1); hold on;
    plot(prs(1,:,1), abs(prs(3,:,1)), plot_styles{i}, 'linewidth', 4, 'Color', colors{i});
    nexttile(3); hold on;
    plot(prs(1,:,1), prs(2,:,1), plot_styles{i}, 'linewidth', 4, 'Color', colors{i});
    
    nexttile(2); hold on;
    plot(prs(1,:,2), abs(prs(3,:,2)), plot_styles{i}, 'linewidth', 4, 'Color', colors{i});
    nexttile(4); hold on;
    plot(prs(1,:,2), prs(2,:,2), plot_styles{i}, 'linewidth', 4, 'Color', colors{i});
end
nexttile(1);
ylabel('Pitch ($\\mu$m)');
l = legend(legends);
set(gca, 'fontsize', 24);
axis([0.04 0.08 1.88 1.96])
ax = gca;
xax = ax.XAxis;
xax.TickValues = [0.04 0.05 0.06 0.07 0.08];
yax = ax.YAxis;
yax.TickValues = [1.89 1.92 1.95];

nexttile(3);
xlabel('Time (s)');
ylabel('Radius ($\\mu$m)');
%legend(legend_strs);
set(gca, 'fontsize', 24);
axis([0.04 0.08 0.156 0.1688])
ax = gca;
xax = ax.XAxis;
xax.TickValues = [0.04 0.05 0.06 0.07 0.08];
yax = ax.YAxis;
yax.TickValues = [0.156 0.16 0.164 0.168];

nexttile(2);
%ylabel('Pitch ($\\mu$m)');
set(gca, 'fontsize', 24);
axis([0.04 0.08 2.106 2.12])
ax = gca;
xax = ax.XAxis;
xax.TickValues = [0.04 0.05 0.06 0.07 0.08];
yax = ax.YAxis;
yax.TickValues = [2.108 2.113 2.118];

nexttile(4);
xlabel('Time (s)');
%ylabel('Radius ($\\mu$m)');
%legend(legend_strs);
set(gca, 'fontsize', 24);
axis([0.04 0.08 0.2015 0.2042])
ax = gca;
xax = ax.XAxis;
xax.TickValues = [0.04 0.05 0.06 0.07 0.08];
yax = ax.YAxis;
yax.TickValues = [0.202 0.203 0.204];

l.Layout.Tile = 'east';