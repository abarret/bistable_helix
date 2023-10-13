fontsize = 38;
%set(groot, 'defaulttextinterpreter', 'tex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
%set(groot, 'defaultLegendInterpreter', 'tex');

color_list = {[0 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};

root_path = "/scratch/bistable_helix/data";
paths = [root_path + "/De=0.0"];
legends = ["N1"];
colors = {[0 0 0]};
plot_styles = {"-."};

%alpha_strs = ["0.0", "0.001", "0.01", "0.05", "0.1", "0.3"];
%De_strs = ["0.5"];

De_strs = ["0.01", "0.1", "0.3", "0.5", "1.0", "2.0"];
alpha_strs = ["0.3"];

i = 1;
for alpha_str = alpha_strs
    for De_str = De_strs
        paths = [paths, root_path + "/giesekus/alpha=" + alpha_str + "/De=" + De_str];
        legends = [legends, "$alpha= " + alpha_str + ",De=" + De_str + "$"];
        colors{end+1} = color_list{i};
        plot_styles{end+1} = "-";
        %legends = [legends, "$\alpha = $ " + alpha_str];
        %legends = [legends, "De = " + De_str];
        i = i+1;
    end
end

paths = [paths, root_path + "/De=0.0_scaled"];
legends = [legends, "N2"];
colors{end+1} = [0 0 0];
plot_styles{end+1} = ":";

ang_vel = 2.0*pi*100;

f = figure(1); clf; hold on;
f.Renderer = 'painters';
tcl = tiledlayout(1,2);
i = 1;
for path = paths
    base_path = path;
    fprintf("getting torque from " + base_path + "\n");
    torque = getTorques(base_path + "/torque");
    fprintf("Read " + num2str(length(torque(1,:))) + " time values\n");
    torque(4,:) = -torque(4,:)*1000;
    torque(4,:) = smoothdata(torque(4,:), 'gaussian', 5);
    
    flux = getAdvectiveFlux(base_path);
    flux(2,:) = smoothdata(flux(2,:), 'gaussian', 5);
    
    figure(1); hold on;
    nexttile(1); hold on;
    plot(torque(1,:),torque(4,:), plot_styles{i}, 'Color', colors{i}, 'linewidth', 6, 'MarkerSize', 8);
    nexttile(2); hold on;
    plot(flux(1,:), flux(2,:), plot_styles{i}, 'Color', colors{i}, 'linewidth', 6, 'MarkerSize', 8);
    i = i+1;
end

figure(1);
nexttile(1);
xlabel('Time (s)');
ylabel('Torque (units)');
%legend(legends);
set(gca, 'fontsize', fontsize);
%axis([0.04 0.095 0.75e-3 2.5e-3]); % De = 1.0, vary alpha
%axis([0.04 0.095 0.8e-3*1000 1.2e-3*1000]); % De = 0.5, vary alpha
axis([0.04 0.095 0.8e-3*1000 1.2e-3*1000]); % vary De, alpha = 0.3

nexttile(2);
xlabel('Time (s)');
ylabel('Volumetric Flow Rate (units)');
%legend(legends);
set(gca, 'fontsize', fontsize);
%axis([0.04 0.1 18 55]) % De=1.0, alpha varies
%axis([0.04 0.1 18 28]) % De=0.5, alpha varies
axis([0.04 0.1 18 24]) % alpha=0.3, De varies

h = legend(legends);
h.Layout.Tile = 'East';