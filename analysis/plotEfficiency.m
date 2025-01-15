fontsize = 38;
%set(groot, 'defaulttextinterpreter', 'tex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
%set(groot, 'defaultLegendInterpreter', 'tex');

alpha_sets = {["0.0", "0.001", "0.01", "0.05", "0.1", "0.3"], ["0.0", "0.001", "0.01", "0.05", "0.1", "0.3"], ["0.3"]};
De_sets = {["0.5"], ["1.0"], ["0.01", "0.1", "0.3", "0.5", "1.0", "2.0"]};

color_list = {[0 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]};

f = figure(1); clf; hold on;
tcl = tiledlayout(3,1);
for l = 1:length(alpha_sets)
    alpha_strs = alpha_sets{l};
    De_strs = De_sets{l};
    root_path = "/scratch/bistable_helix/data";
    paths = [root_path + "/De=0.0"];
    colors = {[0 0 0]};
    plot_styles = {"-."};
    legends = ["N1"];
    i = 1;
    for alpha_str = alpha_strs
        for De_str = De_strs
            paths = [paths, root_path + "/giesekus/alpha=" + alpha_str + "/De=" + De_str];
            legends = [legends, "$\alpha= " + alpha_str + ",\De=" + De_str + "$"];
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
    
    % Get efficiency of N1
    torque_N1 = getTorques(root_path + "/De=0.0/torque");
    torque_N1(4,:) = -torque_N1(4,:);
    torque_N1(4,:) = smoothdata(torque_N1(4,:),'gaussian', 5);
    flux_N1 = getAdvectiveFlux(root_path + "/De=0.0");
    flux_N1(2,:) = smoothdata(flux_N1(2,:),'gaussian', 5);
    ke_N1 = getKE(root_path + "/De=0.0");
    ke_N1(2,:) = smoothdata(ke_N1(2,:), 'gaussian', 5);
    % Old wrong definition
    %E_N1 = torque_N1(4,:).*flux_N1(2,:);
    % New definition
    %E_N1 = ke_N1(2,:) ./ torque_N1(4,:);
    E_N1 = flux_N1(2,:) ./ torque_N1(4,:);
    
    nexttile(l);
    i = 1;
    for path = paths
        base_path = path;
        fprintf("getting torque from " + base_path + "\n");
        torque = getTorques(base_path + "/torque");
        fprintf("Read " + num2str(length(torque(1,:))) + " time values\n");
        torque(4,:) = -torque(4,:);
        torque(4,:) = smoothdata(torque(4,:), 'gaussian', 5);
        
        flux = getAdvectiveFlux(base_path);
        flux(2,:) = smoothdata(flux(2,:), 'gaussian', 5);

        ke = getKE(base_path);
        ke(2,:) = smoothdata(ke(2,:), 'gaussian', 5);
    
        % Old definition
        %E_curr = torque(4,:).*flux(2,:);
        %E = E_N1./E_curr;
        % New definition
        %E_curr = ke(2,:) ./ torque(4,:);
        E_curr = flux(2,:) ./ torque(4,:);
        E = E_curr./E_N1;
        
        hold on;
        plot(torque(1,:), E, plot_styles{i}, 'Color', colors{i}, 'linewidth', 5, 'MarkerSize', 8);
        i = i+1;
    end
    if (l == 3)
        xlabel('Time (s)');
    else
        set(gca, 'XTick', []);
    end
    if (l == 2)
        ylabel('Efficiency');
    end
    %legend(legends);
    set(gca, 'fontsize', fontsize);
    axis([0.05 0.1 0.95 1.6]);
    if (l == 2)
        axis([0.05 0.1 0.95 2]);
    end
    %axis([0.04 0.095 0.75e-3 2.5e-3]); % De = 1.0, vary alpha
    %axis([0.04 0.095 0.8e-3*1000 1.2e-3*1000]); % De = 0.5, vary alpha
    %axis([0.04 0.095 0.8e-3 1.2e-3]); % vary De, alpha = 0.3
    h = legend(legends);
    h.Layout.Tile = 'East';

end