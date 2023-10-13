clear;
fontsize = 38;
%set(groot, 'defaulttextinterpreter', 'tex');
%set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
%set(groot, 'defaultLegendInterpreter', 'tex');

alpha_sets = {["0.01"], ["0.05"], ["0.1"], ["0.2"], ["0.3"]};
De_sets = {["0.01", "0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"], ["0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"], ["0.01", "0.1", "0.3", "0.5", "0.8", "1.0", "1.5", "2.0"], ["0.01", "0.1", "0.3", "0.5", "0.8", "1.5", "2.0"], ["0.01", "0.1", "0.3", "0.5", "1.0", "2.0"]};
%De_sets = {["0.5"], ["1.0"], ["0.01"], ["0.01", "0.1", "0.3", "0.5", "1.0", "2.0"]};
root_path = "/scratch/bistable_helix/data";

% Get efficiency of N1
torque_N1 = getTorques(root_path + "/De=0.0/torque");
torque_N1(4,:) = -torque_N1(4,:);
torque_N1(4,:) = smoothdata(torque_N1(4,:),'gaussian', 5);
flux_N1 = getAdvectiveFlux(root_path + "/De=0.0");
flux_N1(2,:) = smoothdata(flux_N1(2,:),'gaussian', 5);
E_N1 = flux_N1(2,:) ./ torque_N1(4,:);
E_N1 = mean(E_N1(100:end));


torque = getTorques(root_path + "/De=0.0_scaled/torque");
torque(4,:) = -torque(4,:);
torque(4,:) = smoothdata(torque(4,:),'gaussian', 5);
flux = getAdvectiveFlux(root_path + "/De=0.0_scaled");
flux(2,:) = smoothdata(flux(2,:),'gaussian', 5);
E = flux(2,:) ./ torque(4,:);
EN2 = E(100:end)./E_N1;
EN2 = mean(EN2);
%alpha_vals = [0.0];
%De_vals = [0.0];
%E_vals = [mean(E)];
alpha_vals = [];
De_vals = [];
E_vals = [];

for l = 1:length(alpha_sets)
    alpha_strs = alpha_sets{l};
    De_strs = De_sets{l};
    
    for alpha_str = alpha_strs
        for De_str = De_strs
            base_path = root_path + "/giesekus/alpha=" + alpha_str + "/De=" + De_str;
            fprintf("getting torque from " + base_path + "\n");
            torque = getTorques(base_path + "/torque");
            fprintf("Read " + num2str(length(torque(1,:))) + " time values\n");
            torque(4,:) = -torque(4,:);
            torque(4,:) = smoothdata(torque(4,:), 'gaussian', 5);
        
            flux = getAdvectiveFlux(base_path);
            flux(2,:) = smoothdata(flux(2,:), 'gaussian', 5);
            E_curr = flux(2,:) ./ torque(4,:);
            E = E_curr(100:end)./E_N1;
            
            alpha_vals = [alpha_vals, str2num(alpha_str)];
            De_vals = [De_vals, str2num(De_str)];
            E_vals = [E_vals, mean(E)];
        end
    end
end

figure(1); clf;
scatter(alpha_vals, De_vals, 200, E_vals, 'MarkerFaceColor', 'flat');
colorbar;
colormap jet
xlabel('alpha');
ylabel('De');
set(gca, 'fontsize', 36);

F = scatteredInterpolant(alpha_vals', De_vals', E_vals');
plot_alphas = linspace(0.01, 0.3, 100);
plot_Des = linspace(0.01, 2.0, 100);
[aa, dd] = meshgrid(plot_alphas, plot_Des);
FF = F(aa, dd);
figure(2); clf;
h = pcolor(aa, dd, FF);
set(h, 'Edgecolor', 'none');
hold on;
contour(aa, dd, FF, [EN2, EN2], 'k', 'linewidth', 8);
scatter(alpha_vals, De_vals, 200, 'black', 'MarkerFaceColor', 'flat');
xlabel('alpha'); ylabel('De');
colorbar();
set(gca, 'fontsize', 32);
