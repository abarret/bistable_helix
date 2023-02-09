set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

root_path = "/scratch/bistable_helix/data";
De_strs = ["De=0.0", "De=0.01", "De=0.1", "De=0.3", "De=0.5", "De=0.7", "De=1.0", "De=2.0"];
legend_strs = ["0.0", "0.01", "0.1", "0.3", "0.5", "0.7", "1.0", "2.0"];
ang_vel = 2.0*pi*100;

figure(1); clf; hold on;
figure(2); clf; hold on;
for De_str = De_strs
    base_path = root_path + '/' + De_str + '/torque';
    fprintf("getting torque from " + base_path + "\n");
    torque = getTorques(base_path);
    fprintf("Read " + num2str(length(torque(1,:))) + " time values\n");
    torque(4,:) = smoothdata(torque(4,:), 'gaussian', 10);
    figure(1); hold on;
    plot(torque(1,:),torque(4,:),'linewidth', 4);
    
    figure(2); hold on;
    power = torque(4,:) * ang_vel;
    plot(torque(1,:), power, 'linewidth', 4);
end

figure(1);
xlabel('Time (s)');
ylabel('Torque');
legend(legend_strs);
set(gca, 'fontsize', 18);
axis([0.04 0.095 2.5e-5 4.5e-5]);

figure(2);
xlabel('Time (s)');
ylabel('Power');
legend(legend_strs);
set(gca, 'fontsize', 18);
axis([0.04 0.095 0.018 0.028]);