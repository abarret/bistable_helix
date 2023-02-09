set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

root_path = "/scratch/bistable_helix/data";
De_strs = ["De=0.0", "De=0.01", "De=0.1", "De=0.3", "De=0.5", "De=0.7", "De=1.0", "De=2.0"];
legend_strs = ["0.0", "0.01", "0.1", "0.3", "0.5", "0.7", "1.0", "2.0"];

figure(1); clf; hold on;
figure(2); clf; hold on;
for De_str = De_strs
    base_path = root_path + '/' + De_str + '/pitch';
    fprintf("getting pitch from " + base_path + "\n");
    prs = getPitchAndRadius(base_path, [50, 100]);
    prs(3,:) = smoothdata(prs(3,:), 'gaussian', 6);
    prs(2,:) = smoothdata(prs(2,:), 'gaussian', 6);
    
    figure(1);
    subplot(2,1,1); hold on;
    plot(prs(1,:,1), prs(3,:,1), 'linewidth', 4);
    subplot(2,1,2); hold on;
    plot(prs(1,:,1), prs(2,:,1), 'linewidth', 4);
    
    figure(2);
    subplot(2,1,1); hold on;
    plot(prs(1,:,2), prs(3,:,2), 'linewidth', 4);
    subplot(2,1,2); hold on;
    plot(prs(1,:,2), prs(2,:,2), 'linewidth', 4);
end
figure(1);
subplot(2,1,1);
ylabel('Pitch ($\mu$m)');
legend(legend_strs);
set(gca, 'fontsize', 24);
axis([0.04 0.08 -1.94 -1.88])
subplot(2,1,2);
xlabel('Time (s)');
ylabel('Radius ($\mu$m)');
%legend(legend_strs);
set(gca, 'fontsize', 24);
axis([0.04 0.08 0.155 0.168])

figure(2);
subplot(2,1,1);
ylabel('Pitch ($\mu$m)');
legend(legend_strs);
set(gca, 'fontsize', 24);
axis([0.04 0.08 -2.12 -2.11])
subplot(2,1,2);
xlabel('Time (s)');
ylabel('Radius ($\mu$m)');
%legend(legend_strs);
set(gca, 'fontsize', 24);
axis([0.04 0.08 0.2025 0.204])
