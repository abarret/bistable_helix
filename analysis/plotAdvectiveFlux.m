root_path = "/scratch/bistable_helix/data";
De_strs = ["0.0", "0.01", "0.1", "0.3", "1.0", "2.0"];
legend_strs = ["0.0", "0.01", "0.1", "0.3", "1.0", "2.0"];

figure(1); clf; hold on;
for De_str = De_strs
    base_path = root_path + '/De=' + De_str;
    fprintf("getting flux from " + base_path + "\n");
    flux = getAdvectiveFlux(base_path);
    flux(2,:) = smoothdata(flux(2,:), 'gaussian', 6);
    
    figure(1); hold on;
    plot(flux(1,:), flux(2,:), 'linewidth', 4);
end
figure(1);
xlabel('time');
ylabel('Volumetric Flow Rate');
legend(legend_strs);
set(gca, 'fontsize', 18);
axis([0.04 0.1 18 24])