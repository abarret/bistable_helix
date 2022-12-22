base_path = 'pitch';

prs = getPitchAndRadius(base_path, [50, 100]);

figure(1); clf; hold on;
plot(prs(1,:,1), prs(3,:,1), 'linewidth', 4);
plot(prs(1,:,2), prs(3,:,2), 'linewidth', 4);
xlabel('time');
ylabel('Pitch');
legend('half', 'end');
set(gca, 'fontsize', 18);

figure(2); clf; hold on;
plot(prs(1,:,1), prs(2,:,1), 'linewidth', 4);
plot(prs(1,:,2), prs(2,:,2), 'linewidth', 4);
xlabel('time');
ylabel('radius');
legend('half', 'end');
set(gca, 'fontsize', 18);