De_strs = ["0.0", "0.01", "0.05", "0.1", "0.3"];
%De_strs = ["0.3"];

figure(1); clf; hold on;
figure(2); clf; hold on;
for De_str = De_strs
    base_path = 'De_' + De_str + '/pitch';
    fprintf("getting pitch from " + base_path + "\n");
    prs = getPitchAndRadius(base_path, [50, 100]);
    prs(3,:) = smoothdata(prs(3,:), 'rloess', 8);
    prs(2,:) = smoothdata(prs(2,:), 'rloess', 8);
    
    figure(1); hold on;
    plot(prs(1,:,1), prs(3,:,1), 'linewidth', 4);
    %plot(prs(1,:,2), prs(3,:,2), 'linewidth', 4);
    

    figure(2); hold on;
    plot(prs(1,:,1), prs(2,:,1), 'linewidth', 4);
    %plot(prs(1,:,2), prs(2,:,2), 'linewidth', 4);
end
figure(1);
xlabel('time');
ylabel('Pitch');
legend(De_strs);
set(gca, 'fontsize', 18);
axis([0.02 0.08 -1.9 -1.8])

figure(2); 
xlabel('time');
ylabel('radius');
legend(De_strs);
set(gca, 'fontsize', 18);
axis([0.02 0.08 0.144 0.158])