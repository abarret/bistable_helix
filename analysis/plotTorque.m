De_strs = ["0.0", "0.01", "0.05", "0.1", "0.3"];
%De_strs = ["0.3"];
ang_vel = 2.0*pi*100;

figure(1); clf; hold on;
figure(2); clf; hold on;
for De_str = De_strs
    base_path = 'De_' + De_str + '/torque';
    fprintf("getting torque from " + base_path + "\n");
    torque = getTorques(base_path);
    fprintf("Read " + num2str(length(torque(1,:))) + " time values\n");
    torque(4,:) = smoothdata(torque(4,:), 'rloess', 12);
    figure(1); hold on;
    plot(torque(1,:),torque(4,:),'linewidth', 4);
    
    figure(2); hold on;
    power = torque(4,:) * ang_vel;
    plot(torque(1,:), power, 'linewidth', 4);
end

figure(1);
xlabel('Time');
ylabel('Torque');
legend(De_strs);
set(gca, 'fontsize', 18);
axis([0.02 0.08 3e-5 6e-5]);

figure(2);
xlabel('Time');
ylabel('Power');
legend(De_strs);
set(gca, 'fontsize', 18);
axis([0.02 0.08 0.024 0.036]);