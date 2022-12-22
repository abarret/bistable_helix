base_path = 'torque';

torque = getTorques(base_path);

figure(1); clf;
plot(torque(1,:),torque(4,:),'linewidth', 4);