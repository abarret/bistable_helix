f = @(f,s,t) (0.25 + s).*sin(2*pi*f*(s - t));

for t = 0:0.02:10
    figure(1); clf; hold on;
    svals = linspace(0,3,100);
    plot(svals, f(1, svals, t), 'linewidth', 4);
    plot(svals, f(0.8, svals, t), 'linewidth', 4);
    plot(svals, f(1.2, svals, t), 'linewidth', 4);
    legend('f = 1', 'f = 0.8', 'f = 1.2');
    axis([0, 3, -4, 4]);
end