function f2_test_angles

angles = linspace(0, rad2deg(2*pi), 4);

opts.simulation.frames = 50;
opts.simulation.showFrames = 1;
opts.simulation.showStamps = 0;
opts.simulation.saveVideo = 0;
opts.simulation.VideoName = 'temp_test';

opts.worm.mutant = 'odf2'; 
% 'vfl3' - left
% 'odf2' - right
% 'wt'
opts.worm.bodyRatio = 3;
opts.worm.length = 8; %mm
opts.worm.initialPosition = [0, 0, 0];
opts.worm.strokeForce = 0.2; %arbitrary unit
opts.worm.mass = 10; %arbitrary
opts.worm.dt = 1; %time step for integration    




%% show bodies:
% figure('position', [ 94         421        1627         417]);
% ylabel('length (mm)')
% hold on;


allTrajectories = {};
allAngles = {};
plotOn = 0;

for i = 1:length(angles) 
    
    opts.worm.orientation = angles(i); %90 = upwards
    w = gen_worm(opts);    
       
    [traj_i, angles_i] = runStuff(w, opts, plotOn);
    
    allTrajectories = [allTrajectories, {traj_i}];
    allAngles = [allAngles, {angles_i}];
end




%% show results

% close all;
figure;
subplot(1,2,1)
title('trajectories')
axis equal
hold on;
subplot(1,2,2)
hold on
ylim([0 370])
title('angle')

cols = jet(numel(angles));

for i = 1:length(angles)
    traj_i = allTrajectories{i};
    angles_i = allAngles{i};
subplot(1,2,1)
    plot(traj_i(:, 1), traj_i(:, 2), '-', 'color', cols(i,:), 'LineWidth', 2)
    %text(traj_i(end, 1), traj_i(end,2 ), num2str(bodyRatios(i)))
subplot(1,2,2)
    plot(angles_i, '-', 'color', cols(i,:), 'LineWidth', 2)
    %text(length(angles_i), angles_i(end), num2str(bodyRatios(i)))
end

subplot(1,2,1)
legend(strread(num2str(angles), '%s'),'location','southeast')
subplot(1,2,2)
legend(strread(num2str(angles), '%s'))



