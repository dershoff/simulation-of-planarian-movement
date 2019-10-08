function  f2_test_stroke_force_effect

strokes     = [0.2    0.4   0.6   0.8   1.0];
frames      = [10000  5000  3330  2500  2000];


opts.worm.length = 8; %mm
opts.worm.bodyRatio = 3;

opts.simulation.frames = 5000;
opts.simulation.showFrames = 0;
opts.simulation.showStamps = 0;
opts.simulation.saveVideo = 0;
opts.simulation.VideoName = 'temp_test';

opts.worm.mutant = 'wt'; 
% 'vfl3' - left
% 'odf2' - right
% 'wt'
opts.worm.initialPosition = [0, 0, 0];
opts.worm.orientation = 90; %90 = upwards
opts.worm.strokeForce = 0.2; %arbitrary unit
opts.worm.mass = 10; %arbitrary
opts.worm.dt = 1; %time step for integration    




allTrajectories = {};
allAngles = {};


%% make simulations
for i = 1:length(strokes)    
  
    opts.worm.strokeForce = strokes(i); %mm
    opts.simulation.frames = frames(i);
    
    w = gen_worm(opts);
    
    [traj_i, angles_i] = runStuff(w, opts);
    allTrajectories = [allTrajectories, {traj_i}];
    allAngles = [allAngles, {angles_i}];
end


% sort trajes by length
all_y = [];
for i = 1:length(strokes)
    traj = allTrajectories{i};
    all_y = [all_y, allTrajectories{i}(end,2)];
end
[val, ind] = sort(all_y);
thick_set = [length(strokes) : -1 : 1];
allTrajectories = allTrajectories(ind);
strokes = strokes(ind);






% close all;
figure;
subplot(1,2,1)
title('trajectories')
axis equal
hold on;
subplot(1,2,2)
hold on
ylim([0 180])
title('angle')

cols = jet(length(strokes));

for i = 1:length(strokes)
    traj_i = allTrajectories{i};
    angles_i = allAngles{i};
    subplot(1,2,1)
    plot(traj_i(:, 1), traj_i(:, 2), '-', 'color', cols(i,:), 'lineWidth', thick_set(i))
    %text(traj_i(end, 1), traj_i(end,2 ), num2str(bodyRatios(i)))
    subplot(1,2,2)
    plot(angles_i, '-', 'color', cols(i,:), 'LineWidth', 2)
    %text(length(angles_i), angles_i(end), num2str(bodyRatios(i)))
end

subplot(1,2,1)
legend(strread(num2str(strokes),'%s'),'location','southeast')
for i = 1:length(strokes)
    traj = allTrajectories{i};
    plot(traj(end,1), traj(end,2), '^', 'color', cols(i,:), 'lineWidth', thick_set(i))
end

