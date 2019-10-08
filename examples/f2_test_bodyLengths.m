function f2_test_bodyLengths

lengths = [4:10];

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


%% show bodies:
figure('position', [ 94         421        1627         417]);
ylabel('length (mm)')
hold on;

c_temp = [0 0 0];
for i = 1:length(lengths)   
    opts.worm.length = lengths(i); %mm
    opts.worm.initialPosition = c_temp;

    w = gen_worm(opts);
    w.singleShow;
    
    c_temp = c_temp + [2 0 0];
end



%% make simulations
for i = 1:length(lengths)    
  
    opts.worm.length = lengths(i); %mm
    w = gen_worm(opts);
    
    [traj_i, angles_i] = runStuff(w, opts);
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
ylim([0 180])
title('angle')

cols = jet(numel(lengths));

for i = 1:length(lengths)
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
legend(strread(num2str(lengths), '%s'),'location','southeast')
subplot(1,2,2)
legend(strread(num2str(lengths), '%s'))



if 0 
% close all;
tr_fig = figure;
% subplot(1,2,1)
title('trajectories')
axis equal
hold on;

ang_fig = figure;
% subplot(1,2,2)
hold on
ylim([45 135])
title('angle')

cols = jet(numel(lengths));

for i = 1:length(lengths)
    traj_i = allTrajectories{i};
    angles_i = allAngles{i};
    figure(tr_fig)
    plot(traj_i(:, 1), traj_i(:, 2), '-', 'color', cols(i,:), 'LineWidth', 2)
    %text(traj_i(end, 1), traj_i(end,2 ), num2str(bodyRatios(i)))
    figure(ang_fig)
    plot(angles_i, '-', 'color', cols(i,:), 'LineWidth', 2)
    %text(length(angles_i), angles_i(end), num2str(bodyRatios(i)))
end

figure(tr_fig)
legend(strread(num2str(lengths), '%s'))
figure(ang_fig)
legend(strread(num2str(lengths), '%s'))
end

