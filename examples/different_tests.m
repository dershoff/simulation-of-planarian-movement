function different_tests

strokes = [0.01 0.05 0.1 0.5];
strokes = [0.01 0.5];
tIter = [30000, 800];

trajCollection = {};

for i = 1:2
s = strokes(i);
t = tIter(i);
    
%% INITIALIZE A WORM:
opts.simulation.frames = t;
opts.simulation.showFrames = 0;
opts.simulation.saveVideo = 0;

opts.worm.length = 8; %mm
opts.worm.bodyRatio = 2.8; %length to width, normal value = 2.8
opts.worm.initialPosition = [0, 0, 0];
opts.worm.orientation = 90; %90 = upwards

opts.worm.strokeForce = s; %arbitrary unit
opts.worm.mass = 10; %arbitrary
opts.worm.dt = 1; %time step for integration

opts.worm.mutant = 'vfl3'; %'vfl3'; % 'odf2'; %'wt'
% vfl - left
% odf - right

w = gen_worm(opts);
tr = runStuff(w, opts);


trajCollection = [trajCollection, {tr(:,1:2)}];

end

save('trajes_vs_strokes_2_vfl3.mat','trajCollection','strokes')

if 1
close all
figure; 
axis equal
hold on
cols = jet(length(strokes));
title([w.mutant '; trajes, parameter - stroke force'])
for i = 2:-1:1
tr = trajCollection{i};
    plot(tr(:,1), tr(:,2),'-','color',cols(i,:),'LineWidth', i)
end
end
