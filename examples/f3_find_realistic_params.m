function  f3_find_realistic_params

masses = linspace(2, 16, 30);
forces = linspace(0.01, 0.2, 30);


targetTrajLength = 50; %mm
targetThresh = 5; %mm
targetTime = 150; %sec


phaseMap = zeros(length(forces),length(masses));
    
% INITIALIZE A WORM:
opts.worm.mutant = 'wt'; 
% 'vfl3' - left
% 'odf2' - right
% 'wt'

% constant parameters:
opts.worm.length = 8; %mm
opts.worm.bodyRatio = 5.0; %length to width, normal value = 2.8
opts.worm.dt = 1; %time step for integration
opts.worm.initialPosition = [0, 0, 0];
opts.worm.orientation = 90; %90 = upwards

opts.simulation.frames = targetTime;%sec
opts.simulation.showFrames = 0;
opts.simulation.showStamps = 0;
opts.simulation.saveVideo = 0;
opts.simulation.VideoName = 'temp_test';


close all;
counter = 0;
N = length(forces)*length(masses);

for i = 1:length(forces);
    for j = 1:length(masses)
        counter = counter+1;
        disp(['-------------  ' num2str(counter)  ' out of '   num2str(N)   '  ---------------'])
        
        opts.worm.strokeForce = forces(i); %arbitrary unit
        opts.worm.mass =  masses(j); %arbitrary

        w = gen_worm(opts);
        plotGraph = 0;
        [traj_i, angles_i] = runStuff(w, opts, plotGraph);

        L=0;
        for k = 1:size(traj_i, 1) - 1
            p1 = traj_i(k, [1 2]);
            p2 = traj_i(k+1, [1 2]);
            dl = sqrt(sum((p2-p1).^2, 2));
            L = L + dl;
        end
        phaseMap(i,j) = L;
        
    end
end
        
%reference_plane = (phaseMap./phaseMap)*targetTrajLength;
%result_surf = abs(phaseMap - reference_plane);
close all;
figure;
% surf(phaseMap - reference_plane)
imagesc(phaseMap)
hold on;

% shading interp
masses_labels = 0.01*round(masses*100);
forces_labels = 0.001*round(forces*1000);

set(gca,'XTick', [1 : length(masses)])
set(gca,'XTickLabel', masses_labels)
set(gca,'YTick', [1 : length(forces)])
set(gca,'YTickLabel', forces_labels)
ylabel('forces')
xlabel('masses')
% view(2)
% zlims = zlim;
% zlim([0 zlims(2)])


checkMat = abs(phaseMap - targetTrajLength) < targetThresh;
[r, c] = find(checkMat == 1);
ll = polyfit(c,r,1);
mm = c(1):c(end);
ff = polyval(ll, mm);
plot(c, r, 'wx')
plot(mm, ff,'b-','LineWidth', 2)

fileName = ['L' num2str(targetTrajLength) '_t' num2str(targetTime) ];
save([fileName '.mat'],'targetTrajLength','targetThresh','targetTime','mm','ff','phaseMap','r','c')

