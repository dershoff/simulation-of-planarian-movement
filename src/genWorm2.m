classdef genWorm2 < handle

    properties
        faces;
        vertices;
        outlineAbs;     
        outlineRel;     
        center;
        nodeList; %centers of the faces formed by vertices
        
        bodyAngle; % body orientation
        bodyRatio;
        bodyLength;
        
        mass;
        nodeMass;
        inertia;
        strokeForce;
        totalTorque;
        totalForce;
        dt;

        mutant;       
    end
    
    
    methods     
        
        function wrm = genWorm2(opts)

            wrm.bodyLength = opts.worm.length; %mm
            wrm.bodyRatio = opts.worm.bodyRatio;
            wrm.center = opts.worm.initialPosition'; %column
            wrm.bodyAngle = opts.worm.orientation;
            wrm.mass = opts.worm.mass;
            wrm.dt = opts.worm.dt;
            wrm.mutant = opts.worm.mutant;            
            wrm.strokeForce = opts.worm.strokeForce;           

            %get the nodes; the coordinates rescaled to proer lengths and
            % body-ratios. Faces are sorted according to the Angluar data.
            [verts, facesMat, outlineCoord] =  get_nodes_and_springs(wrm.bodyRatio, wrm.bodyLength, wrm.center);
            
            if 0 
                figure; axis equal; hold on;
                plot(outlineCoord(:,1), outlineCoord(:,2), 'k-')
                plot(verts(:,1), verts(:,2), 'r.')
                plot(facesMat(:,1), facesMat(:,2), 'bo')
                plot(wrm.center(1), wrm.center(2), 'kx')
            end
            
            
            outlineR = sqrt((outlineCoord(:,1) - wrm.center(1)) .^2 + ...
                (outlineCoord(:,2) - wrm.center(2)) .^2);
            outlinePhi = [];
            for i = 1:size(outlineCoord, 1)
                dx = outlineCoord(i,1) -  wrm.center(1);
                dy = outlineCoord(i,2) -  wrm.center(2);
                angleDeg = rad2deg(atan2(dy, dx));                
                outlinePhi = [outlinePhi; angleDeg];
            end

            outlineCoord_rotated    = rotation(outlineCoord(:, [1:2]), [wrm.center(1), wrm.center(2)], - wrm.bodyAngle);
            verts_rotated           = rotation(verts(:, [1:2]), [wrm.center(1), wrm.center(2)], - wrm.bodyAngle);
            faces_rotated           = rotation(facesMat(:, [1:2]), [wrm.center(1), wrm.center(2)], - wrm.bodyAngle);
            
            wrm.outlineRel          = [outlineR(:), outlinePhi(:)];
            wrm.outlineAbs          = [outlineCoord_rotated, zeros(size(outlineCoord_rotated,1),1)];
            wrm.vertices            = [verts_rotated, zeros(size(verts_rotated,1),1)];            
            wrm.faces               = [faces_rotated, zeros(size(faces_rotated,1),1)];
                        
           
             if 0%debugging figures
                 figure; axis equal; hold on;
                 plot(wrm.center(1)         , wrm.center(2)         , 'kx')
                 plot(wrm.outlineAbs(:,1)   , wrm.outlineAbs(:,2)   , 'k-')
                 plot(wrm.vertices(:,1)     , wrm.vertices(:,2)     , 'r.')
                 plot(wrm.faces(:,1)        , wrm.faces (:,2)       , 'bo')
             end
            
            % make the nodes:
            wrm.nodeMass = wrm.mass/(size(wrm.faces, 1));

            for i = 1:size(wrm.faces, 1);           
                p = wrm.faces(i, :)';             
                node_i = genNode2(i, wrm.nodeMass, p,  wrm);
                wrm.nodeList = [wrm.nodeList, node_i];
                
                % recalculate nodes parameters:
                wrm.nodeList(i).update;
            end
            
            % moment of inertia: calculated via |r| (not a vector)
            wrm.inertia = sum(wrm.nodeMass * ([wrm.nodeList.r]).^2);
            
            % wrm.setNodeStrokes(indList, StrokeForceList, StrokeAngleList)
            wrm.setNodesStrokesFromXls;            
               
            % no forces at initiation:
            wrm.totalTorque = wrm.getTotalTorque; %column
            wrm.totalForce = wrm.getTotalForce; %column
            
            % wrm.singleShow;

        end
        
        
        
        
        
        function setNodesStrokesFromXls(wrm)
            % Set default strokes:
            A = get_angles_from_xls(wrm.mutant);            
            StrokeAngleList = A(:);
            StrokeForceList = wrm.strokeForce * ones(1, length(StrokeAngleList));
                        
            for ind = [wrm.nodeList.ind]                
                wrm.nodeList(ind).setStrokeStrength(StrokeForceList(ind));
                wrm.nodeList(ind).setStrokeDirection(StrokeAngleList(ind));
                wrm.nodeList(ind).update;
            end      
        end
        
        
        
        
        function setNodeStrokes(wrm, indList, StrokeStrengthList, StrokeDirList)            
            for i = 1:length(indList)
                ind = indList(i);
                wrm.nodeList(ind).setStrokeStrength(StrokeStrengthList(i));
                wrm.nodeList(ind).setStrokeDirection(StrokeDirList(i));
                wrm.nodeList(ind).update;
            end            
        end
        
        
        function F_total = getTotalForce(wrm)
            xyForces = [wrm.nodeList.generatedForce];
            F_total = sum(xyForces, 2);
        end
         function T_total = getTotalTorque(wrm)
             xyTorques = [wrm.nodeList.generatedTorque];
            T_total = sum(xyTorques, 2);
        end
        
        
        function move(wrm)
            
            DT = wrm.dt;
            
            % movement is due to integrated effect of all strokes:
            F_total = wrm.getTotalForce;
            
            center_accel = F_total / wrm.mass;
            center_speed = [0; 0; 0];%????
            wrm.center = wrm.center + center_speed * DT + 0.5 * center_accel * DT^2;
            wrm.totalForce = F_total;
            
           % torque is due to the total effect of all strokes 
            T_total = wrm.getTotalTorque;
            
            % totalT = I*w; I - moment of inertia; w - angluar acceleration
            angular_accel = norm(T_total) / wrm.inertia;
            if T_total(3) > 0
                angular_accel = +abs(angular_accel);
            else
                angular_accel = -abs(angular_accel);
            end
            angluar_speed = [0];%????
            
            % PHI = PHI0 + wt
            wrm.bodyAngle = wrm.bodyAngle + angluar_speed * DT + 0.5 * angular_accel * DT^2;
            
            wrm.totalTorque = T_total;
            
            %renew node coordinates:
            for i = [wrm.nodeList.ind]
                wrm.nodeList(i).update; 
                %wrm.verts
                %wrm.outline
            end
            
            wrm.outlineUpdate;
        end
        
        
        
        function outlineUpdate(wrm)
            R = wrm.outlineRel(:, 1);
            Phi = wrm.outlineRel(:, 2);
            wrm.outlineAbs = [wrm.center(1) + R .* cosd(Phi + wrm.bodyAngle), ...
             wrm.center(2) + R .* sind(Phi + wrm.bodyAngle),...
             zeros(size(R,1),  1)];
        end
        
        
        
        
        function bodyMark = stamp(wrm, col)
            % show outline:
            xy = [wrm.outlineAbs];
            %plot3(xy(:,1), xy(:,2), xy(:,3),'-', 'color', col)  
                     
            % show center
            c = wrm.center;
            %plot3(c(1), c(2), c(3), 'rx')
            
            bodyMark = [xy, repmat(c', size(xy,1), 1), repmat(col, size(xy,1), 1)];
        end
        
        
        
        function show(wrm) 
            
            % show outline:
            xy = [wrm.outlineAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'k-')  
            
            % show all verts:
            % xy = [wrm.verts];
            % plot3(xy(:,1), xy(:,2), xy(:,3), 'b.') 
            
            % show center
            c = wrm.center;
            plot3(c(1), c(2), c(3), 'ro')
            
            % show torque direction (normalized to 1)
            T_to_show = c + wrm.totalTorque/norm(wrm.totalTorque);            
            line([c(1), T_to_show(1)], [c(2), T_to_show(2)], [c(3), T_to_show(3)], 'color', 'r','LineWidth', 1)
%             quiver3(c(1), c(2), c(3),  wrm.torque(1),  wrm.torque(2),  wrm.torque(3), 'color', 'r','LineWidth', 1)

            % show Force direction (normalized to 1)
            F_to_show = c + wrm.totalForce/norm(wrm.totalForce);            
%             line([c(1), F_to_show(1)], [c(2), F_to_show(2)], [c(3), F_to_show(3)], 'color', 'g','LineWidth', 1)
            quiver3(c(1), c(2), c(3),  wrm.totalForce(1),  wrm.totalForce(2),  wrm.totalForce(3), 'color', 'g','LineWidth', 1)
            
            % show nodes:
            xyPos = [wrm.nodeList.pAbs];
            plot3(xyPos(1,:), xyPos(2,:), xyPos(3,:), 'b.')            
           
            % show strokes:
            xyStroke = [wrm.nodeList.generatedForce];
            xyEnds  = xyPos + xyStroke;
            strXX = [xyPos(1,:); xyEnds(1,:)];
            strYY = [xyPos(2,:); xyEnds(2,:)];
            strZZ = [xyPos(3,:); xyEnds(3,:)];
            if any(any(xyStroke ~= 0))   
%                 quiver3(xyPos(1,:), xyPos(2,:), xyPos(3,:), xyStroke(1,:), xyStroke(2,:), xyStroke(3,:), 'color','m')
                line(strXX, strYY, strZZ, 'color', 'm')
            end
            
            % show torques:
            xyTorque = [wrm.nodeList.generatedTorque];
            xyEnds  = xyPos + xyTorque;
            torXX = [xyPos(1,:); xyEnds(1,:)];
            torYY = [xyPos(2,:); xyEnds(2,:)];
            torZZ = [xyPos(3,:); xyEnds(3,:)];            
            if any(any(xyTorque ~= 0))   
                line(torXX, torYY, torZZ, 'color', 'g')
%                 quiver3(xyPos(1,:), xyPos(2,:), xyPos(3,:), xyTorque(1,:), xyTorque(2,:), xyTorque(3,:), 'color','k')
            end
            
        end
        
        
        
        function singleShow(wrm)  
            %close all;
            fig_handle = figure('position', [680    74   785   904]);
            title({['Mutant: ' wrm.mutant],['Stroke Force: ' num2str(wrm.strokeForce)]})
            axis equal
            hold on
            
            % show outline:
            xy = [wrm.outlineAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'k-')  
            
            % show all verts:
            xy = [wrm.vertices];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'b.') 
            
            % show center
            center_xyz = wrm.center;
            plot3(center_xyz(1), center_xyz(2), center_xyz(3), 'rx','markersize', 14, 'linewidth', 2)
            
            % show all nodes:
            xyPos = [wrm.nodeList.pAbs];
            plot3(xyPos(1,:), xyPos(2,:), xyPos(3,:), 'ro') 
            
            % show strokes:
            xyStroke = [wrm.nodeList.generatedForce];
            xyEnds  = xyPos + xyStroke;
            strXX = [xyPos(1,:); xyEnds(1,:)];
            strYY = [xyPos(2,:); xyEnds(2,:)];
            strZZ = [xyPos(3,:); xyEnds(3,:)];
            if any(any(xyStroke ~= 0))   
                %quiver(xyPos(1,:), xyPos(2,:), xyEnds(1,:), xyEnds(2,:), 0)
                line(strXX, strYY, strZZ, 'color', 'm')
            end
            
            % show torques:
            xyTorque = [wrm.nodeList.generatedTorque];
            xyEnds  = xyPos + xyTorque;
            torXX = [xyPos(1,:); xyEnds(1,:)];
            torYY = [xyPos(2,:); xyEnds(2,:)];
            torZZ = [xyPos(3,:); xyEnds(3,:)];            
            if any(any(xyTorque ~= 0))   
                line(torXX, torYY, torZZ, 'color', 'g')
            end
            
            % show total torque, normalized
            T = center_xyz + wrm.totalTorque/norm(wrm.totalTorque);            
            line([center_xyz(1), T(1)], [center_xyz(2), T(2)], [center_xyz(3), T(3)],...
                'color', 'r','LineWidth',2) 
            
            % show Force direction, normalized.
            F_to_show = center_xyz + wrm.totalForce/norm(wrm.totalForce);            
            line([center_xyz(1), F_to_show(1)], [center_xyz(2), F_to_show(2)], [center_xyz(3), F_to_show(3)],...
                'color', 'g','LineWidth', 1)            
            
            
        end
        
             

    end
end



function rotated_coords = rotation(input_XY, center, anti_clockwise_angle, varargin)


%                               431-400 Year Long Project 
%                           LA1 - Medical Image Processing 2003
%  Supervisor     :  Dr Lachlan Andrew
%  Group Members  :  Alister Fong    78629   a.fong1@ugrad.unimelb.edu.au
%                    Lee Siew Teng   102519  s.lee1@ugrad.unimelb.edu.au
%                    Loh Jien Mei    103650  j.loh1@ugrad.unimelb.edu.au
% 
%  File and function name : rotation
%  Version                : 1.0
%  Date of completion     : 5 September 2003
%  Written by    :   Alister Fong    78629   a.fong1@ugrad.unimelb.edu.au
%
%  Inputs        :
%           input_XY    -   The X and Y coordinates to be rotated. [X,Y]
%           center      -   A 1D matrix with 2 formats either 
%                           [centerX,centerY] or [centerX;centerY]
%                           this is the center of rotation
%           anti-clockwise_angle - The angle of rotation about the center 
%                                  from which the new coordinates will be produced.
%           'degree','radians' (Optional)
%                       -   Description of what the "anti-clockwise_angle" is.
%                           Default = 'degree'
%
%  Outputs       :   
%           rotated_coords     - The rotated coordinates "input_XY" about "center" 
%                                by angle "anti-clockwise_angle". [X,Y]
%
%  Description   :
%       Rotates the input_XY coordinates by a given angle about a center.
%
%  To Run >>    final_coords = rotation(input_XY,center,anti-clockwise_angle,varargin)
%
%  Example>>    X = [1:1:10]'; Y = [2;4;1;2;5;7;2;7;8;3];
%               final_coords = rotation([X,Y],[0,0],90);
%               figure;
%               plot(X,Y,'b+-');
%               hold on;
%               plot(final_coords(:,1),final_coords(:,2),'r*-');
%
% See "test_rotation.m" for more examples

degree = 1; %Radians : degree = 0; Default is calculations in degrees

% Process the inputs
if length(varargin) ~= 0
for n = 1:1:length(varargin)
    if strcmp(varargin{n},'degree') 
        degree = 1;
    elseif strcmp(varargin{n},'radians')
        degree = 0;
    end
end
clear n;
end
[r,c] = size(input_XY);
if c ~= 2
error('Not enough columns in coordinates XY ');
end
[r,c] = size(center);
if (r~=1 & c==2) | (r==1 & c~=2)
error('Error in the size of the "center" matrix');
end

% Format the coordinate of the center of rotation
center_coord = input_XY;
center_coord(:,1) = center(1);
center_coord(:,2) = center(2);

% Turns the angles given to be such that the +ve is anti-clockwise and -ve is clockwise
anti_clockwise_angle = -1*anti_clockwise_angle;
% if in degrees, convert to radians because that's what the built-in functions use. 
if degree == 1 
anti_clockwise_angle = deg2rad(anti_clockwise_angle);
end

%Produce the roation matrix
rotation_matrix = [cos(anti_clockwise_angle),-1*sin(anti_clockwise_angle);...
               sin(anti_clockwise_angle),cos(anti_clockwise_angle)];
%Calculate the final coordinates
rotated_coords = ((input_XY-center_coord) * rotation_matrix)+center_coord;

end