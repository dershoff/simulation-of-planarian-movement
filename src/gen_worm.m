classdef gen_worm < handle

    properties
        
        facesRelCoord;
        vertsRelCoord;
        outlineRelCoord;
        outlineAbs;
        verticesAbs;
        facesAbs;
        
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
        
        function wrm = gen_worm(opts)

            wrm.bodyLength      = opts.worm.length; %mm
            wrm.bodyRatio       = opts.worm.bodyRatio;
            wrm.center          = opts.worm.initialPosition'; %column 
            wrm.bodyAngle       = opts.worm.orientation;
            wrm.mass            = opts.worm.mass;
            wrm.dt              = opts.worm.dt;
            wrm.mutant          = opts.worm.mutant;            
            wrm.strokeForce     = opts.worm.strokeForce;           
 
            %get the nodes; the coordinates rescaled to proer lengths and
            % body-ratios. Faces are sorted according to the Angluar data.
            [vertsCoord, facesCoord, outlineCoord] =  get_nodes_and_springs(wrm.bodyRatio, wrm.bodyLength, wrm.center);            
                 
            % get the relative Polar coordinates for all items:
            wrm.facesRelCoord       = getRelPolarCoord(facesCoord, wrm.center);
            wrm.vertsRelCoord       = getRelPolarCoord(vertsCoord, wrm.center);
            wrm.outlineRelCoord     = getRelPolarCoord(outlineCoord, wrm.center); 
            
            %rotate to 0:
            vertsCoord              = rotation(vertsCoord(:, [1:2]), [wrm.center(1), wrm.center(2)], - 90);
            facesCoord              = rotation(facesCoord(:, [1:2]), [wrm.center(1), wrm.center(2)], - 90);
            outlineCoord            = rotation(outlineCoord(:, [1:2]), [wrm.center(1), wrm.center(2)], - 90);
            
            
            vertsZeros              = zeros(size(vertsCoord, 1), 1);
            facesZeros              = zeros(size(facesCoord, 1), 1);
            outlineZeros            = zeros(size(outlineCoord, 1), 1);
       
            % Rotate the coordinates according to the body angle:
            wrm.outlineAbs          = [rotation(outlineCoord(:, [1:2]), [wrm.center(1), wrm.center(2)], wrm.bodyAngle),   outlineZeros];
            wrm.verticesAbs         = [rotation(vertsCoord(:,   [1:2]), [wrm.center(1), wrm.center(2)], wrm.bodyAngle),   vertsZeros];
            wrm.facesAbs            = [rotation(facesCoord(:,   [1:2]), [wrm.center(1), wrm.center(2)], wrm.bodyAngle),   facesZeros];
                        
           
            if 0 %debug
                nn = 2;
                                
                figure('position',[    680         124        1003         854]);
                title('input data')
                axis equal; hold on;
                plot(opts.worm.initialPosition(1), opts.worm.initialPosition(2), 'kx')                
                plot(outlineCoord(:,1), outlineCoord(:,2), 'k-')
                plot(vertsCoord(:,1),   vertsCoord(:,2), 'r.')
                plot(facesCoord(:,1),   facesCoord(:,2), 'bo')  
                %show some angles and R:
                plot(outlineCoord(nn,1), outlineCoord(nn,2),  'ks')
                plot( vertsCoord(nn,1), vertsCoord(nn,2), 'rs')
                plot( facesCoord(nn,1),  facesCoord(nn,2),   'bs')  
                text( outlineCoord(nn,1),   outlineCoord(nn,2),  num2str(wrm.outlineRelCoord(nn,:)))
                text( vertsCoord(nn,1),     vertsCoord(nn,2),    num2str(wrm.vertsRelCoord(nn,:)))
                text( facesCoord(nn,1),     facesCoord(nn,2),    num2str(wrm.facesRelCoord(nn,:))) 
                
                
                figure('position',[    680         124        1003         854]);
                title('rotated data')
                axis equal; hold on;
                plot(opts.worm.initialPosition(1), opts.worm.initialPosition(2), 'kx')                
                plot(wrm.outlineAbs(:,1),  wrm.outlineAbs(:,2), 'k-')
                plot( wrm.verticesAbs(:,1),    wrm.verticesAbs (:,2), 'r.')
                plot( wrm.facesAbs (:,1),    wrm.facesAbs (:,2), 'bo')  
                
                plot( wrm.outlineAbs(nn,1),     wrm.outlineAbs(nn,2),   'ks')
                plot( wrm.verticesAbs(nn,1),    wrm.verticesAbs(nn,2),  'rs')
                plot( wrm.facesAbs (nn,1),      wrm.facesAbs (nn,2),    'bs')  
                text( wrm.outlineAbs(nn,1),     wrm.outlineAbs(nn,2),   num2str( wrm.outlineRelCoord(nn,:)))
                text( wrm.verticesAbs(nn,1),    wrm.verticesAbs(nn,2),  num2str( wrm.vertsRelCoord(nn,:)))
                text( wrm.facesAbs(nn,1),       wrm.facesAbs(nn,2),     num2str(wrm.facesRelCoord(nn,:))) 
            end                
      
   
            % generate the nodes:
            Nnodes = size(wrm.facesAbs, 1);
            wrm.nodeMass = wrm.mass/Nnodes;

            for i = 1 : Nnodes         
                p = wrm.facesAbs(i, :);
                %column with z=0:
                p = [p(:); 0];
                
                r   =  wrm.facesRelCoord(i, 1);
                phi =  wrm.facesRelCoord(i, 2);                
                
                % generate node:
                node_i = gen_node(i, wrm.nodeMass, p, r, phi,  wrm);
                wrm.nodeList = [wrm.nodeList, node_i];
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
            
            %update coordinates that changed due to the changes in the
            %center and angle
            wrm.coordsUpdate;
        end
        
        
        
        function coordsUpdate(wrm)
            
            % get the abs coords from Polar coordinates for all items:
            wrm.facesAbs       = getAbsFromPolar(wrm.facesRelCoord,     wrm.facesAbs,       wrm.center,     wrm.bodyAngle);
            wrm.verticesAbs    = getAbsFromPolar(wrm.vertsRelCoord,     wrm.verticesAbs,    wrm.center,     wrm.bodyAngle);
            wrm.outlineAbs     = getAbsFromPolar(wrm.outlineRelCoord,   wrm.outlineAbs,     wrm.center,     wrm.bodyAngle); 
                               
              % renew node coordinates:
            for i = [wrm.nodeList.ind]
                wrm.nodeList(i).update; 
            end
        
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
        
        
        
        function addToFig(wrm, figHandle)  
            %close all;
%             fig_handle = figure('position', [680    74   785   904]);
%             title({['Mutant: ' wrm.mutant],['Stroke Force: ' num2str(wrm.strokeForce)]})
%             axis equal
%             hold on
            figure(figHandle)
            
            % show outline:
            xy = [wrm.outlineAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'k-')  
            
            % show all verts:
            xy = [wrm.verticesAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'r.') 
            
            % show all faces:
            xy = [wrm.facesAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'bo') 
            % show all nodes:
            xyPos = [wrm.nodeList.pAbs];
            plot3(xyPos(1,:), xyPos(2,:), xyPos(3,:), 'm+') 
            
            % show center
            center_xyz = wrm.center;
            plot3(center_xyz(1), center_xyz(2), center_xyz(3), 'rx','markersize', 14, 'linewidth', 2)
            
           
            
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
            
            % show  Total, normalized
            T = center_xyz + wrm.totalTorque/norm(wrm.totalTorque);            
            line([center_xyz(1), T(1)], [center_xyz(2), T(2)], [center_xyz(3), T(3)],...
                'color', 'r','LineWidth',2) 
            
            % show Force direction, normalized.
            F_to_show = center_xyz + wrm.totalForce/norm(wrm.totalForce);            
            line([center_xyz(1), F_to_show(1)], [center_xyz(2), F_to_show(2)], [center_xyz(3), F_to_show(3)],...
                'color', 'g','LineWidth', 1)            
            
            
        end
        
             
        
        
        
        function singleShow(wrm)  
            fig_handle = figure('position', [680    74   785   904]);
            title({['Mutant: ' wrm.mutant],['Stroke Force: ' num2str(wrm.strokeForce)]})
            axis equal
            hold on
            
            % show outline:
            xy = [wrm.outlineAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'k-')  
            
            % show all verts:
            xy = [wrm.verticesAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'r.') 
            
            % show all faces:
            xy = [wrm.facesAbs];
            plot3(xy(:,1), xy(:,2), xy(:,3), 'bo') 
            % show all nodes:
            xyPos = [wrm.nodeList.pAbs];
            plot3(xyPos(1,:), xyPos(2,:), xyPos(3,:), 'm+') 
            
            % show center
            center_xyz = wrm.center;
            plot3(center_xyz(1), center_xyz(2), center_xyz(3), 'rx','markersize', 14, 'linewidth', 2)
            
           for n = [wrm.nodeList]
               text(n.pAbs(1), n.pAbs(2), num2str(n.ind))
           end
            
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
            
            % show  Total, normalized
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




function RelPolCoord = getRelPolarCoord(xyzCoord, center)
R = [];
phis = [];
xc = center(1);
yc = center(2);
for i = 1:size(xyzCoord, 1)
    xi = xyzCoord(i, 1);
    yi = xyzCoord(i, 2);
    dx = xi-xc; 
    dy = yi-yc;   
    ri = sqrt(dx^2 + dy^2);
    angleDeg = rad2deg(atan2(dy, dx));
    phis = [phis; angleDeg];
    R = [R; ri];
end
RelPolCoord = [R, phis];
end


function  xyzNew = getAbsFromPolar(RelPolCoord, xyzOld, c, a)
    
a = 90 + a;
r           = RelPolCoord(:, 1);
phi         = RelPolCoord(:, 2);
xyzNew    = [c(1) + r .* cosd(phi + a), ...
               c(2) + r .* sind(phi + a),...
               zeros(size(r, 1),  1)];
           
% plot(xyzOld(:,1), xyzOld(:,2), 'k.')
% plot(xyzNew(:,1), xyzNew(:,2), 'g.')

end


