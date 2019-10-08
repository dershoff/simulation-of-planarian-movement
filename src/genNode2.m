classdef genNode2 < handle
    
    
    properties  (SetAccess = private)%, Hidden = true)
        % cannot change:
        ind;
        m;
        strokeStrength;
        strokeDirection;  
        parent; 
        pAbs; % array 3x1 
        generatedForce; %array 2x1
        generatedTorque;
    end
    
%     properties         % CAN CHANGE
%     end
    
     methods
         function nd = genNode2(ind, node_mass, xyz, parent)
             nd.ind = ind;
             nd.m = node_mass;   
             nd.parent = parent;    
             nd.pAbs = xyz;
             nd.strokeStrength = 0;
             nd.strokeDirection = 0; %0 is parallel to the body orientation
             nd.generatedForce = [0; 0; 0]; %vector 2x1
             nd.generatedTorque = [0; 0; 0];
         end 
         
         function setStrokeStrength(nd, str_stren)
             nd.strokeStrength = str_stren;
         end
         
         function setStrokeDirection(nd, str_dir)
             nd.strokeDirection = str_dir;
         end
         
         
         function  update(nd)
             % Recalculate node position and direciton of generated force.
             % Recalculate after each parent movement!
             
             % recalculate Absolute node position from Relative-to-center position:
             nd.pAbs = nd.parent.faces(nd.ind, :)';
           
             % STROKE:
             % recalculate direction of the stroke according to the new parent
             % position
             parentMainOrientaton = nd.parent.bodyAngle;
             nodeRelativeOrientation = nd.strokeDirection;
             nodeAbsoluteOrientation = parentMainOrientaton + nodeRelativeOrientation;

             strokeStr = nd.strokeStrength;
             nd.generatedForce = strokeStr*[ cosd(nodeAbsoluteOrientation); sind(nodeAbsoluteOrientation); 0];
             
             % TORQUE:
             force = nd.generatedForce;
             try
                 radius = nd.pAbs - nd.parent.center;
             catch
                 disp(' ')
             end
             T = cross(radius, force);
             nd.generatedTorque = T;
         end
         
     end
end