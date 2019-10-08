function [verts_shifted, faces_shifted, outline] = get_nodes_and_springs(bodyRatio, bodyLength, bodyCenter)
% generate worm from blender:
file_path = 'worm_blender.txt';

f_id = fopen(file_path, 'r');
what_read = '';

vertices = [];
edges = [];
faces = [];

while 1
    tline = fgetl(f_id);
    if ~ischar(tline), break, end
    
    if strfind(tline, 'vertices')
        what_read = 'vertices';
    elseif strfind(tline, 'edges')
         what_read = 'edges';  
     elseif strfind(tline, 'faces')
         what_read = 'faces';
    elseif strfind(tline, 'bodies')
        break;
    end
    
    if ischar(tline)
        nums_i = str2num(tline);
        if ~isempty(nums_i)
            eval([what_read ' = [' what_read '; nums_i];'])        
        end
        %         disp(tline)
    end
    
end

fclose(f_id);

% remove numbering in vertices and remove z coord:
% vertices = vertices(:, 1:3);


% shift y to 0, rescale to 1:
col_ind = 3; % ind for y column
miny = min(vertices(:, col_ind));
vertices(:, col_ind) = vertices(:, col_ind) - miny;
maxy = max(vertices(:, col_ind));
vertices(:, col_ind) = vertices(:, col_ind)/maxy;

% rescale x to 1: 
col_ind = 2; % ind for y column
minx = min(vertices(:, col_ind));
vertices(:, col_ind) = vertices(:, col_ind) - minx;
maxx = max(vertices(:, col_ind));
vertices(:, col_ind) = vertices(:, col_ind)/maxx;
% scale x to proper width:      
vertices(:, col_ind) = vertices(:, col_ind)/bodyRatio; 



% apply bodyLength
%for x
vertices(:, 2) = vertices(:, 2)*bodyLength;
%for y
vertices(:, 3) = vertices(:, 3)*bodyLength; 


faces_centers = [];
for i = 1 : size(faces,1)
    
    edges_i = abs(faces(i, 2:end));
    verts_i = unique(edges(edges_i, 2:3));%has to be 4 nodes for square
     
    coords_i = vertices(verts_i,2:3);
    face_mean_coord = mean(coords_i, 1);
    
    faces_centers = [faces_centers; face_mean_coord];
end


% figure; hold on
% axis equal;
% plot(vertices(:,2),vertices(:,3),'k.')
% plot(faces_centers(:,1),faces_centers(:,2),'r.')


faces_centers_sorted = nan(size(faces_centers));
 
fx = faces_centers(:, 1);
fy = faces_centers(:, 2); 

x_unique = unique(fx)';
% 
% [val, indx] = sort(fx);
counter = 0;
for x = x_unique
    [indx, ~] = find(fx == x);
    y = fy(fx == x);
    [valy, indy] = sort(y);
    indx = indx(indy);    
    faces_centers_sorted(counter + 1 : counter + 5, :) = faces_centers(indx,:);
    counter = counter + 5;
end



% add z-dimension to them:
verts = [vertices(:, 2:3), zeros(size(vertices, 1), 1)];% add z=0;
faces = [faces_centers_sorted, zeros(size(faces_centers_sorted, 1), 1)];



 
% first shift all data so that the mean is around the desired position (= bodyCenter):
old_center = [mean(verts(:, 1));  mean(verts(:, 2)); 0];            

% find the required shift:
general_shift  = bodyCenter - old_center;

% shift all verts:
verts_shift = repmat(general_shift', size(verts, 1), 1);
verts_shifted = verts + verts_shift;

% shift faces:
faces_shift = repmat(general_shift', size(faces, 1), 1);
faces_shifted = faces + faces_shift; 

% outline shifted:
x = verts_shifted(:, 1); y = verts_shifted(:, 2); z = zeros(size(x));
k = convhull(x,y); 
outline = [x(k), y(k), z(k)];






