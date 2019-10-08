function A = get_angles_from_xls(what_mutant)

data_path  = 'angles98.xls';

switch what_mutant
    case 'wt'
        data_sheet = 'wt';
    case 'vfl3'
        data_sheet = 'vfl3';
    case 'odf2'     
        data_sheet = 'odf2';
    otherwise
        warning('unknown mutant. Using WT.')
        data_sheet = 'wt';
end


% tail - top; head - bottom;
% -0.9 - left side; +0.9 - right side.
num = xlsread(data_path, data_sheet);

A = num(2:end, 2:end);
%copy the side extremeties (they are normally unavailable)
A(:,1) = A(:,2);
A(:,end) = A(:,end-1);

A = flipud(A); %top - head; bottom - tail;
% plot(A');


