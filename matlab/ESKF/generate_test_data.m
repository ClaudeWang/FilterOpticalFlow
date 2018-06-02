clear ESKF estimate_vel_research
data = load('data/studentdata1.mat');
all_x = zeros([864, 1]);
all_y = zeros([864, 1]);
all_z = zeros([864, 1]);

all_o = zeros([864, 3]);
all_vx = zeros([864, 1]);
all_vy = zeros([864, 1]);
all_vz = zeros([864, 1]);

all_accx = zeros([864, 1]);
all_accx_body = zeros([864, 1]);

all_rpy = zeros([864, 3]);
all_rpy_estimated = zeros([864, 3]);

ts = zeros([864, 1]);


all_vx_nom = zeros([864, 1]);
all_vy_nom = zeros([864, 1]);
all_vz_nom = zeros([864, 1]);

fileID = fopen('data_with_imu.txt','w');
fprintf(fileID, '%d\n', 864);

til = 864;
for i = 1:til

    [v, omg, ~, flow, points, depth] = estimate_vel_research(data.data(i));
      
    if size(points, 1) == 0
        % ignore the frame
        'no flow'
        fprintf(fileID, '%d:\n', 0);
        fprintf(fileID, '%f %f %f', data.data(i).acc);
        fprintf(fileID, '%f %f %f', data.data(i).omg);
        fprintf(fileID, '\n');
        continue
    end
    
    fprintf(fileID, '%d:\n', size(points, 2));
    data_print = [flow; points; depth;]';
    
    for j = 1:size(data_print, 1)
        one_row = data_print(j, :);
        fprintf(fileID, '%f %f %f %f %f', one_row(1), one_row(2), one_row(3), one_row(4), one_row(5));
        fprintf(fileID, '\n');
    end
    
    fprintf(fileID, '%f %f %f', data.data(i).acc);
    fprintf(fileID, '%f %f %f', data.data(i).omg);
    fprintf(fileID, '\n');

end

fclose(fileID);
