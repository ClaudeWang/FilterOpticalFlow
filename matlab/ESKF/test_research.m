clear ESKF estimate_vel_research
data = load('data/studentdata1.mat');
all_x = zeros([864, 1]);
all_y = zeros([864, 1]);
all_z = zeros([864, 1]);

all_o = zeros([864, 3]);
all_vx = zeros([864, 1]);

all_accx = zeros([864, 1]);
all_accx_body = zeros([864, 1]);

all_rpy = zeros([864, 3]);
all_rpy_estimated = zeros([864, 3]);

ts = zeros([864, 1]);


fileID = fopen('data.txt','w');
fprintf(fileID, '%d\n', 864);

til = 864;
for i = 1:til

%     [v, omg, flow, points, depth] = estimate_vel_research(data.data(i));
     [state, R] = ESKF(data.data(i));
%     all_x(i) = state(1);
%     all_y(i) = state(2);
%     all_z(i) = state(3);
    
%     all_vx(i) = state(8);
%       all_accx(i) = acc(1);
%       all_accx_body(i) = data.data(i).acc(1);
      
%      all_o(i, :) = data.data(i).omg;
      all_rpy(i, :) = data.data(i).rpy';
      all_rpy_estimated(i, :) = rotm2eul(R, 'XYZ');
%     if size(points, 1) == 0
%         % ignore the frame
%         'no flow'
%         continue
%     end
    
%     fprintf(fileID, '%d:\n', size(points, 2));
%     data_print = [flow; points; depth]';
%     
%     for j = 1:size(data_print, 1)
%         one_row = data_print(j, :);
%         fprintf(fileID, '%f %f %f %f %f %f', one_row(1), one_row(2), one_row(3), one_row(4), one_row(5));
%         fprintf(fileID, '\n');
%     end

end

% fclose(fileID);

figure(2);
% 
% subplot(3, 1, 1);
% plot(all_o(:, 1));
% 
% subplot(3, 1, 2);
% plot(all_o(:, 2));
% 
% subplot(3, 1, 3);
% plot(all_o(:, 3));






subplot(3, 1, 1);
plot(all_rpy(1:til, 1));
hold on;
plot(all_rpy_estimated(1:til, 1));

subplot(3, 1, 2);
plot(all_rpy(1:til, 2));
hold on;
plot(all_rpy_estimated(1:til, 2));

subplot(3, 1, 3);
plot(all_rpy(1:til, 3));
hold on;
plot(all_rpy_estimated(1:til, 3));


% figure(1)
% subplot(3,1,1)
% plot(all_x)
% 
% subplot(3,1,2)
% plot(all_y)
% % 
% subplot(3,1,3)
% plot(all_z)


% figure(2)
% subplot(3,1,1)
% plot(ts, all_ox)
% 
% subplot(3,1,2)
% plot(ts,all_oy)
% 
% subplot(3,1,3)
% plot(ts,all_oz)



% subplot(2,1,2)
% plot(num_points_left)


% minQuality = 0.2;
% nPyrLevels = 4;
% maxBiError = 1;
% 
% img0 = sensor.img;
% points = detectMinEigenFeatures(img0,'MinQuality',minQuality);
% point_tracker = vision.PointTracker('NumPyramidLevels',nPyrLevels,...
%                                    'MaxBidirectionalError',maxBiError);
% initialize(point_tracker, points.Location, img0);
% last_points = points.Location;
% last_t = sensor.t;
% vel = [0, 0, 0];
% omg = [0, 0, 0];





