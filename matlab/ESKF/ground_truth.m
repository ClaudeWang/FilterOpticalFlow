load('data/studentdata1.mat');
all_vicon = vicon;
num_pts = size(all_vicon, 2);

vx_next = all_vicon(7, 2: num_pts);
vx_prev = all_vicon(7, 1: num_pts - 1);
t_diff = time(1, 2: num_pts) - time(1, 1:num_pts - 1);

acc_x = (vx_next - vx_prev) ./ t_diff;

figure(1);
% subplot(2, 1, 1);
% plot(acc_x);
% 
% subplot(2, 1, 2);
% plot(all_vicon(7, :));
subplot(3, 1, 1);
plot(all_vicon(4, :));

subplot(3, 1, 2);
plot(all_vicon(5, :));

subplot(3, 1, 3);
plot(all_vicon(6, :));

