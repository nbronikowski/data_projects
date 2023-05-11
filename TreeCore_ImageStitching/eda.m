

% treethick = [3 3 2 1 6 3];
% rain = [0.5 .2 .6 1.3 0.6 0.2];
% temp = rand(6,1)
% 
% scatter3(rain,temp,20,'m','filled')


depth=[0.5,0.7,1,2,3,4,5,6,7,8,9,10,9,10,9,8,7,6,5,4,3,2,1.5,0.9,0];
time = linspace(0,length(depth)-1,length(depth));
time=time(:);
depth=depth(:);
[profile_index,profile_dir]=findProfiles(time,depth,'length',9,'stall',0,'shake',0)

figure()
idx_down = profile_dir == 1;
idx_up   = profile_dir == -1;
hold on;
plot(time(idx_down),-depth(idx_down),'*b')
plot(time(idx_up),-depth(idx_up),'*r')
xlabel('time(s)')
ylabel('depth(m)')
title('Algorithm for Profile Splitting Code')
grid on
legend('down cast','up cast','location','best')
