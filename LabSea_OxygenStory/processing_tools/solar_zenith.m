function [theta_deg,theta] = solar_zenith(time, lat)
    % use UTC time!!! 
    time = time(:);
    phi  = deg2rad(lat(:));

    % fract time hours
    [~,~,~,HH,MM,SS] = datevec(time);
    frac_hours = HH + MM/60 + SS/3600;

    % Hour angle
    hour_angle = deg2rad(15 * (frac_hours - 12));

    % Declination angle
    di = deg2rad(23.45 * sin(2*pi*(doy(time) - 81)/365));

    % Cosine of the solar zenith angle
    cos_theta = sin(phi).*sin(di) + ...
                cos(phi).*cos(di).*cos(hour_angle);

    % Solar zenith angle
    theta = acos(cos_theta);
    theta_deg = rad2deg(theta);
end
% 
% function d = doy(time)
%     % Function to compute day of the year
%     ref = datenum(year(time),1,1);
%     d = time - ref + 1;
% end
