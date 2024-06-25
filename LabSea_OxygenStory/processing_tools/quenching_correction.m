function thomalla_correction = quenching_correction(angle,dt,chlorophyll,bbp)
%The thomalla function corrects the non-photochemical quenching (NPQ)
%signal for vertical profiles taken during the day time based upon the relationship
%between backscatter and chlorophyll laid out in Thomalla et al 2017 : npq_corrected_day_chl = chl_night/bbp_night*bbp_day . 
%This correction requires both day and night chlorophyll and backscatter profiles.
%
% Inputs: 
% angle       - Zenith angle of the sun when the profile was taken
% dt          - time of the profile, this should either be inputted as a datetime or
%              a datenum
% chlorophyll - Matrix of chlorophyll observations where each row is a
%              unique depth and each column is a single profile
% bbp         - Matrix of backscatter observations where each row is a
%               unique depth and each column is a single profile.
%               bbp should be the same size as chlorophyll. 
%
% Outputs:
% Structure which contains profiles seperated out in to day and night,
% paired by date and the thomalla correction for the daytime observations

    %check to make sure the inputted profile times are datetimes. If they are
    %date nums, convert them. Elsewise, return an error. 
    if isdatetime(dt) == 0
        try
            dt = datetime(dt,'ConvertFrom','datenum');
        catch
            error('Please enter cast time as a datetime or a datenum')
        end
    end
    
    %pair each day's profiles with night profiles that are taken 
    count = 1;
    pair = zeros(length(angle),1);
    
    for k = 1:length(angle)
        if k == 1 %first setting
            pair(k) = 1;
            continue;
        end
        if angle(k)< 90 && angle(k-1)>= 90 %night turns to day
            count = count+1;
            pair(k) = count;
        end
        if angle(k)< 90 && angle(k-1)< 90 %day to day
            if day(dt(k)) == day(dt(k-1))
                pair(k) = count;
            else
                count = count+1;
                pair(k) = count;
            end
        end
        if angle(k)>= 90 && angle(k-1)< 90 % day turns to night
            if day(dt(k)) == day(dt(k-1))
                pair(k) = count;
            elseif hour(dt(k))<12 && day(dt(k)) == day(dt(k-1))+1
                pair(k) = count;
            else
                count = count+1;
                pair(k) = count;
            end
        end
        if angle(k)>= 90 && angle(k-1)>= 90 %night to night
            %check if same day or next day
            if day(dt(k)) == day(dt(k-1))
                pair(k) = count;
            elseif day(dt(k)) == day(dt(k-1))+1
                pair(k) = count;
            else
                count = count+1;
                pair(k) = count;
            end
        end
    end
    
    % struct per D/N pair, w/ Backscatter for N & D, Chlorophyll for N/D,
    % multiple distinct day profiles, night profiles aggregated; vars exist
    % Need corrected profile based on the adjust formula below
    % D/N pairs are made up some number of profiles, each profile is a column vector
    miss = [];
    fn = 0;
    
    for k = 1:length(pair)
        name = dt(k);
        name.Format = 'dd-MM-yyyy';
         if fn == 0 && angle(k) > 90
            fn = 1;
         end
        if k == 1
            thomalla_correction(pair(k)).name = string(name); % grab the date
            thomalla_correction(pair(k)).bbp_day = bbp(:,k);
            thomalla_correction(pair(k)).bbp_night = NaN(size(bbp,1),1);
            thomalla_correction(pair(k)).chl_day = chlorophyll(:,k);
            thomalla_correction(pair(k)).chl_night = NaN(size(chlorophyll,1),1);
            continue
        end
        if pair(k) == pair(k-1) % Same day night pair
            if angle(k) < 90 % day
                thomalla_correction(pair(k)).bbp_day = [thomalla_correction(pair(k)).bbp_day bbp(:,k)];
                thomalla_correction(pair(k)).chl_day = [thomalla_correction(pair(k)).chl_day chlorophyll(:,k)];
            elseif fn == 1 % first night
                thomalla_correction(pair(k)).bbp_night = [bbp(:,k)];
                thomalla_correction(pair(k)).chl_night = [chlorophyll(:,k)];
                fn = 2;
            else % night
                thomalla_correction(pair(k)).bbp_night = [thomalla_correction(pair(k)).bbp_night bbp(:,k)];
                thomalla_correction(pair(k)).chl_night = [thomalla_correction(pair(k)).chl_night chlorophyll(:,k)];
            end
        else % New day night pair
            % let's summarize the previous day/night pair
            try
                 m_bbp =  mean(thomalla_correction(pair(k)-1).bbp_night, 2, 'omitnan');
                 m_chl =  mean(thomalla_correction(pair(k)-1).chl_night, 2, 'omitnan');
                thomalla_correction(pair(k)-1).thomallacorrect = NaN(size(thomalla_correction(pair(k)-1).chl_night,1),size(thomalla_correction(pair(k)-1).bbp_day,2));
                for count = 1:size(thomalla_correction(pair(k)-1).bbp_day,2)
                    thomalla_correction(pair(k)-1).thomallacorrect(:,count) = m_chl ./ m_bbp .* thomalla_correction(pair(k)-1).bbp_day(:,count);
                end
            catch
                    miss = [miss, pair(k-2)];
                    disp(["Night Paired Backwards", pair(k)-1])
            end
            % now move onto the new data
            thomalla_correction(pair(k)).name = string(name);
            
            if angle(k) < 90 % day
                thomalla_correction(pair(k)).bbp_day = bbp(:,k);
                thomalla_correction(pair(k)).chl_day = chlorophyll(:,k);
            else % night
                thomalla_correction(pair(k)).bbp_night = bbp(:,k);
                thomalla_correction(pair(k)).chl_night = chlorophyll(:,k);
            end
        end
        
        if k == length(pair) % last time through the loop
            try
                m_bbp =  mean(thomalla_correction(pair(k)).bbp_night, 2, 'omitnan');
                m_chl =  mean(thomalla_correction(pair(k)).chl_night, 2, 'omitnan');
                for count = 1:size(thomalla_correction(pair(k)).bbp_day,2)
                    thomalla_correction(pair(k)).thomallacorrect(:,count) = m_chl ./ m_bbp .* thomalla_correction(pair(k)).bbp_day(:,count);
                end
            catch
                disp(["Manually pair day/night pair number", pair(k),k])
                miss = [miss, pair(k-1)];
            end
        end
    end 
    tt = 1;
    for k = 1:length(miss)
        j = miss(k);
        if size(thomalla_correction(j-1).chl_night,2) == 1
            thomalla_correction(j).thomallacorrect = thomalla_correction(j-1).chl_night./thomalla_correction(j-1).bbp_night.*thomalla_correction(j).bbp_day;
        else
            XXfact = mean(thomalla_correction(j-1).chl_night,2,'omitnan')./mean(thomalla_correction(j-1).bbp_night,2,'omitnan');
            if isempty(XXfact)
                thomalla_correction(j).thomallacorrect = NaN*thomalla_correction(j).bbp_day;
            else
                thomalla_correction(j).thomallacorrect = XXfact.*thomalla_correction(j).bbp_day;
        
            end
        end
    end
end