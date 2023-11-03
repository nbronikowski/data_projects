function RMSprof  = calcUpDownRMS(VAR,time,pres,prof_idx,prof_dir)
    p_bins = [0:2.5:ceil(max(pres,[],'omitnan'))]';
    p_bins_mid = p_bins+mean(diff(p_bins))/2;
    p_grid = [0:1:ceil(max(pres,[],'omitnan'))]';
    [b, a] = butter(1, 0.05);  % Geomar used 3rd order Butterworth, 0.05 normalized cutoff frequency
    x1avg = NaN(size(p_bins_mid));
    x2avg = NaN(size(p_bins_mid));
    uidx = unique(prof_idx(~isnan(prof_idx)));
    RMSprof = NaN(size(uidx));
    for i = 1:length(uidx)
        idx1 = find(prof_idx == uidx(i) & prof_dir  == 1);
        idx2 = find(prof_idx == uidx(i)+1 & prof_dir == -1);
        if ~isempty(idx1) && ~isempty(idx2)
            x1raw = VAR(idx1); % downcast
            p1raw = pres(idx1);
            x2raw = VAR(idx2); % upcast
            p2raw = pres(idx2);
            % average to be of p_bins_mod
            for j = 1:length(p_bins)-1
                bidx1 = find(p1raw>=p_bins(j) & p1raw<=p_bins(j+1));
                bidx2 = find(p2raw>=p_bins(j) & p2raw<=p_bins(j+1));
                if length(bidx1(~isnan(bidx1)))>1 && length(bidx2(~isnan(bidx2)))>1
                    
                    x1avg(j) = mean(x1raw(bidx1),'omitnan');
                    x2avg(j) = mean(x2raw(bidx2),'omitnan');
                else
                    x1avg(j) = NaN;
                    x2avg(j) = NaN;
                end
            end
            % interpolate var to 1 dbar grid
            x1intp = interp1(p_bins_mid,x1avg,p_grid);
            x2intp = interp1(p_bins_mid,x2avg,p_grid);
            
            if length(x1intp(~isnan(x1intp)))>3 && length(x2intp(~isnan(x2intp)))>3
                % filter to remove deconv. effects
                x1filt = safe_filtfilt(b,a,x1intp);
                x2filt = safe_filtfilt(b,a,x2intp);
                %  figure(1)
                %  plot(p_grid,x1intp,':r'); hold on
                %  plot(p_grid,x2intp,':b')
                %  plot(p_grid,x1filt,'r');
                %  plot(p_grid,x2filt,'b');

                % compute RMS
                RMSprof(i) = sqrt(mean((x1intp-x2intp).^2, 'omitnan'));
            else
                RMSprof(i) = NaN;
            end
        end
    end
end

