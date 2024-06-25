function VAR_FILT  = filtCorrectedProfiles(VAR,b,a,prof_idx,prof_dir)
    uidx = unique(prof_idx(~isnan(prof_idx)));
    VAR_FILT=VAR*NaN;
    for i = 1:length(uidx)
        idx1 = find(prof_idx == uidx(i) & prof_dir  == 1);
        idx2 = find(prof_idx == uidx(i)+1 & prof_dir == -1);
        if ~isempty(idx1) && ~isempty(idx2)
            idx = [idx1;idx2];
            xo = VAR(idx); % downcast

            % filter to remove deconv. effects
            VAR_FILT(idx) = safe_filtfilt(b,a,xo);
        end
    end
end

function filtered_data = safe_filtfilt(b, a, data)
    % Find the NaNs
    nanIdx = isnan(data);
    nonNaNData = data(~nanIdx);

    pad_size = 0; % Define a suitable pad size based on your data
    nonNaNData = [flipud(nonNaNData(1:pad_size)); nonNaNData; flipud(nonNaNData(end-pad_size+1:end))];
        
    % Apply filtfilt
    w = window(@rectwin,length(nonNaNData));
    nonNaNData = nonNaNData .* w;
    nonNaNFiltered = filtfilt(b, a, nonNaNData);
    
    % Remove the padding
    nonNaNFiltered = nonNaNFiltered(pad_size+1:end-pad_size);

    % Restore the NaNs
    filtered_data = nan(size(data));
    filtered_data(~nanIdx) = nonNaNFiltered;
end
