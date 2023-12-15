function [var_intpy,var_avg] = gridcolumns(var,y,y_qc,yq,y_scale)
% Author Nicolai Bronikowski (MUN)
% nvob37@mun.ca

    % Force column
    yq = yq(:);

    y(y_qc>1)=NaN;
    
    % Reduce y to binned size vectors of same size
    y_scaled = scale_var(y,y_scale);
    
    % Define query lengths
    [~,Nx]=size(var);
    Ny = length(yq);

    % Allocate NaNs
    var_avg = NaN(Ny,Nx); 
    var_intpy = var_avg; 
    
    %%% bin average along bin grid defined above (2D lookup)
    for i = 1:Nx
        flags = find(~isnan(var(:,i)));
        if length(flags)>3
            var_avg(:,i)=mean_interp(y_scaled(:,i),var(:,i),yq);
        else
            var_avg(:,i)=yq*NaN;
        end
    end

    %%%  interpolate along first dimension
    for j = 1:Nx
        flags = find(~isnan(var_avg(:,j)));
        if length(flags)>10
            min_y = min(yq(flags))-y_scale*0;
            max_y = max(yq(flags))+y_scale*0;
            [~,id1] = min(abs(yq-min_y));
            [~,id2] = min(abs(yq-max_y));
            var_intpy(id1:id2,j) = interp1(yq(flags),var_avg(flags,j),yq(id1:id2),'linear');
        end
    end
end

