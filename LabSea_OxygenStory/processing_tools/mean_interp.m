function temp_intp=mean_interp(t,v,ti,extrap_flag)
    if nargin<4
        extrap_flag = 0;
    end
    temp_intp = NaN*ti;
    t = t(:); v = v(:); ti = ti(:);
    idx = ~isnan(v);
    temp_var  = v(idx);
    temp_x    = t(idx);
    [~,~,loc]=histcounts(temp_x,ti);
    temp_x(loc==0)=[]; temp_var(loc==0)=[]; loc(loc==0)=[];           
    x_mean   = accumarray(loc(:),temp_x(:))  ./accumarray(loc(:),1);
    var_mean = accumarray(loc(:),temp_var(:))./accumarray(loc(:),1);
    idnan = find(~isnan(var_mean));
    

%     id1 = loc(1);
%     id2 = loc(end);
%     
    if extrap_flag>0.5 && length(idnan)>2
        %     [~,id1]=nanmin(abs(ti-nanmin(x_mean(idnan))));
        %     [~,id2]=nanmin(abs(ti-nanmax(x_mean(idnan))));
        temp_intp = interp1(x_mean(idnan),var_mean(idnan),ti,'linear','extrap');
    else,if length(idnan)>2
            temp_intp = interp1(x_mean(idnan),var_mean(idnan),ti,'linear');
        else
            temp_intp = ti*NaN;
        end
    end
end