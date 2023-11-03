function save_figure(h_fig,name,dims,extension,resolution)
    set(h_fig,'units','inches','position',[0 0 dims(1) dims(2)]);
%     set(h_fig,'color',[.89 .9 .5]);
    
    gcfcol=get(h_fig,'color');
    if gcfcol == [1 1 1]*0.94
        set(h_fig,'color','w');
    end
    set(h_fig,'InvertHardcopy','off');
% 
    if isempty(resolution)
        resolution = '-r300';
    else 
        resolution = ['-r',num2str(resolution)];
    end
    
    if strcmp(extension,'.png')
        print([name,'.png'],resolution,'-dpng');
    else if strcmp(extension,'.eps')
            print([name,'.eps'],resolution,'-depsc');
        end
    end
end