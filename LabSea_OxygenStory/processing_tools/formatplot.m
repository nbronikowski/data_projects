function formatplot()
%     set(gca,'DataAspectRatioMode','manual');
%     set(gca,'Xcolor','k','Ycolor','k','ZColor','k')
    set(gca,'LineWidth',0.5);
    box on; 
%     grid on
    set(gca,'XMinorTick','on','YMinorTick','on','ZMinorTick','on')
    set(gca,'TickLength',[.015 .015])
    set(gca,'Clipping','on')
%     set(gca,'Color',[200 170 0]./255)
    set(gca,'YColor','k');
    set(gca,'ClippingStyle','3dbox')
    set(gca,'TickDir','in','FontSize',11,'FontAngle','normal')
    set(gca,'BoxStyle','back');
    set(gca,'Layer','top');
end


