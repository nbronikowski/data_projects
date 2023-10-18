figure()
[X,~] = meshgrid(ds_subset.time,adat.range_cells);
pcolor(X,-(ds_subset.TrueDepth+ds_subset.depth),ds_subset.UVelocity); 
shading flat; cb=colorbar; ylabel(cb,'UV Vel (m/s)')
title('UV ENU Velocity over time and depth for cells')
ylabel('depth (m)')
formatplot
datetick('x','dd.mm HH:MM','keeplimits')