figure()
[X,~] = meshgrid(ds_subset.time,adat.range_cells);
pcolor(X,-(ds_subset.TrueDepthBeam1+ds_subset.depth),ds_subset.VelocityBeam1); 
shading flat; cb=colorbar; ylabel(cb,'UV Vel (m/s)')
title('Beam1 Velocity over time and depth for cells')
ylabel('depth (m)')
formatplot
datetick('x','dd.mm HH:MM','keeplimits')