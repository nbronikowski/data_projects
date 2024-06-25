figure() 
subplot(311);
scatter(X,-Y,20,Z)
cb=colorbar;
xlabel('time (days)')
ylabel('depth (m)')
xlim([min(X) max(X)])
ylim([-max(Y) -min(Y)])

subplot(312);hold on;
imagescn(gridX,-gridY,Z_grid)
scatter(X,-Y,20,Z)
cb2=colorbar;
xlabel('time (days)')
ylabel('depth (m)')
xlim([min(X) max(X)])
ylim([-max(Y) -min(Y)])



[dZdx,dZdy]=gradient(Z_grid,gridX(1,:),gridY(:,1));
subplot(313);hold on;
imagescn(gridX,-gridY,dZdx)
% scatter(X,-Y,20,Z)
cb2=colorbar;
xlabel('time (days)')
ylabel('depth (m)')
xlim([min(X) max(X)])
ylim([-max(Y) -min(Y)])
