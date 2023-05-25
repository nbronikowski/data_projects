import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import contextily as ctx
from pyproj import CRS, Transformer
from matplotlib.colors import Normalize, ListedColormap
from shapely.geometry import Point

def set_plot_limits(ax, lon_min, lon_max, lat_min, lat_max):
	# Transformer for converting between WGS84 and Web Mercator
	transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857")
	# Transform input coordinates to Web Mercator
	x_min, y_min = transformer.transform(lat_min, lon_min)
	x_max, y_max = transformer.transform(lat_max, lon_max)
	ax.set_xlim(x_min, x_max)
	ax.set_ylim(y_min, y_max)

def findControlPoint(gdf_mercator, element, opts):
	# Create a GeoDataFrame with control point
	control_point = gpd.GeoDataFrame(geometry=[Point(opts['control_lon'], opts['control_lat'])], crs='EPSG:4326').to_crs(epsg=3857)
	# Calculate distance to control point
	gdf_mercator['distance_to_control'] = gdf_mercator.geometry.distance(control_point.geometry[0])

	# Subset to points within max distance
	least_contaminated = gdf_mercator.loc[gdf_mercator['distance_to_control'] <= opts['max_distance']]
	control_x, control_y, control_concentration = (
		least_contaminated['x'].mean(),
		least_contaminated['y'].mean(),
		least_contaminated[element].mean()
	)
	return control_x, control_y, control_concentration
	
def process_file(file_name, element, opts):
	df = pd.read_csv(file_name, na_values='< LOD')
	#df = df[(df['y'] <= opts['ignore_lat_max']) & (df['y'] >= opts['ignore_lat_min']) & (df['x'] <= opts['ignore_lon_max']) & (df['x'] >= opts['ignore_lon_min'])]
	
	gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x'], df['y']), crs='EPSG:4326')
	gdf_mercator = gdf.to_crs(epsg=3857) # Convert coordinate reference system
	gdf_mercator = gdf_mercator.dropna(subset=[element])
	control_x, control_y, control_concentration  = findControlPoint(gdf_mercator,element,opts)
	control_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([control_x],
															[control_y]), crs='EPSG:4326').to_crs(epsg=3857)
	element_deltaX = gdf_mercator[element]-control_concentration
	return gdf_mercator, control_gdf, element_deltaX

def plot_contamination_map(file_name, element, opts):
	gdf_mercator, control_gdf, element_deltaX = process_file(file_name, element, opts)
	N = 10
	cmap = ListedColormap([plt.cm.viridis(i / N) for i in range(N)])
	fig, ax = plt.subplots(figsize=(10, 10))

	# Normalize the element_deltaX values for colormap
	element_deltaX_normalized = Normalize()(element_deltaX)
	gdf_mercator['color'] = [cmap(x) for x in Normalize()(element_deltaX)]

	# Plot the points with color based on relative contamination
	#ax.scatter(gdf_mercator['x'], gdf_mercator['y'], c=element_deltaX, cmap=cmap, s=50)
	gdf_mercator.plot(column='color', ax=ax, markersize=50)
	control_gdf.plot(ax=ax, color='red', alpha=1, edgecolor='none', markersize=50)

	# Add a colorbar
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=element_deltaX.min(), vmax=element_deltaX.max()))
	sm.set_array([])
	colorbar = fig.colorbar(sm, ax=ax, orientation='vertical', shrink=0.4)
	colorbar.set_label(f'Difference in {element} Concentration Relative to Control Point (ppm)')

	# Set plot limits
	# Calculate longitude and latitude for tick locations, set labels, etc.
	set_plot_limits(ax, opts['map_lon_min'], opts['map_lon_max'], opts['map_lat_min'], opts['map_lat_max'])
	xmin, xmax = ax.get_xlim()
	ymin, ymax = ax.get_ylim()
	xticks = np.linspace(xmin, xmax, num=5)
	yticks = np.linspace(ymin, ymax, num=5)
	transformer = Transformer.from_crs("epsg:3857", "epsg:4326")
	xticklabels = [round(transformer.transform(x, ymin)[1], 4) for x in xticks]
	yticklabels = [round(transformer.transform(xmin, y)[0], 4) for y in yticks]
	ax.set_xticks(xticks)
	ax.set_yticks(yticks)
	ax.set_xticklabels(xticklabels)
	ax.set_yticklabels(yticklabels)
	ax.set_xlabel('Longitude')
	ax.set_ylabel('Latitude')
	
	plt.title('Soil Sample Locations and Contamination Relative to Control Point')
	ctx.add_basemap(ax, source=ctx.providers.OpenTopoMap)
	# add_inset_map(ax, control_gdf, lon_min, lat_min, lon_max, lat_max) # Add inset map
	plt.savefig(f"./plots/contamination_map_{element}.png", dpi=300)

'''
def add_inset_map(ax, gdf, lon1, lat1, lon2, lat2):
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	transformer_inset = Transformer.from_crs('EPSG:4326', 'EPSG:3857')
	iax = inset_axes(ax, width="30%", height="33%", loc=1)
	inset_xmin, inset_ymin = transformer_inset.transform(lat2, lon2)
	inset_xmax, inset_ymax = transformer_inset.transform(lat1, lon1)
	iax.set_xlim(min(inset_xmin, inset_xmax), max(inset_xmin, inset_xmax))
	iax.set_ylim(min(inset_ymin, inset_ymax), max(inset_ymin, inset_ymax))
	ctx.add_basemap(iax, source=ctx.providers.OpenTopoMap, zoom=12, attribution_size=0)
	main_xmin, main_xmax =  ax.get_xlim()
	main_ymin, main_ymax =  ax.get_ylim()
	iax.plot([main_xmin, main_xmin, main_xmax, main_xmax, main_xmin],
			 [main_ymin, main_ymax, main_ymax, main_ymin, main_ymin],
			 color='red', linewidth=1.5, linestyle='-')
	iax.set_xticks([]), iax.set_yticks([])

def findControlPoint(gdf_mercator, element, opts):
	from scipy.cluster.hierarchy import dendrogram, linkage
	from scipy.cluster.hierarchy import fcluster
	# First, lets drop rows that have NaN values in the specific element column
	gdf_mercator = gdf_mercator.dropna(subset=[element])
	# normalize the element values to be on a similar scale as the coordinates
	gdf_mercator['norm_element'] = (gdf_mercator[element] - gdf_mercator[element].min()) / (gdf_mercator[element].max() - gdf_mercator[element].min())
	# generate the linkage matrix
	Z = linkage(gdf_mercator[['x', 'y', 'norm_element']], 'ward')
	# set a distance threshold to form the flat clusters
	max_d = opts['max_distance'] # adjust this value based on the dendrogram
	clusters = fcluster(Z, max_d, criterion='distance')
	gdf_mercator['cluster'] = clusters
	# find the average concentration for each cluster
	cluster_averages = gdf_mercator.groupby('cluster')[element].mean()
	lowest_cluster = cluster_averages.idxmin()
	least_contaminated = gdf_mercator[gdf_mercator['cluster'] == lowest_cluster]
	control_x, control_y, control_concentration = (
		least_contaminated['x'].mean(),
		least_contaminated['y'].mean(),
		least_contaminated[element].mean()
	)
	return control_x, control_y, control_concentration, least_contaminated

def process_file(file_name, element, opts):
	df = pd.read_csv(file_name, na_values='< LOD')
	control_x, control_y, control_concentration,_ = findControlPoint(gdf_mercator,element,opts)
	control_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([control_x],
															[control_y]), crs='EPSG:4326').to_crs(epsg=3857)
	df = df[(df['y'] <= opts['ignore_lat_max']) & (df['y'] >= opts['ignore_lat_min']) & (df['x'] <= opts['ignore_lon_max']) & (df['x'] >= opts['ignore_lon_min'])]
	gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x'], df['y']), crs='EPSG:4326')
	gdf_mercator = gdf.to_crs(epsg=3857) # Convert coordinate reference system

	# Compute differences in coordinates and concentrations
	dx = gdf_mercator.geometry.x - control_gdf.geometry.x.iloc[0]
	dy = gdf_mercator.geometry.y - control_gdf.geometry.y.iloc[0]
	
	# Compute distance between points and control points
	distance = np.sqrt(dx**2 + dy**2)
	d_element = df[element] - control_concentration  # Change element here
	
	# Compute the change in contamination per unit distance
	delta_element_deltaX = np.abs(d_element / distance)
	# Normalize for color mapping
	delta_element_deltaX_normalized = (
		(delta_element_deltaX - delta_element_deltaX.min()) /
		(delta_element_deltaX.max() - delta_element_deltaX.min())
	)

	dx_normalized = dx / distance
	dy_normalized = dy / distance
	
	# Multiply the normalized direction vector by the absolute gradient magnitude to create a vector
	dx_gradient = dx_normalized * delta_element_deltaX
	dy_gradient = dy_normalized * delta_element_deltaX
	return gdf_mercator, control_gdf, dx_gradient, dy_gradient, delta_element_deltaX, delta_element_deltaX_normalized

def plot_contamination_map(file_name, element, opts):
	gdf_mercator,control_gdf,dx_grad,dy_grad,element_deltaX,element_deltaX_normalized=process_file(file_name,element, opts)
	N = 10
	cmap = ListedColormap([plt.cm.viridis(i/N) for i in range(N)])
	fig, ax = plt.subplots(figsize=(10,10))
	
	#gdf_mercator.plot(ax=ax, alpha=1, edgecolor='none')
	control_gdf.plot(ax=ax,color='red' ,alpha=1, edgecolor='none', markersize=50)
	q = ax.quiver(gdf_mercator.geometry.x, gdf_mercator.geometry.y, dx_grad, dy_grad,color=cmap(element_deltaX_normalized), linewidths=1)
	
	# Add a colorbar
	#norm = Normalize(vmin=abs(element_deltaX).min(), vmax=abs(element_deltaX).max())
	norm = Normalize(vmin=element_deltaX.min(), vmax=element_deltaX.max())
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	sm.set_array([])
	colorbar = fig.colorbar(sm, ax=ax, orientation='vertical', shrink=0.4)
	colorbar.set_ticks(np.linspace(np.round(abs(element_deltaX).min(), 2),
							   np.round(abs(element_deltaX).max(), 2), num=10))
	colorbar.set_ticks(np.linspace(np.round(element_deltaX.min(), 2),
								   np.round(element_deltaX.max(), 2), num=10))
	colorbar.set_label(f'Change in {element} concentration (ppm/m)')  # Replaced 'Mn' with {element}

	# Set plot limits
	set_plot_limits(ax, opts['map_lon_min'],opts['map_lon_max'],opts['map_lat_min'],opts['map_lat_max'])
	# Calculate longitude and latitude for tick locations, set labels, etc.
	xmin, xmax = ax.get_xlim()
	ymin, ymax = ax.get_ylim()
	xticks = np.linspace(xmin, xmax, num=5)
	yticks = np.linspace(ymin, ymax, num=5)
	transformer = Transformer.from_crs("epsg:3857", "epsg:4326")
	xticklabels = [round(transformer.transform(x, ymin)[1], 4) for x in xticks]
	yticklabels = [round(transformer.transform(xmin, y)[0], 4) for y in yticks]
	ax.set_xticks(xticks)
	ax.set_yticks(yticks)
	ax.set_xticklabels(xticklabels)
	ax.set_yticklabels(yticklabels)
	ax.set_xlabel('Longitude')
	ax.set_ylabel('Latitude')
	plt.title('Soil Sample Locations and Gradients of Contamination Relative to Local Minimum')
	# Update basemap
	ctx.add_basemap(ax, source=ctx.providers.OpenTopoMap)
	#add_inset_map(ax, control_gdf, lon_min, lat_min, lon_max, lat_max) # Add inset map
	plt.savefig(f"./plots/contamination_map_{element}.png", dpi=300)
'''

