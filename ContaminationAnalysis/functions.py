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

def add_text_box(ax,control_gdf):
	# Create a transformer for the coordinate conversion
	transformer_mercator_to_latlon = Transformer.from_crs("epsg:3857", "epsg:4326")

	# Convert Mercator meters to degrees north and east
	lon, lat = transformer_mercator_to_latlon.transform(control_gdf.geometry.x.values[0], control_gdf.geometry.y.values[0])

	# Add text box for control point concentration and location
	control_point_text = f"Control Point (Blue):\nConcentration: {control_gdf['control_concentration'][0]:.1f} ppm\nLocation: ({lon:.2f}, {lat:.2f})"
	ax.text(0.02, 0.98, control_point_text, transform=ax.transAxes, fontsize=10, verticalalignment='top',
			bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
	

def process_file(file_name, element, opts):
	df = pd.read_csv(file_name, na_values='< LOD')
	gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x'], df['y']), crs='EPSG:4326')
	gdf = gdf.to_crs(epsg=3857) # Convert coordinate reference system
	gdf = gdf.dropna(subset=[element])
	
	if opts['control_point']:
		if 'control_concentration' not in opts:
			control_x, control_y, control_concentration = findControlPoint(gdf, element, opts)
			control_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([control_x], [control_y]), crs='EPSG:4326').to_crs(epsg=3857)
			control_gdf['control_concentration'] = control_concentration
		else:
			control_x, control_y, control_concentration = opts['control_lon'], opts['control_lat'], opts['control_concentration']
			control_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([control_x], [control_y]), crs='EPSG:4326').to_crs(epsg=3857)
			control_gdf['control_concentration'] = control_concentration
	else:
		control_x, control_y, control_concentration = 0, 0, 0
		control_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy([control_x], [control_y]), crs='EPSG:4326').to_crs(epsg=3857)
		control_gdf['control_concentration'] = control_concentration
	
	# now remove points in a gdf that is focussed on a specific area: gdf_filt
	df = df[
		(df['y'] <= opts['ignore_lat_max']) &
		(df['y'] >= opts['ignore_lat_min']) &
		(df['x'] <= opts['ignore_lon_max']) &
		(df['x'] >= opts['ignore_lon_min'])
	]
	gdf_filt = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df['x'], df['y']), crs='EPSG:4326')
	gdf_filt = gdf_filt.to_crs(epsg=3857) # Convert coordinate reference system
	gdf_filt = gdf_filt.dropna(subset=[element])
	gdf_filt['element_deltaX'] = gdf_filt[element]-control_concentration
	 
	return gdf, gdf_filt, control_gdf


def add_inset_map(ax, lon1, lat1, lon2, lat2, gdf, control_gdf):
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes

	transformer_inset = Transformer.from_crs('EPSG:4326', 'EPSG:3857')
	iax = inset_axes(ax, width="35%", height="35%", loc=1,
					 bbox_to_anchor=(0.2, 0.23, 0.8, 0.8), bbox_transform=ax.transAxes, borderpad=0.2)
	
	inset_xmin, inset_ymin = transformer_inset.transform(lat1, lon1)
	inset_xmax, inset_ymax = transformer_inset.transform(lat2, lon2)
	iax.set_xlim(min(inset_xmin, inset_xmax), max(inset_xmin, inset_xmax))
	iax.set_ylim(min(inset_ymin, inset_ymax), max(inset_ymin, inset_ymax))

	iax.scatter(gdf.geometry.x, gdf.geometry.y, color='black', s=20)
	iax.scatter(control_gdf.geometry.x, control_gdf.geometry.y, color='blue', s=20)

	ctx.add_basemap(iax, source=ctx.providers.OpenTopoMap, zoom=11, attribution_size=0)
	main_xmin, main_xmax =  ax.get_xlim()
	main_ymin, main_ymax =  ax.get_ylim()
	iax.plot([main_xmin, main_xmin, main_xmax, main_xmax, main_xmin],
			 [main_ymin, main_ymax, main_ymax, main_ymin, main_ymin],
			 color='red', linewidth=1.5, linestyle='-')
	iax.set_xticks([]), iax.set_yticks([])
	
def add_inset_map2(ax, lon1, lat1, lon2, lat2, gdf, gdf_filt):
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	cmap = ListedColormap([plt.cm.viridis(i / 10) for i in range(10)])
	transformer_inset = Transformer.from_crs('EPSG:4326', 'EPSG:3857')
	iax = inset_axes(ax, width="35%", height="35%", loc=1,
					 bbox_to_anchor=(0.01, 0.01, 0.8, 0.8), bbox_transform=ax.transAxes, borderpad=0.2)
	
	inset_xmin, inset_ymin = transformer_inset.transform(lat1, lon1)
	inset_xmax, inset_ymax = transformer_inset.transform(lat2, lon2)
	iax.set_xlim(min(inset_xmin, inset_xmax), max(inset_xmin, inset_xmax))
	iax.set_ylim(min(inset_ymin, inset_ymax), max(inset_ymin, inset_ymax))
	
	gdf_filt["normalized_contamination"] = Normalize()(gdf_filt['element_deltaX'])
	# Sort GeoDataFrame based on normalized contamination
	gdf_filt = gdf_filt.sort_values(by="normalized_contamination")
	scatter = iax.scatter(gdf_filt.geometry.x, gdf_filt.geometry.y, c=gdf_filt["normalized_contamination"], cmap=cmap, s=50)

	ctx.add_basemap(iax, source=ctx.providers.OpenTopoMap, zoom=11, attribution_size=0)
	iax.set_xticks([]), iax.set_yticks([])

def plot_contamination_map(file_name, element, opts):
	gdf, gdf_filt, control_gdf = process_file(file_name, element, opts)
	
	cmap = ListedColormap([plt.cm.viridis(i / 10) for i in range(10)])
	fig, ax = plt.subplots(figsize=(7, 6))
	fig.set_facecolor('white')
	# Normalize the element_deltaX values for colormap
	gdf_filt["normalized_contamination"] = Normalize()(gdf_filt['element_deltaX'])

	# Sort GeoDataFrame based on normalized contamination
	gdf_filt = gdf_filt.sort_values(by="normalized_contamination")

	# Plot the points with color based on relative contamination
	scatter = ax.scatter(gdf_filt.geometry.x, gdf_filt.geometry.y, c=gdf_filt["normalized_contamination"], cmap=cmap, s=50)

	# Set plot limits
	set_plot_limits(ax, opts['map_lon_min'], opts['map_lon_max'], opts['map_lat_min'], opts['map_lat_max'])

	# Plot the base map
	ctx.add_basemap(ax, source=ctx.providers.OpenTopoMap)

	# Calculate percentiles for the colorbar limits
	vmin, vmax = np.percentile(gdf_filt['element_deltaX'], [10, 90])
    
	# Add a colorbar
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
	sm.set_array([])
	colorbar = fig.colorbar(sm, ax=ax, orientation='vertical', shrink=0.4)
	colorbar.set_label(f'Difference in {element} Concentration Relative to Control Point (ppm)')
	
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

	# Add text box for control point concentration and location
	add_text_box(ax,control_gdf)
	
	plt.title(opts['plot_title'])
	
	add_inset_map(ax,-60.595,53.25, -60.30, 53.404,gdf,control_gdf) # Add inset map
	add_inset_map2(ax,-60.4062,53.285, -60.3873, 53.296,gdf,gdf_filt) # Add inset map

	plt.savefig(f"{opts['plot_dir']}contamination_map_{element}.png", dpi=300)
