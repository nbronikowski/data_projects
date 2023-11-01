import numpy as np
import pandas as pd
import astral


def pgrid_columns(x, y2, var, yq):
	x = np.array(x)
	y = np.array(y2)
	var = np.array(var)
	temp_t = np.arange(len(var))
	id = ~np.isnan(var)

	# Interpolation
	y = np.interp(temp_t, temp_t[id], y2[id])
	y = y.reshape((-1, 1))

	xu = np.unique(x[~np.isnan(x)])

	temp_sort_var = np.empty((len(yq), len(xu)))
	temp_sort_var.fill(np.nan)

	for i, val in enumerate(xu):
		xidx = np.where(x == val)[0]

		temp_var = var[xidx]
		temp_y = y[xidx]
		idnan = np.where(~np.isnan(temp_var) & ~np.isnan(temp_y))[0]

		if len(idnan) > 3:
			temp_var = temp_var[idnan]
			temp_y = temp_y[idnan]

			loc = np.digitize(temp_y, yq) - 1  # similar to histcounts in MATLAB
			loc = loc.flatten()
			valid_loc = loc >= 0
			temp_y = temp_y.flatten()
			temp_y = temp_y[valid_loc]
			temp_var = temp_var[valid_loc]
			loc = loc[valid_loc]

			counts = np.bincount(loc)
			with np.errstate(divide='ignore', invalid='ignore'):  # Temporarily suppress warnings for the division operation
				y_mean = np.bincount(loc, temp_y) / counts
				var_mean = np.bincount(loc, temp_var) / counts

			y_mean[counts == 0] = np.nan  # Assign NaN where there are no data
			var_mean[counts == 0] = np.nan

			id = ~np.isnan(var_mean)
			if sum(id) > 2:
				temp_sort_var[:, i] = np.interp(yq, y_mean[id], var_mean[id], left=np.nan, right=np.nan)

	return temp_sort_var, xu

def OLDpgrid_columns(x, y2, var, yq):
	x = np.array(x)
	y = np.array(y2)
	var = np.array(var)
	temp_t = np.arange(len(var))
	id = ~np.isnan(var)
	
	# Interpolation
	y = np.interp(temp_t, temp_t[id], y2[id])
	y = y.reshape((-1, 1))

	xu = np.unique(x[~np.isnan(x)])

	temp_sort_var = np.empty((len(yq), len(xu)))
	temp_sort_var.fill(np.nan)

	for i, val in enumerate(xu):
		xidx = np.where(x == val)[0]

		temp_var = var[xidx]
		temp_y = y[xidx]
		idnan = np.where(~np.isnan(temp_var) & ~np.isnan(temp_y))[0]

		if len(idnan) > 3:
			temp_var = temp_var[idnan]
			temp_y = temp_y[idnan]
			
			loc = np.digitize(temp_y, yq) - 1  # similar to histcounts in MATLAB
			loc = loc.flatten()
			valid_loc = loc >= 0
			temp_y = temp_y.flatten()
			temp_y = temp_y[valid_loc]
			temp_var = temp_var[valid_loc]
			loc = loc[valid_loc]

			y_mean = np.bincount(loc, temp_y) / np.bincount(loc)
			var_mean = np.bincount(loc, temp_var) / np.bincount(loc)

			id = ~np.isnan(var_mean)
			if sum(id) > 2:
				temp_sort_var[:, i] = np.interp(yq, y_mean[id], var_mean[id], left=np.nan, right=np.nan)
				
	return temp_sort_var, xu


def delete_almost_empty_columns(dat, depthg):
	"""
	Removes columns from `dat` where the span of non-NaN data in depthg is less than 50.

	Parameters:
	- dat (numpy.ndarray): 2D array from which columns will be removed
	- depthg (numpy.ndarray): 1D array representing depth values

	Returns:
	- subsetData (numpy.ndarray): Data after removing columns
	- cols_to_keep (numpy.ndarray): Boolean array indicating which columns were kept
	"""
	
	N = dat.shape[1]  # Number of columns
	cols_to_keep = np.zeros(N, dtype=bool)  # Initialize a boolean array to identify columns to keep

	# Loop through each column to find the data span
	for t in range(N):
		non_nan_idx = np.where(~np.isnan(dat[:, t]))[0]  # Find non-NaN indices

		# Check if there are any non-NaN entries
		if len(non_nan_idx) > 0:
			data_span = depthg[non_nan_idx].max() - depthg[non_nan_idx].min() + 1  # Compute data span

			# Check if data span is at least 50
			if data_span >= 50:
				cols_to_keep[t] = True

	# Keep only columns with sufficient data span
	subset_data = dat[:, cols_to_keep]

	return subset_data, cols_to_keep

def interp1gap(x, y, xi, gap_thresh):
	"""
	Interpolate data with a gap threshold.
	
	x : array-like
		Array of x values.
	y : array-like
		Array of y values. Contains NaNs where gaps are present.
	xi : array-like
		Array of x values where interpolated data is needed.
	gap_thresh : float
		Gap threshold for interpolation. If a gap is larger than this value, it's not interpolated over.
	"""
	
	# Identify gaps
	gap_starts = np.where(np.isnan(y[:-1]) & ~np.isnan(y[1:]))[0]
	gap_ends = np.where(~np.isnan(y[:-1]) & np.isnan(y[1:]))[0]
	
	# Adjust the y data to introduce NaNs where gaps are larger than the threshold
	for start, end in zip(gap_starts, gap_ends):
		if (x[end] - x[start]) > gap_thresh:
			y[start:end] = np.nan
			
	# Now interpolate over the modified y data
	valid_mask = ~np.isnan(y)
	yi = np.interp(xi, x[valid_mask], y[valid_mask])
	
	return yi
	
	
def quenching_correction(fluorescence, backscatter, lat, lon, photic_layer,
						 night_day_group=True,
						 surface_layer=7,
						 rolling_window=3):

	"""
	Calculates the quenching depth and performs the quenching correction
	based on backscatter. The compulsory inputs must all be

	INPUT:
		fluorescence - pandas.DataFrame(index=depth, columns=surface_time, dtype=float) despiked
		backscatter  - pandas.DataFrame(index=depth, columns=surface_time, dtype=float) despiked
		lat          - pandas.DataFrame(index=depth, columns=surface_time, dtype=float)
		lon          - pandas.DataFrame(index=depth, columns=surface_time, dtype=float)
		photic_layer - pandas.DataFrame(index=depth, columns=surface_time, dtype=bool)
					   1% surface PAR True/False mask.
		night_day_group - True: quenching corrected with preceding night
						  False: quenching corrected with following night
		rolling_window  - We use a rolling window to find the quenching depth
						  the data may be spikey, this smooths it a little.
		surface_layer   - The surface fluorescence data may be quite noisy/spikey
						  hence we leave out the surface layer (meters)

	OUTPUT:
		corrected fluorescence
		quenching layer - boolean mask of quenching depth
		number of profiles per night used to correct quenching.

	METHOD:
		Correct for difference between night and daytime fluorescence.

		QUENCHING DEPTH
		===============
		The default setting is for the preceding night to be used to
		correct the following day's quenching (`night_day_group=True`).
		This can be changed so that the following night is used to
		correct the preceding day. The quenching depth is then found
		from the differnece between the night and daytime fluorescence.
		We use the steepest gradient of the {5 minimum differences and
		the points the differnece changes sign (+ve/-ve)}.

		BACKSCATTER / CHLOROPHYLL RATIO
		===============================
		1. Get the ratio between quenching depth and fluorescence
		2. Find the mean nighttime ratio for each night
		3. Get the ratio between nighttime and daytime quenching
		4. Apply the day/night ratio to the fluorescence
		5. If the corrected value is less than raw return to raw
	"""

	def sunset_sunrise(time, lat, lon):
		"""
		Uses the Astral package to find sunset and sunrise times.
		The times are returned rather than day or night indicies.
		More flexible for quenching corrections.
		"""

		ast = astral.Astral()

		df = pd.DataFrame.from_items([
			('time', time),
			('lat', lat),
			('lon', lon),
		])
		# set days as index
		df = df.set_index(df.time.values.astype('datetime64[D]'))

		# groupby days and find sunrise for unique days
		grp = df.groupby(df.index).mean()
		date = grp.index.to_pydatetime()

		grp['sunrise'] = list(map(ast.sunrise_utc, date, df.lat, df.lon))
		grp['sunset'] = list(map(ast.sunset_utc, date, df.lat, df.lon))

		# reindex days to original dataframe as night
		df[['sunrise', 'sunset']] = grp[['sunrise', 'sunset']].reindex(df.index)

		# set time as index again
		df = df.set_index('time', drop=False)
		cols = ['time', 'sunset', 'sunrise']
		return df[cols]

	def quench_nmin_grad(diff_ser, window, surface_layer):
		"""
		Quenching depth for a day/night fluorescence difference

		INPUT:   pandas.Series indexed by depth
				 window [4] is a rolling window size to remove spikes
				 skip_n_meters [5] skips the top layer that is often 0
		OUPUT:   estimated quenching depth as a float or int
				 note that this can be NaN if no fluorescence measurement
				 OR if the average difference is less than 0
		"""

		# When the average is NaN or less than 0 don't give a depth
		# Average difference of < 0 is an artefact
		if not (diff_ser.loc[:surface_layer].mean() > 0):
			return np.NaN

		# The rolling window removes spikes creating fake shallow QD
		x = diff_ser.rolling(window, center=True).mean()
		# We also skip the first N meters as fluorescence is often 0
		x_abs = x.loc[surface_layer:].abs()

		# 5 smallest absolute differences included
		x_small = x_abs.nsmallest(5)
		# Points that cross the 0 difference and make nans False
		sign_change = (x.loc[surface_layer:] > 0).astype(int).diff(1).abs()
		sign_change.iloc[[0, -1]] = False
		x_sign_change = x_abs[sign_change]
		# subset of x to run gradient on
		x_subs = pd.concat([x_small, x_sign_change])

		# find the steepest gradient from largest difference
		x_ref = x_subs - x.loc[:surface_layer].max()
		x_ref = x_ref[x_ref.notnull()]
		x_grad = (x_ref / x_ref.index.values).astype(float)
		# index of the largest negative gradient
		x_grad_min = x_grad.idxmin()

		return x_grad_min

	# ######################## #
	# DAY / NIGHT TIME BATCHES #
	# ######################## #
	# get the coordinates of the top 20 meters of the dives (surface)
	surf_lat = lat.loc[:20].mean()
	surf_lon = lon.loc[:20].mean()
	surf_time = fluorescence.columns.values

	# get the sunrise sunset times
	sun = sunset_sunrise(surf_time, surf_lat, surf_lon).astype('datetime64[ns]')
	# calculate day night times
	day = (sun.time > sun.sunrise) & (sun.time < sun.sunset)

	# creating quenching correction batches, where a batch is a
	# night and the following day
	if type(night_day_group) is not bool:
		raise TypeError("`night_day_group` must be boolean.")
	batch = (day.astype(int).diff().abs().cumsum() + night_day_group) // 2
	batch[0] = 0

	# Group the fluorescence by daytime and quenching batch
	grouped = fluorescence.groupby([day.values, batch.values], axis=1)
	fluo_night_median = grouped.median()[False]  # get the night values
	dives_per_night = grouped.count()[False].iloc[0]

	# Calculate the nighttime fluorescence and extrapolate this to day
	# so that the difference between night and day can be calculated
	fluo_night = fluo_night_median.reindex(columns=batch.values)
	fluo_night.columns = fluorescence.columns

	# #################################################### #
	# QUENCHING DEPTH BASED ON GRADIENT OF MIN DIFFERENCES #
	# #################################################### #
	# find the depth at which mean-nighttime and daytime fluorescence cross
	diff = (fluo_night - fluorescence).where(photic_layer)
	quench_depth = diff.apply(quench_nmin_grad, args=(rolling_window, surface_layer))

	quench_layer = diff.copy()
	idx = quench_layer.index.values
	quench_layer = quench_layer.apply(lambda s: idx < quench_depth[s.name])

	# #################################################################### #
	# QUENCHING CORRECTION FROM BACKSCATTER / FLUORESCENCE NIGHTTIME RATIO #
	# #################################################################### #
	# find mean fluorescence to backscatter ratio for nigttime
	flr_bb_ratio = backscatter / fluorescence
	flr_bb_grps = flr_bb_ratio.groupby([day.values, batch.values], axis=1)
	flr_bb_night = flr_bb_grps.mean()[False]
	flr_bb_night = flr_bb_night.reindex(columns=batch.values)
	flr_bb_night.columns = fluorescence.columns

	# quenching ratio for nighttime
	quench_ratio = (flr_bb_night * fluorescence / backscatter)
	quench_ratio = quench_ratio.where(quench_layer)

	# apply the quenching ratio to the fluorescence
	quench_correction = fluorescence / quench_ratio
	mask = quench_correction.notnull().values
	fluorescence_corrected = fluorescence.copy()
	fluorescence_corrected.values[mask] = quench_correction.values[mask]

	# if corrected fluorescence is lower than raw, return to raw
	mask = (fluorescence_corrected < fluorescence).values
	fluorescence_corrected.values[mask] = fluorescence.values[mask]

	return fluorescence_corrected, quench_layer, dives_per_night