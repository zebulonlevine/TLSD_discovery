from trackingCodeFunctions import *
import os
import sys
import re

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
from PIL import Image, ImageDraw

plt.rcParams['font.family'] = 'DejaVu Sans'

import pandas as pd
import numpy as np
import glob

import imageio
import pims
from tifffile import imwrite

import warnings
warnings.filterwarnings("ignore")


def find_outlier_filenames(posterior_df , mult_on_sd = 2):
    outliers = []
    grouped_density = [
        group["density"].to_numpy()
        for _, group in posterior_df.groupby("file")
    ]
    group_names = posterior_df["file"].drop_duplicates().tolist()
    outliers.append(find_outlier_positions(grouped_density  , group_names , mult_on_sd = mult_on_sd))
    outliers = [item for sublist in outliers for item in sublist]
    return outliers




## this does the statistical *work* underlying the mark_outliers function
## starting from the posterior populations... calculate the median of the
## posterior probabilty for each diffusion rate. for each XY position, calculate
## the KS-distance from the median distribution, and keep any XY position within mult_on_sd
## of that distribution
def find_outlier_positions(grouped_density , group_names , mult_on_sd = 2):

    # Assume `pdfs` is a list or array of 1D arrays (n_pdfs x n_points)
    pdfs = np.array(grouped_density)

    # Step 1: Compute the empirical CDFs
    cdfs = np.cumsum(pdfs, axis=1)
    cdfs = cdfs / cdfs[:, -1][:, np.newaxis]  # Normalize each to 1

    # Step 2: Compute the median CDF
    median_cdf = np.median(cdfs, axis=0)

    # Step 3: KS test between each PDF's CDF and the median CDF
    def ks_statistic(empirical, reference):
        return np.max(np.abs(empirical - reference))

    ks_stats = [ks_statistic(cdf, median_cdf) for cdf in cdfs]

    # Step 4: Flag outliers — e.g., using robust threshold
    threshold = np.median(ks_stats) + mult_on_sd * np.std(ks_stats) ## SD is weird statistic because distribution is one-sided
    outliers = [i for i, stat in enumerate(ks_stats) if stat > threshold]
    outlier_files = [group_names[i] for i in outliers]

#     print("Outlier indices:", outliers)
#     # Print the result
#     for file in outlier_files:
#         print(f"  - {file}:")
    return(outlier_files)

## this calculates the weighted average posterior diffusion rate, i.e. it marginalizes over the whol
## population of molecules, then scales by the population size to make it a PDF
## returns a plottable table of the posterior MLE diffusion rate, with a column to store the condition name
def generate_posterior_plotting_df(posterior_df):
    posterior_plot_df = posterior_df.copy()
    posterior_plot_df = posterior_plot_df.groupby("D")["density"].mean().reset_index()
    posterior_plot_df['density'] = posterior_plot_df['density'] / posterior_plot_df['density'].sum()
    posterior_plot_df['condition'] = np.unique([re.sub('_[0-9]*_traj.csv' , '' , f) for f in posterior_df['file']])[0]
    return(posterior_plot_df)




## take MLE_df and scaled posterior_plot_df and generate a histogram and bar graph
## The safe usage of this function is to pass only one condition. could likely be modified
## to accept many conditions, and plot many plots into a single PDF, but that hasn't been implemented yet
## lightgray lines correspond to file-wise PDFs, black line is mean PDF
def generate_plots(MLE_df , posterior_plot_df , all_posteriors , settings , showPlot = True , conditions = []):

    n = len(conditions)
    if n == 0:
        conditions = MLE_df['condition'].unique()

    bin_labels = settings['plot']['barGraphLabels']
    barGraphs = pd.DataFrame(index=bin_labels)

    #print(f"Working on conditions: {conditions}")

    for condition in conditions:
        sm_df = posterior_plot_df[posterior_plot_df['condition'] == condition]
        sm_MLE = MLE_df[MLE_df['condition'] == condition]
        sm_all_posterior = all_posteriors[all_posteriors['condition']==condition]

        normalized_df = sm_all_posterior.copy()
        normalized_df['density_norm'] = normalized_df.groupby('file')['density'].transform(lambda x: x / x.sum())

        # Create a new figure for each condition
        #fig, ax = plt.subplots(figsize=(8, 2.5))
        fig, axes = plt.subplots(
            1, 2,                   # 1 row, 2 columns
            figsize=(9, 2.5 * len(conditions)),        # total figure size = 8 + 1 width, 2.5 height
            gridspec_kw={'width_ratios': [8, 1]}  # relative widths
        )

        ax_hist, ax_bar = axes

        plot_xlim = settings['plot']['xlim']  # default (-3, 2)

        # MLE_D distribution: KDE (smooth) or histogram (raw bins)
        log_mle = np.log10(sm_MLE['MLE_D'].values)
        log_mle = log_mle[np.isfinite(log_mle)]

        if settings['plot'].get('smooth_mle', True):
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(log_mle, bw_method=0.1)
            x_kde = np.linspace(plot_xlim[0], plot_xlim[1], 500)
            y_kde = kde(x_kde)
            ax_hist.fill_between(x_kde, y_kde, alpha=0.45, color='#D28F8F')
            ax_hist.plot(x_kde, y_kde, color='#B05050', linewidth=1.2)
        else:
            ax_hist.hist(log_mle,
                         bins=np.linspace(plot_xlim[0], plot_xlim[1],
                                          settings['plot']['nbins']),
                         density=True, color="#D28F8F");

        for f, group in normalized_df.groupby('file'):
            ax_hist.plot(np.log10(group['D']), 50 * group['density_norm'], linewidth=1.5, color='darkgray' , alpha = 0.4)

        ax_hist.plot(np.log10(sm_df['D']) , 50 * sm_df['density'] , linewidth = 1.5 , color = 'black')

        ax_hist.set_xlim(plot_xlim)

        tick_positions = list(range(int(plot_xlim[0]), int(plot_xlim[1]) + 1))
        tick_labels = [f"{10**v:.4f}" if 10**v < 0.01 else f"{10**v:.2f}" if 10**v < 1 else f"{10**v:.0f}" for v in tick_positions]
        ax_hist.set_xticks(tick_positions, tick_labels);

        minor_ticks = []
        for base_exp in range(int(plot_xlim[0]), int(plot_xlim[1]) + 1):
            base = 10 ** base_exp
            for i in range(2, 10):
                tick = base * i
                if tick <= 10 ** plot_xlim[1]:
                    minor_ticks.append(tick)
        minor_ticks = np.log10(minor_ticks)
        ax_hist.set_xticks(minor_ticks, minor=True)
        ax_hist.tick_params(axis='x', which='minor', length=4, color='gray')

        ax_hist.set_ylim(settings['plot']['ylim'][0] , settings['plot']['ylim'][1])

        ax_hist.set_title(f"{condition} ; n_jumps = {settings['saspt']['splitsize']} ; {len(sm_MLE['MLE_D'])} particles", pad=20)
        ax_hist.set_xlabel("Value")
        ax_hist.set_ylabel("Probability Density")

        for spine in ['top', 'right']:
            ax_hist.spines[spine].set_visible(False)

        # Draw vertical lines at cutoff boundaries
        for brk in settings['plot']['barGraphBreaks'][1:-1]:
            ax_hist.axvline(np.log10(brk), color='gray', linestyle='--',
                            linewidth=0.8, alpha=0.6)

        counts, _ = np.histogram(sm_MLE['MLE_D'], bins=settings['plot']['barGraphBreaks']);
        barGraphs[condition] = counts

        barGraphs = barGraphs.T
        barGraphs = barGraphs.div(barGraphs.sum(axis=1), axis=0)
        barGraphs.plot(kind="bar", stacked=True, ax=ax_bar, edgecolor='black')
        #ax_bar.set_xlabel("Condition")
        ax_bar.set_ylabel("Fraction")
        if len(barGraphs) == 1:
            ax_bar.set_xticklabels([]) ## if there is only one condition, suppress title
        ax_bar.legend(title="Population", bbox_to_anchor=(1.05, 1), loc="upper left")

        plt.tight_layout()

        # Save to PDF with condition in filename
        fig.savefig(f"{settings['io']['plot_directory']}/{condition}_hist.pdf", dpi = 300 , format="pdf", bbox_inches='tight')

        if showPlot:
            plt.show()

        plt.close(fig)  # Close to avoid memory buildup
      
      



## QC movies
## Bleaching analysis

## from a data frame of trajectories given by all_traj, extract the trajectories matching
## the name ur_trajectory, and return a padded data frame to pass to movie-making code
def isolate_traj(all_traj , ur_trajectory , half_size = 15 , max_frame = 500):
    this_traj = all_traj[all_traj['ur_trajectory'] == ur_trajectory].copy()
        
    this_traj_center = int(np.mean(this_traj['y'])) , int(np.mean(this_traj['x']))
    timeBox = np.min(this_traj['frame']) - 5 , np.max(this_traj['frame']) + 5
    if timeBox[0] < 0:
        timeBox = 0 , timeBox[1]
    if timeBox[1] > 500:
        timeBox = timeBox[0] , 500
    #print(timeBox)

    spaceBox = {'xmin' : (this_traj_center[1] - half_size) , 
                          'xmax' : (this_traj_center[1] + half_size + 1) ,
                          'ymin' : (this_traj_center[0] - half_size) , 
                          'ymax' : (this_traj_center[0] + half_size + 1)}

    this_traj['og_frame'] = this_traj['frame'] ## store the global frame number where the trajectory started
    this_traj['rel_frame'] = this_traj['frame'] -  np.min(timeBox)
        
    first_frame = this_traj['frame'].min()
    last_frame = this_traj['frame'].max()
    # Generate the new frame numbers
    new_frames = list(range(first_frame - 5, first_frame)) + \
                 list(range(last_frame + 1, last_frame + 6))
    new_frames = [x for x in new_frames if x > 0 & x < max_frame]

    # Create a DataFrame of NaN rows with the same columns
    nan_rows = pd.DataFrame(np.nan, index=range(len(new_frames)), columns=this_traj.columns)
    traj_id = this_traj['ur_trajectory'].iloc[0]

    # Assign the new frame numbers
    nan_rows['frame'] = new_frames
    nan_rows['og_frame'] = new_frames
    #nan_rows['rel_frame'] = new_frames
    nan_rows['ur_trajectory'] = traj_id
    nan_rows['condition'] = this_traj['condition'].iloc[0]
    nan_rows['rel_frame'] = nan_rows['frame'] - np.min(timeBox) ## frame relative to the bounding box time interval start
    
    bounding_box = {'timeBox' : timeBox , 'spaceBox' : spaceBox}
    
    # Concatenate and sort
    this_traj = pd.concat([nan_rows, this_traj], ignore_index=True).sort_values('frame', ignore_index=True)
    return(this_traj , bounding_box)



## colorized version
def annotate_trajectory_building_tail(box, traj_df, spaceBox, upscale=8, spot_radius=8, lwd=4):
    """
    Annotate trajectory with subpixel spots and draw all prior jumps between adjacent timepoints.
    Each frame shows the trajectory up to that time point, with line and spot colors based on MLE_D.

    Color is directly mapped to MLE_D in the fixed range [0.01, 100], no normalization.

    Parameters:
    - box: 3D numpy array (T, Y, X)
    - traj_df: DataFrame with 'frame', 'x', 'y', 'ur_trajectory', 'MLE_D'
    - spaceBox: dict with 'xmin', 'xmax', 'ymin', 'ymax'
    - upscale: scale factor for subpixel rendering
    - spot_radius: radius of spot circle to draw
    - lwd: line width for trajectory segments

    Returns:
    - upsized_input: list of upscaled input images
    - annotated_frames: list of RGB arrays with annotations
    """

    annotated_frames = []
    upsized_input = []

    # Sort and group the detections
    traj_df = traj_df.sort_values(['ur_trajectory', 'frame'])
    grouped = traj_df.groupby('ur_trajectory')

    # Precompute upscaled trajectory lines
    trajectory_segments_by_frame = {t: [] for t in range(len(box))}

    # Colormap and fixed value range
    cmap = cm.viridis
    vmin, vmax = -2, 2

    # Helper: map MLE_D to colormap without normalization
    def mle_to_color(val):
        val_clipped = np.clip(np.log10(val), vmin, vmax)
        fraction = (val_clipped - vmin) / (vmax - vmin)
        return tuple(int(255*c) for c in cmap(fraction)[:3])

    # Precompute segments with colors
    for _, group in grouped:
        group = group.reset_index(drop=True)
        for i in range(1, len(group)):
            f_prev = group.loc[i-1, 'rel_frame']
            f_curr = group.loc[i, 'rel_frame']
            if f_curr == f_prev + 1:
                x1, y1 = (group.loc[i-1, 'x'] - spaceBox['xmin']) * upscale, (group.loc[i-1, 'y'] - spaceBox['ymin']) * upscale
                x2, y2 = (group.loc[i, 'x'] - spaceBox['xmin']) * upscale, (group.loc[i, 'y'] - spaceBox['ymin']) * upscale

                val = group.loc[i, 'MLE_D']
                if not np.isnan(val):
                    color = mle_to_color(val)
                    for t in range(f_curr, len(box)):
                        trajectory_segments_by_frame[t].append(((x1, y1), (x2, y2), color))

    # Render each frame
    for t in range(len(box)):
        height, width = box[t].shape  # dynamically get size
        frame = np.zeros((height, width, 3), dtype=np.uint8)
        
        img = Image.fromarray(frame, 'RGB')
        img_up = img.resize((img.width * upscale, img.height * upscale), resample=Image.BICUBIC)
        draw = ImageDraw.Draw(img_up)

        boximg = box[t]
        boximg = Image.fromarray(boximg)
        #boximg = boximg.convert("RGB")
        boximg_up = boximg.resize((boximg.width * upscale, boximg.height * upscale), resample=Image.BICUBIC)

        # Draw all previous displacement segments
        for x1y1, x2y2, color in trajectory_segments_by_frame[t]:
            x1, y1 = x1y1
            x2, y2 = x2y2
            if np.isnan([x1, x2, y1, y2]).any():
                continue
            draw.line([(x1, y1), (x2, y2)], fill=color, width=lwd)

        # Draw current spots colored by MLE_D
        current_detections = traj_df[traj_df['rel_frame'] == t]
        for _, row in current_detections.iterrows():
            x = (row['x'] - spaceBox['xmin']) * upscale
            y = (row['y'] - spaceBox['ymin']) * upscale
            bbox = [x - spot_radius, y - spot_radius, x + spot_radius, y + spot_radius]

            if not np.isnan(row['MLE_D']):
                ellipse_color = mle_to_color(row['MLE_D'])
                draw.ellipse(bbox, outline=ellipse_color, width=lwd)

        annotated_frames.append(np.array(img_up))
        upsized_input.append(np.array(boximg_up))

    return upsized_input, annotated_frames




def save_overlay_movie(box_resized, box_annotated, output_path='overlay_output.mp4', fps=20):
    frames_rgb = []

    ## ensure the image is 8-bit and invert BW
    box_resized = box_resized - np.min(box_resized)
    box_resized = (box_resized / np.max(box_resized) * 255).astype(np.uint8)
    box_resized = [255 - frame for frame in box_resized]
    
    #box_resized = box_resized - np.min(box_resized)
    #box_resized = (box_resized / np.max(box_resized) * 255).astype(np.uint8)
    
    #print(np.percentile(box_resized , 2))
    #print(np.percentile(box_resized , 99))
    box_resized = box_resized - np.percentile(box_resized , 2)
    box_resized = (box_resized / np.percentile(box_resized , 99) * 255)
    box_resized = np.clip(box_resized , 0 , 255).astype(np.uint8)
    #print(np.min(box_resized))
    #print(np.max(box_resized))
    
    for gray, blue in zip(box_resized, box_annotated):
        # Ensure frames are 2D arrays
        gray = np.squeeze(gray)
        blue = np.squeeze(blue)

        ## scale intensity on a per-frame basis
        #gray = gray - np.min(gray)
        #gray = (gray / np.max(gray) * 255).astype(np.uint8)
        
        # Create RGB channels
        rgb = np.zeros((gray.shape[0], gray.shape[1], 3), dtype=np.uint8)

        # Grayscale in all channels
        rgb[..., 0] = gray  # Red
        rgb[..., 1] = gray  # Green
        rgb[..., 2] = gray  # Blue

        # Overlay blue mask where blue > 0
        rgb[..., 0][blue.sum(axis = 2) > 0] = blue[...,0][blue.sum(axis = 2) > 0]     # Remove red
        rgb[..., 1][blue.sum(axis = 2) > 0] = blue[...,1][blue.sum(axis = 2) > 0]      # Remove green
        rgb[..., 2][blue.sum(axis = 2) > 0] = blue[...,2][blue.sum(axis = 2) > 0]    # Add blue at annotation locations

        frames_rgb.append(rgb)

    # Save as MP4
    imageio.mimsave(output_path, frames_rgb, fps=fps, codec='libx264', quality=8)
    print(f"Movie saved to {output_path}")



import numpy as np
import pandas as pd

def bleaching_curve(img, traj_df, signal_size=3, background_size=7, post_frames=10):
    """
    Measure signal and background around detected positions in an image stack.

    Parameters:
    - img: np.ndarray of shape (T, Y, X)
    - traj_df: DataFrame with columns ['x', 'y', 'frame']
    - signal_size: width of signal square (pixels)
    - background_size: width of background square (pixels)
    - post_frames: number of frames to continue measuring after trajectory ends

    Returns:
    - traj_df with new columns ['signal', 'background', 'rel_frame']
    """

    half_signal = signal_size // 2
    half_background = background_size // 2
    results = []

    # Drop NaN positions
    traj_valid = traj_df.dropna(subset=['x', 'y', 'frame']).copy()

    # Get last frame and position for post-measurement
    if len(traj_valid) > 0:
        last_x = traj_valid['x'].iloc[-1]
        last_y = traj_valid['y'].iloc[-1]
        last_frame = int(traj_valid['frame'].iloc[-1])

        # Add extra frames at the same position
        post_rows = []
        for f in range(last_frame + 1, min(last_frame + 1 + post_frames, img.shape[0])):
            post_rows.append({
                'x': last_x, 'y': last_y, 'frame': f,
                'ur_trajectory': traj_valid['ur_trajectory'].iloc[-1],
                'condition': traj_valid['condition'].iloc[-1]
            })
        if post_rows:
            traj_valid = pd.concat([traj_valid, pd.DataFrame(post_rows)], ignore_index=True)

    # Measure intensity per frame
    for _, row in traj_valid.iterrows():
        x = int(round(row['x']))
        y = int(round(row['y']))
        f = int(row['frame'])

        if np.isnan(x) or np.isnan(y) or f >= img.shape[0]:
            signal = np.nan
            background = np.nan
        else:
            frame_img = img[f]
            h, w = frame_img.shape

            # Define region bounds safely (clipped to image)
            y1s, y2s = np.clip([y - half_signal, y + half_signal + 1], 0, h)
            x1s, x2s = np.clip([x - half_signal, x + half_signal + 1], 0, w)
            y1b, y2b = np.clip([y - half_background, y + half_background + 1], 0, h)
            x1b, x2b = np.clip([x - half_background, x + half_background + 1], 0, w)

            # Extract patches
            signal_patch = frame_img[y1s:y2s, x1s:x2s]
            background_patch = frame_img[y1b:y2b, x1b:x2b]

            signal = np.mean(signal_patch)
            background = np.median(background_patch)

        result = row.to_dict()
        result['signal'] = signal
        result['background'] = background
        results.append(result)

    df_out = pd.DataFrame(results)

    # Compute rel_frame relative to trajectory end (0 = last frame)
    if not df_out.empty:
        non_nan = df_out[df_out['I0'].notna()]
        last_frame = non_nan['frame'].max()
        df_out['rel_frame'] = df_out['frame'] - last_frame

    return df_out


def normalize_signal_by_baseline(df, pre_window=10):
    """
    For each trajectory, divide 'signal' by the mean signal in the pre_window frames
    before rel_frame = 0 (i.e., rel_frame < 0).
    """
    df = df.copy()

    def normalize(traj):
        # baseline: frames where rel_frame < 0 and >= -pre_window
        mask = (traj['rel_frame'] < 0) & (traj['rel_frame'] >= -pre_window)
        baseline = traj.loc[mask, 'corrected_signal'].mean()
        if pd.isna(baseline) or baseline == 0:
            # avoid division by zero
            traj['signal_norm'] = traj['corrected_signal']
        else:
            traj['signal_norm'] = traj['corrected_signal'] / baseline
        return traj

    df = df.groupby('ur_trajectory', group_keys=False).apply(normalize)
    return df


def aggregate_bleaching_data(rollingMLE , settings , n_particles = 200 , min_obs = 7):

    counts = rollingMLE['ur_trajectory'].value_counts()
    valid_trajectories = counts[counts >= min_obs].index
    rollingMLE_filtered = rollingMLE[rollingMLE['ur_trajectory'].isin(valid_trajectories)].copy()
    bc = pd.DataFrame()
    prev_open_image = ''

    unique_trajs = np.unique(rollingMLE_filtered["ur_trajectory"])
    n = min(n_particles, len(unique_trajs))
    #n = n_particles
    sampled_ids = np.random.choice(unique_trajs, size=n, replace=False)
    sampled_ids = sorted(sampled_ids)
    sampled_trajs = rollingMLE_filtered[rollingMLE_filtered['ur_trajectory'].isin(sampled_ids)]
    sampled_trajs = sampled_trajs.sort_values(by=['ur_trajectory', 'frame'], ignore_index=True)

    # iterate over each sampled trajectory and run your code
    for traj_id in sampled_ids:
        this_traj , boundingBox = isolate_traj(rollingMLE_filtered , traj_id)

        try:
            base_pattern = this_traj['ur_trajectory'].replace(
                r'_traj::\d+', '', regex=True).iloc[0]
            # Try multiple image extensions (nd2 first, then tif/tiff)
            image_file = None
            for ext in ['.nd2', '.tif', '.tiff']:
                matches = glob.glob(os.path.join(
                    settings['io']['data_directory'], base_pattern + ext))
                if matches:
                    image_file = matches[0]
                    break
            if image_file is None:
                continue

            ## only reload the image if neccessary
            if image_file == prev_open_image:
                pass
            else:
                img = pims.open(image_file)
                img_stack = np.array(img)
                prev_open_image = image_file
            bc = pd.concat([bc , bleaching_curve(img_stack , this_traj)])

        except:
           print('Caught an error')
           pass
    return bc


## can handle many conditions, but as currently implemented is only called on one condition at a time
def plot_bleaching_curves(settings , bc ):

    if bc.empty or 'signal' not in bc.columns:
        print('  No bleaching data to plot (empty DataFrame)')
        return

    plot_xlim = settings['plot']['bleach_xlim'] * 1000
    plot_xlim = (settings['saspt']['frame_interval'] * -30 * 1000) , (settings['saspt']['frame_interval'] * 10 * 1000)

    output_dirname = os.path.join(settings['io']['plot_directory'] , 'bleaching_curves')
    os.makedirs(output_dirname, exist_ok=True)

    bc['corrected_signal'] = bc['signal'] - bc['background']
    # Apply to bc
    bc_norm = normalize_signal_by_baseline(bc, pre_window=10)
    bc_norm['signal_smooth'] = bc_norm.groupby('ur_trajectory')['signal_norm'] \
                                      .transform(lambda x: x.rolling(3, min_periods=1).mean())
    median_fluorescence = (
        bc_norm
        .groupby(['condition', 'rel_frame'], as_index=False)
        .agg(median_signal=('signal_norm', 'median'))
    )

    for cond in bc_norm['condition'].unique():
        output_filename = os.path.join(output_dirname, f"{cond}.pdf")
        print(f"Writing to {output_filename}")

        # Filter data for this condition
        df_cond = bc_norm[bc_norm['condition'] == cond]
        df_median = median_fluorescence[median_fluorescence['condition'] == cond]
        
        df_cond['rel_frame'] = df_cond['rel_frame'] * settings['saspt']['frame_interval'] * 1000
        df_median['rel_frame'] = df_median['rel_frame'] * settings['saspt']['frame_interval'] * 1000

        # Create figure
        plt.figure(figsize=(8, 6))

        # Plot each trajectory separately
        for traj_id, traj_data in df_cond.groupby('ur_trajectory'):
            plt.plot(
                traj_data['rel_frame'],
                traj_data['signal_smooth'] ,
                color = 'darkgray' , 
                alpha = 0.25
            )
        plt.plot(
            df_median['rel_frame'],
            df_median['median_signal'] ,
            color = 'black'
        )   

        # Axis limits
        plt.xlim(plot_xlim[0], plot_xlim[1])
        plt.ylim(-0.5, 3.5)

        # Labels and title
        plt.xlabel("Time rel bleaching (ms)")
        plt.ylabel("Signal")
        plt.title(f"Condition: {cond}")

        # Save figure
        plt.tight_layout()
        plt.savefig(output_filename, dpi=300)
        plt.close()




## take a rollingMLE object and display movies of the tracked particles
def generate_movies(rollingMLE , settings , min_obs = 8 , overlay = True):
    counts = rollingMLE['ur_trajectory'].value_counts()
    valid_trajectories = counts[counts >= min_obs].index

    rollingMLE_filtered = rollingMLE[rollingMLE['ur_trajectory'].isin(valid_trajectories)].copy()
    condition = rollingMLE_filtered['ur_trajectory'].str.replace(r'_\d+_traj::\d+', '', regex=True).unique()[0]
    
    unique_trajs = np.unique(rollingMLE_filtered["ur_trajectory"])
    n = min(100, len(unique_trajs))

    sampled_trajs = np.random.choice(unique_trajs, size=n, replace=False)

    output_dirname = os.path.join(settings['io']['movie_directory'] , condition)
    os.makedirs(output_dirname, exist_ok=True)
    prev_open_image = ''

    # iterate over each sampled trajectory and run your code
    for traj_id in sampled_trajs:
        this_traj , boundingBox = isolate_traj(rollingMLE_filtered , traj_id)
        output_filename = f"{this_traj['ur_trajectory'][0]}.mp4"
        output_filename = os.path.join(output_dirname , output_filename)

        #print(f"Writing to     {output_filename}")

        try:
            base_pattern = this_traj['ur_trajectory'].replace(
                r'_traj::\d+', '', regex=True).iloc[0]
            # Try multiple image extensions (nd2 first, then tif/tiff)
            image_file = None
            for ext in ['.nd2', '.tif', '.tiff']:
                matches = glob.glob(os.path.join(
                    settings['io']['data_directory'], base_pattern + ext))
                if matches:
                    image_file = matches[0]
                    break
            if image_file is None:
                continue

             ## only reload the image if neccessary
            if image_file == prev_open_image:
                pass
                #print('skipping open new image')
            else:
                #print(f'Old image:: {prev_open_image} ... new image:: {image_file}')
                img = pims.open(image_file)
                #half_size = 15  # since 15x15 has half-width 7
                img_stack = np.array(img)
                prev_open_image = image_file


            ## isolte bounding box around the particle
            box = img_stack[boundingBox['timeBox'][0] : boundingBox['timeBox'][1], boundingBox['spaceBox']['ymin'] : boundingBox['spaceBox']['ymax'] ,
                                 boundingBox['spaceBox']['xmin'] : boundingBox['spaceBox']['xmax']]

            ## generate a render-able image (resized as needed) and a printable annotation layer matching the resolution of the image
            box_resized , box_resized_anotations = annotate_trajectory_building_tail(box , this_traj , boundingBox['spaceBox'] , upscale=8, spot_radius = 12 , lwd = 3)

            if overlay == True:
                save_overlay_movie(box_resized , box_resized_anotations , output_filename)
            else:
                box_resized_annotations = np.zeros(box_resized_annotations)
                save_overlay_movie(box_resized , box_resized_anotations , output_filename)

                #output_movie = os.path.join(settings['io']['movie_directory'] , 'movie.tif')
                #save_overlay_tiff(box_resized , box_resized_anotations , output_movie) ## this is 24mb... save these sparingly
        except Exception as e:
            print(f'error generating movie for {traj_id}: {e}')