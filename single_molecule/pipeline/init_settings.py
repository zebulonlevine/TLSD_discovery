import numpy as np
from saspt import RBME
import os
import pims
import glob


def get_default_settings():
    settings = {}
    settings['io'] = {}
    settings['quot'] = {}
    settings['saspt'] = {}
    settings['plot'] = {}
    
    settings['io']['home_directory'] = ""
    settings['io']['data_filepath'] = "raw_data"
    settings['io']['split_dimension'] = 0 ## pixel dimension to split image into -- this is useful for dense / large images, i.e. for slow SPT
    settings['io']['forceQUOT'] = True 
    settings['io']['plotOnly'] = False 
    # settings['io']['traj_filestring_suffix'] = '_splitsize10'  ## optional suffix appended to tracking_output_q=X dir name; lets parallel runs (e.g., different splitsize) coexist without clobbering

    settings['quot']['filter'] = {}
    settings['quot']['filter']['start'] = 0 
    settings['quot']['filter']['method'] = 'identity'
    settings['quot']['filter']['chunk_size'] = 100    
    
    settings['quot']['detect'] = {}
    settings['quot']['detect']['method'] = 'llr'
    settings['quot']['detect']['k'] = 2.0 ## this is related to expected spot size i think
    settings['quot']['detect']['w'] = 15 ## window dimension over which to look for spots
    settings['quot']['detect']['t'] = 15.00 ## likelihood ratio testing the model that the window contains a spot vs doesn't
    
    '''
    ## try LoG filtereing
    settings['quot']['detect'] = {}
    settings['quot']['detect']['method'] = 'log'
    settings['quot']['detect']['k'] = 1 ## this is related to expected spot size i think
    settings['quot']['detect']['w'] = 11 ## window dimension over which to look for spots
    settings['quot']['detect']['t'] = 20 ## likelihood ratio testing the model that the window contains a spot vs doesn't
    
    settings['quot']['detect'] = {}
    settings['quot']['detect']['method'] = 'hess_det_var'
    settings['quot']['detect']['k'] = 1.5 ## this is related to expected spot size i think
    settings['quot']['detect']['w0'] = 15 ## window dimension over which to look for spots
    settings['quot']['detect']['w1'] = 9 ## window dimension over which to look for spots
    settings['quot']['detect']['t'] = 20.0 ## likelihood ratio testing the model that the window contains a spot vs doesn't
    '''
    
    ### note that i often get error flags for my spot localization, which means the spot intensity is below the background.
    ### that happens when spots have locally high background, or when the background varies a lot between files.
    ### i dont think it matters... the "error_flag=1" spots still seem real. The actual spot calling is all based on shape, not intensity
    settings['quot']['localize'] = {}
    settings['quot']['localize']['method'] = 'ls_int_gaussian'
    settings['quot']['localize']['window_size'] = 15 ## window size of 7 was bad. 11 is decent. 15 was best by a hair.
    settings['quot']['localize']['sigma'] = 2.5 
    settings['quot']['localize']['ridge'] = 0.00001
    settings['quot']['localize']['max_iter'] = 50
    settings['quot']['localize']['damp'] = 1 # higher damp means faster convergence
    #settings['quot']['localize']['camera_gain'] = 100
    settings['quot']['localize']['camera_bg'] = 100
    
    settings['quot']['track'] = {}
    settings['quot']['track']['method'] = 'euclidean' ## conservative is too conservative. euclidean should match trackmate
    settings['quot']['track']['pixel_size_um'] = 0.11
    settings['quot']['track']['frame_interval'] = 0.01 ## frame interval in seconds
    settings['quot']['track']['search_radius'] = 1.2 ## linking distance in microns
    settings['quot']['track']['max_blinks'] = 0
    settings['quot']['track']['min_I0'] = 100.0 ## "spot intensity above background." might be a hard parameter to tune
    settings['quot']['track']['scale'] = 1.0
    
    '''
    settings['trackmate'] = {}
    settings['trackmate']['radius'] = 0.2
    settings['trackmate']['quality'] = 14.0
    settings['trackmate']['linking_distance'] = 1.6
    settings['trackmate']['channel'] = 1
    '''

    settings['saspt']['likelihood_type'] = RBME
    settings['saspt']['pixel_size_um'] = 0.11
    settings['saspt']['frame_interval'] = 0.01  ## frame interval in seconds — must match quot/track/frame_interval default
    settings['saspt']['focal_depth'] = 0.1
    settings['saspt']['progress_bar'] = 'False'
    settings['saspt']['splitsize'] = 10
    settings['saspt']['sample_size'] = 1000000
    settings['saspt']['num_workers'] = 4
    settings['saspt']['diff_coefs'] = np.power(10,(np.linspace(-3,2,150)))
    settings['saspt']['loc_errors'] = np.linspace(0.025 , 0.032 , 5)
    settings['saspt']['start_frame'] = 0
    settings['saspt']['save_subtrack_posteriors'] = False  # opt-in: save per-subtrack posterior NPZ files

    settings['plot']['xlim'] = (-3,2)
    settings['plot']['bleach_xlim'] = (-300,50)
    settings['plot']['ylim'] = [ 0 , 2 ]
    settings['plot']['nbins'] = 50
    settings['plot']['tick_positions'] = (-3, -2, -1, 0, 1, 2)
    settings['plot']['mult_on_sd'] = 5.0
    settings['plot']['barGraphBreaks'] = [0 , 0.05, 0.15, 0.5, 2, 20]
    settings['plot']['barGraphLabels'] = ['immobile', 'membrane', 'polysome', 'free', 'non-ribosome']

    settings['extrack'] = {}
    settings['extrack']['enabled'] = True               # run ExTrack state analysis in pipeline
    settings['extrack']['nb_states'] = None              # None = len(barGraphLabels); 'auto' = BIC scan; int = fixed
    settings['extrack']['min_states'] = 2               # BIC scan minimum states
    settings['extrack']['max_states'] = 5               # BIC scan maximum states
    settings['extrack']['n_iter'] = 3                   # fitting iterations
    settings['extrack']['frame_len'] = 4                # exact probability window (nb_states^frame_len sequences; 4 is safe for 6 states)
    settings['extrack']['workers'] = 4                  # parallel workers (match SLURM cpus-per-task)
    settings['extrack']['max_length'] = 25              # truncate tracks (25 frames ≈ 0.15-0.25s; sufficient for D + transitions)
    settings['extrack']['max_tracks'] = 5000            # subsample to this many tracks per condition (None = all)
    settings['extrack']['fit_max_tracks'] = 1000        # tracks for parameter fitting (subsample); None = use max_tracks
    settings['extrack']['predict_max_length'] = 100     # max trajectory length for state prediction (longer than fit)
    settings['extrack']['fix_loc_error'] = None         # fix LocErr to this value (um), or None
    settings['extrack']['cell_dims'] = None             # FOV dimensions in um, or None
    settings['extrack']['loc_error_um'] = 0.030         # initial localization error estimate
    settings['extrack']['min_frames'] = 6               # minimum trajectory length for ExTrack
    settings['extrack']['per_category_auto_states'] = True  # BIC per category sub-analysis

    return(settings)


def update_default_settings_for_analysis(settings , home_directory):
    
    settings['io']['home_directory'] = home_directory
    settings['io']['data_directory'] = os.path.join(settings['io']['home_directory'] , settings['io']['data_filepath'])

    traj_filestring = f"tracking_output_q={str(float(settings['quot']['detect']['t'])).replace('.','p')}"
    # Optional suffix to coexist with sibling runs (different splitsize, etc.)
    suffix = settings['io'].get('traj_filestring_suffix', '')
    if suffix:
        traj_filestring = f"{traj_filestring}{suffix}"
    settings['io']['analysis_directory'] = os.path.join(settings['io']['home_directory'] , traj_filestring)
    settings['io']['traj_directory'] = os.path.join(settings['io']['analysis_directory'] , 'trajectories')
    settings['io']['post_directory'] = os.path.join(settings['io']['analysis_directory'] , 'posterior')
    settings['io']['MLE_directory'] = os.path.join(settings['io']['analysis_directory'] , 'MLE')
    settings['io']['plot_directory'] = os.path.join(settings['io']['analysis_directory'] , 'plots')
    #settings['io']['split_directory'] = os.path.join(settings['io']['analysis_directory'] , 'split_directory')
    #settings['io']['split_traj_directory'] = os.path.join(settings['io']['analysis_directory'] , 'split_trajectories')
    settings['io']['rolling_window_directory'] = os.path.join(settings['io']['analysis_directory'] , 'rolling_windows')
    settings['io']['movie_directory'] = os.path.join(settings['io']['analysis_directory'] , 'movies')
    settings['io']['subtrack_posterior_directory'] = os.path.join(settings['io']['analysis_directory'] , 'subtrack_posteriors')
    # Combination-stage I/O. Used by combination_functions.redirect_settings_for_combo
    # to redirect saSPT outputs without per-call argument plumbing.
    settings['io']['traj_combinations_directory'] = os.path.join(settings['io']['analysis_directory'] , 'traj_combinations')
    settings['io']['post_combinations_directory'] = os.path.join(settings['io']['analysis_directory'] , 'posterior_combinations')
    settings['io']['MLE_combinations_directory'] = os.path.join(settings['io']['analysis_directory'] , 'MLE_combinations')
    settings['io']['rolling_window_combinations_directory'] = os.path.join(settings['io']['analysis_directory'] , 'rolling_windows_combinations')
    settings['io']['plot_combinations_directory'] = os.path.join(settings['io']['analysis_directory'] , 'plots_combinations')

    # Megapool-stage I/O. Used by megapool_functions.redirect_settings_for_megapool
    # to redirect saSPT outputs without per-call argument plumbing.
    settings['io']['traj_megapools_directory'] = os.path.join(settings['io']['analysis_directory'] , 'traj_megapools')
    settings['io']['post_megapools_directory'] = os.path.join(settings['io']['analysis_directory'] , 'posterior_megapools')
    settings['io']['MLE_megapools_directory'] = os.path.join(settings['io']['analysis_directory'] , 'MLE_megapools')


    settings['quot']['filter']['start'] = 0

    os.makedirs(settings['io']['analysis_directory'], exist_ok=True)
    os.makedirs(settings['io']['traj_directory'], exist_ok=True)
    os.makedirs(settings['io']['post_directory'], exist_ok=True)
    os.makedirs(settings['io']['MLE_directory'], exist_ok=True)
    os.makedirs(settings['io']['plot_directory'], exist_ok=True)
    #os.makedirs(settings['io']['split_directory'] , exist_ok=True)
    #os.makedirs(settings['io']['split_traj_directory'] , exist_ok=True)
    os.makedirs(settings['io']['rolling_window_directory'] , exist_ok=True)
    os.makedirs(settings['io']['movie_directory'] , exist_ok=True)
    os.makedirs(settings['io']['subtrack_posterior_directory'] , exist_ok=True)
    os.makedirs(settings['io']['traj_combinations_directory'] , exist_ok=True)
    os.makedirs(settings['io']['post_combinations_directory'] , exist_ok=True)
    os.makedirs(settings['io']['MLE_combinations_directory'] , exist_ok=True)
    os.makedirs(settings['io']['rolling_window_combinations_directory'] , exist_ok=True)
    os.makedirs(settings['io']['plot_combinations_directory'] , exist_ok=True)
    os.makedirs(settings['io']['traj_megapools_directory'] , exist_ok=True)
    os.makedirs(settings['io']['post_megapools_directory'] , exist_ok=True)
    os.makedirs(settings['io']['MLE_megapools_directory'] , exist_ok=True)

    return(settings)


def update_settings_with_image_metadata(settings , cond = None):
    
    
    #arbitrary_image = os.listdir(settings['io']['data_directory'])[0]
    
    #    if cond == None:
    #        arbitrary_image = glob.glob(f"{settings['io']['data_directory']}/*.nd2")[0]
    #    else:
    #        arbitrary_image = glob.glob(f"{settings['io']['data_directory']}/*{cond}*.nd2")[0]
    
    # Skip metadata extraction for TIFF files (already have frame_interval in settings)
    data_dir = settings['io']['data_directory']
    tiff_files = glob.glob(f"{data_dir}/*{cond}*.tif")
    if tiff_files:
        # Sync saspt frame_interval from quot value (no ND2 metadata available)
        settings['saspt']['frame_interval'] = settings['quot']['track']['frame_interval']
        print(f"TIFF input detected - using frame_interval from settings: {settings['saspt']['frame_interval']}")
        return settings
    
    
    if cond == None:
        arbitrary_image = glob.glob(f"{settings['io']['data_directory']}/*.nd2")
        if not arbitrary_image:
                arbitrary_image = glob.glob(f"{settings['io']['data_directory']}/*.tif")
        arbitrary_image = arbitrary_image[0]
    else:
        arbitrary_image = glob.glob(f"{settings['io']['data_directory']}/*{cond}*.nd2")
        if not arbitrary_image:
                arbitrary_image = glob.glob(f"{settings['io']['data_directory']}/*{cond}*.tif")
        arbitrary_image = arbitrary_image[0]
    
    #img = pims.open(os.path.join(settings['io']['data_directory'] ,arbitrary_image))
    img = pims.open(arbitrary_image)
    frame = img[0]
    bg_level = np.median(frame)

    dt = 1 / img.frame_rate
    
    settings['quot']['localize']['camera_bg'] = bg_level
    settings['quot']['track']['frame_interval'] = dt
    settings['saspt']['frame_interval'] = dt
    
    return(settings)
