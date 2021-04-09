%%
%% ======================================================================
%%
%% Sample configurations file for WPE
%%
%% Copyright (c) 2015 Nippon Telegraph and Telephone corporation (NTT).
%% All rights reserved.
%% By Takuya Yoshioka, Marc Delcroix 24-06-2015.
%% ======================================================================


% Basic parameters
%----------------------------------------------------------------------

fs    = 16000;     %% Sampling frequency

num_mic = 2;       %% Number of channels

num_out = num_mic; %% Number of outputs (should be <= microphones).
		   %% set to 1 if output a single channel 

blk_len = 100;      %% Block length (in sec). Set to a large value for
                   %% utterance-batch processing
		   
opt_blk_sz = 1;    %% Optimize the block size for block batch processing
                   %%  to speedup computations

% Signal analysis parameters
%----------------------------------------------------------------------

analy_param = struct('win_size'  , 1024, ...
		     'shift_size', 512, ...
		     'win'       , hanning(1024));


% Dereverberation parameters 
%----------------------------------------------------------------------

%% Parameters of prediction filter
%% [number of filter coefficient; prediction delay; upper frequency]
%% It is possible to set different filter parameters for different
%% frequency bands e.g. filter length of 10 up to 500 Hz and 6 for the rest 
%% channel_setup = [10, 6;
%%                 3, 3;
%%                 500, inf]
channel_setup = [5; ...
		 2; ...
		 inf];

%% Dereverberation filter configuration
%% 'channel_setup' consists of the prediction filter settings set above
%% 'p_channel'     sets the index of the target channel for prediction
%% 'speech_order'  sets the speech lpc order
%%                 "speech_order = inf" disables LPC speech modeling
ssd_param = struct('channel_setup', channel_setup, ...
		   'p_channel'    , [1 : num_out], ...
		   'speech_order' , 20);

%% Optimization parameters
%% 'max_iter'  number of iterations
%% 'spcorr'    structure of the speech covariance matrix ('scaleye'
%%             corresponds to diagonal)
%% 'scaling'   gain between input and output
%% 'forget'    forgetting factor for the correlation matrix (the larger
%%             the more the past observations are remembered)
ssd_conf = struct('max_iter', 1, ...
		  'spcorr'  , 'scaleye', ...
		  'scaling' , 1, ...
		  'forget'  , 0.7);

%% Enhancement configuration
%% 'method'   sets the method used to suppress late reverberation (either
%%            linear filtering 'lti' or spectral subtraction 'ss'
%% 'osub'     oversubtraction (for 'ss' only)
%% 'scal'     scaling (for 'ss' only)
%% 'floor'    flooring (for 'ss' only)
enh_conf = struct('method', 'lti', ...
		  'osub'  , 1.0, ...
		  'scal'  , 2, ...
		  'floor' , -80);

