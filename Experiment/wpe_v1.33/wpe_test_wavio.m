%%
%% ======================================================================
%% wpe.m
%%
%% Perform dereverberation using WPE.
%%
%% Copyright (c) 2015  Nippon Telegraph and Telephone corporation (NTT).
%% All rights reserved
%% By Takuya Yoshioka, Marc Delcroix 24-06-2015.
%% ======================================================================
%%



cfgs      = 'settings/local.m'
scp       = 'settings/sample.scp'
arrayname = 'settings/arrayname.lst';


% Check the args
%----------------------------------------------------------------------

if ~exist('cfgs', 'var')
  disp('Matlab config file must be given');
  return;
end
if ~exist('scp', 'var')
  disp('SCP file must be given');
  return;
end
if exist('arrayname', 'var')
  use_multi = 1;
else
  use_multi = 0;
end

if use_multi
  disp('WPE with multi-channel mode\n');
  fprintf('ARRAYNAME = %s\n', arrayname);
  fprintf('CFGS      = %s\n', cfgs);
  fprintf('SCP       = %s\n', scp);
else
  disp('WPE with single-channel mode\n');
  fprintf('CFGS = %s\n', cfgs);
  fprintf('SCP  = %s\n', scp);
end

% Create the list of the target files.
%----------------------------------------------------------------------

fid = fopen(scp, 'r');
if fid == -1
  fprintf('Open failed: %s\n', scp);
  return;
end

num_file = 0;
ilist    = cell(10000, 1);
olist    = cell(10000, 1);

while ~feof(fid)
  num_file = num_file + 1;

  l      = strtrim(fgetl(fid));
  [l, r] = strtok(l);
  
  ilist{num_file} = strtrim(l);
  olist{num_file} = strtrim(r(2 : end));
end

fclose(fid);

ilist = ilist(1 : num_file);
olist = olist(1 : num_file);


% Multi-channel file name expansion
%----------------------------------------------------------------------

if use_multi
  num_line = 0;
  fid = fopen(arrayname, 'r');
  
  while ~feof(fid)
    num_line = num_line + 1;
    
    [l, r] = strtok(strtrim(fgetl(fid)));
    renames{num_line, 1} = strtrim(l);
    renames{num_line, 2} = strtrim(r);
  end
  fclose(fid);
end


% Process each utterance.
%----------------------------------------------------------------------

for fidx = 1 : num_file  
  ifname{1} = ilist{fidx};
  ofname{1} = olist{fidx};
  
  if use_multi
    for ren = 1 : size(renames, 1)
      ifname{1 + ren} = strrep(ifname{1}, renames{ren, 1}, renames{ren, 2});
      ofname{1 + ren} = strrep(ofname{1}, renames{ren, 1}, renames{ren, 2});
    end
  end
  
  odir = fileparts(ofname{1});
  cmd = ['mkdir -p -v ', odir];
  system(cmd);
  
  %% Show progress.
  fprintf('Processing %d-th file (%d files to process)\n', fidx, num_file);
  for ch = 1 : length(ifname)
    fprintf('  %s -> %s \n', ifname{ch}, ofname{ch});
  end
  
  
  
  % Perform dereverberation with WPE
  %----------------------------------------------------------------------
  wpe_wavio(ifname, ofname, cfgs);  
  
end

