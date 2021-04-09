%
% The main script.
%
clear all;
addpath('bss_eval');

% sample rate
fs=16000;
% reverberation time
RT60_SET=[0.6, 0.3, 0.8];
% SIR (dB) of the test data
SIR_SET=[0];
% SER (dB) of the test data
SER_SET=[0, -10];
% no. of experiments
numrpts=20;
% for algorithm convergence
UTT_REPEAT=2;
% seed
rngseed=0;

for i=1:length(RT60_SET)
for j=1:length(SIR_SET)
for k=1:length(SER_SET)
    rt60=RT60_SET(i);
    sirdb=SIR_SET(j);
    serdb=SER_SET(k);
    savedir=['output/rt' num2str(rt60) '_sir' num2str(sirdb) '_ser' num2str(serdb)];
    if ~exist(savedir)
        mkdir(savedir);
    end
    
for rptcount=1:numrpts
    fprintf(1, 'rptcount = %d ...\n', rptcount);
    
    %% generate test data
    fprintf(1, 'generate test data\n');
    
    rngseed=(i-1)*1e4+(j-1)*1e3+(k-1)*1e2+rptcount;
    [testdata, target_interf_echo] = simulateDataImageRIR(sirdb, serdb, rt60, UTT_REPEAT, rngseed);
    siglen=size(testdata, 1)/UTT_REPEAT;
    audiowrite([savedir '/input_data_', num2str(rptcount), '.wav'], testdata(end-siglen+1:end, :), fs);
    audiowrite([savedir '/target_interf_echo', num2str(rptcount), '.wav'], target_interf_echo(end-siglen+1:end, :), fs);
    
    %% perform separation
    nummics=2;
    numrefs=1;

    fprintf(1, '1. WPE-NLMS-BSS\n');
    out=modeWPE_NLMS_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite([savedir '/01_out_wpe_nlms_bss_', num2str(rptcount), '.wav'], out(end-siglen+1:end, :), fs);
    
    fprintf(1, '2. NLMS-WPE-BSS\n');
    out=modeNLMS_WPE_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite([savedir '/02_out_nlms_wpe_bss_', num2str(rptcount), '.wav'], out(end-siglen+1:end, :), fs);
    
    fprintf(1, '3. DR-AEC-BSS\n');
    out=modeDR_AEC_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite([savedir '/03_out_dr_aec_bss_', num2str(rptcount), '.wav'], out(end-siglen+1:end, :), fs);
    
    fprintf(1, '4. AEC-DR-BSS\n');
    out=modeAEC_DR_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite([savedir '/04_out_aec_dr_bss_', num2str(rptcount), '.wav'], out(end-siglen+1:end, :), fs);
    
    fprintf(1, '5. DRAEC-BSS\n');
    out=modeDRAEC_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite([savedir '/05_out_draec_bss_', num2str(rptcount), '.wav'], out(end-siglen+1:end, :), fs);
    
    fprintf(1, '6. Joint-SS\n');
    out=modeBSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite([savedir '/06_out_joint_ss_', num2str(rptcount), '.wav'], out(end-siglen+1:end, :), fs);
end

Main_score;

end
end
end
