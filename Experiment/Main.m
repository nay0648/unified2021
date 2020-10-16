%
% The main script.
%
clear;
addpath('bss_eval');

% sample rate
fs=16000;
% SNR (dB) of the test data
snrdb=0;
% SDR (dB) of the test data
serdb=0;
% no. of experiments
numrpts=100;

sdri_wpe_nlms_bss=zeros(numrpts, 1);
sdri_nlms_wpe_bss=zeros(numrpts, 1);
sdri_dr_aec_bss=zeros(numrpts, 1);
sdri_aec_dr_bss=zeros(numrpts, 1);
sdri_draec_bss=zeros(numrpts, 1);
sdri_bss=zeros(numrpts, 1);

siri_wpe_nlms_bss=zeros(numrpts, 1);
siri_nlms_wpe_bss=zeros(numrpts, 1);
siri_dr_aec_bss=zeros(numrpts, 1);
siri_aec_dr_bss=zeros(numrpts, 1);
siri_draec_bss=zeros(numrpts, 1);
siri_bss=zeros(numrpts, 1);

sisdri_wpe_nlms_bss=zeros(numrpts, 1);
sisdri_nlms_wpe_bss=zeros(numrpts, 1);
sisdri_dr_aec_bss=zeros(numrpts, 1);
sisdri_aec_dr_bss=zeros(numrpts, 1);
sisdri_draec_bss=zeros(numrpts, 1);
sisdri_bss=zeros(numrpts, 1);

for rptcount=1:numrpts
    fprintf(1, 'rptcount = %d ...\n', rptcount);
    
    %% generate test data
    fprintf(1, 'generate test data\n');
    
    rt60=0.4;
    [testdata, target_interf_echo] = simulateDataImageRIR(snrdb, serdb, rt60);
    audiowrite(['output/input_data_', num2str(rptcount), '.wav'], testdata, fs);

    %% performance indices for input signals
    [sdr, sir, sar, perm]=bss_eval_sources([testdata(:, 1:2), target_interf_echo(:, 3)]', target_interf_echo');
    sdr_input=sdr(1);
    sir_input=sir(1);
    sisdr_input=sisdr(target_interf_echo(:, 1), testdata(:, 1:2));
    
    %% perform separation
    nummics=2;
    numrefs=1;

    fprintf(1, '1. WPE-NLMS-BSS\n');
    out=modeWPE_NLMS_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite(['output/01_out_wpe_nlms_bss_', num2str(rptcount), '.wav'], out, fs);
    [sdr, sir, sar, perm]=bss_eval_sources([out, target_interf_echo(:, 3)]', target_interf_echo');
    sdri_wpe_nlms_bss(rptcount)=sdr(1)-sdr_input;
    siri_wpe_nlms_bss(rptcount)=sir(1)-sir_input;
    sisdri_wpe_nlms_bss(rptcount)=sisdr(target_interf_echo(:, 1), out)-sisdr_input;
    
    fprintf(1, '2. NLMS-WPE-BSS\n');
    out=modeNLMS_WPE_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite(['output/02_out_nlms_wpe_bss_', num2str(rptcount), '.wav'], out, fs);
    [sdr, sir, sar, perm]=bss_eval_sources([out, target_interf_echo(:, 3)]', target_interf_echo');
    sdri_nlms_wpe_bss(rptcount)=sdr(1)-sdr_input;
    siri_nlms_wpe_bss(rptcount)=sir(1)-sir_input;
    sisdri_nlms_wpe_bss(rptcount)=sisdr(target_interf_echo(:, 1), out)-sisdr_input;
    
    fprintf(1, '3. DR-AEC-BSS\n');
    out=modeDR_AEC_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite(['output/03_out_dr_aec_bss_', num2str(rptcount), '.wav'], out, fs);
    [sdr, sir, sar, perm]=bss_eval_sources([out, target_interf_echo(:, 3)]', target_interf_echo');
    sdri_dr_aec_bss(rptcount)=sdr(1)-sdr_input;
    siri_dr_aec_bss(rptcount)=sir(1)-sir_input;
    sisdri_dr_aec_bss(rptcount)=sisdr(target_interf_echo(:, 1), out)-sisdr_input;
    
    fprintf(1, '4. AEC-DR-BSS\n');
    out=modeAEC_DR_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite(['output/04_out_aec_dr_bss_', num2str(rptcount), '.wav'], out, fs);
    [sdr, sir, sar, perm]=bss_eval_sources([out, target_interf_echo(:, 3)]', target_interf_echo');
    sdri_aec_dr_bss(rptcount)=sdr(1)-sdr_input;
    siri_aec_dr_bss(rptcount)=sir(1)-sir_input;
    sisdri_aec_dr_bss(rptcount)=sisdr(target_interf_echo(:, 1), out)-sisdr_input;
    
    fprintf(1, '5. DRAEC-BSS\n');
    out=modeDRAEC_BSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite(['output/05_out_draec_bss_', num2str(rptcount), '.wav'], out, fs);
    [sdr, sir, sar, perm]=bss_eval_sources([out, target_interf_echo(:, 3)]', target_interf_echo');
    sdri_draec_bss(rptcount)=sdr(1)-sdr_input;
    siri_draec_bss(rptcount)=sir(1)-sir_input;
    sisdri_draec_bss(rptcount)=sisdr(target_interf_echo(:, 1), out)-sisdr_input;
    
    fprintf(1, '6. BSS\n');
    out=modeBSS(nummics, numrefs, testdata);
    out(size(testdata, 1)+1:end, :)=[];
    audiowrite(['output/06_out_bss_', num2str(rptcount), '.wav'], out, fs);
    [sdr, sir, sar, perm]=bss_eval_sources([out, target_interf_echo(:, 3)]', target_interf_echo');
    sdri_bss(rptcount)=sdr(1)-sdr_input;
    siri_bss(rptcount)=sir(1)-sir_input;
    sisdri_bss(rptcount)=sisdr(target_interf_echo(:, 1), out)-sisdr_input;
end

figure(1);
boxplot([sdri_wpe_nlms_bss, sdri_nlms_wpe_bss, sdri_dr_aec_bss, sdri_aec_dr_bss, sdri_draec_bss, sdri_bss],...
    'labels', {'WPE-NLMS-BSS', 'NLMS-WPE-BSS', 'DR-AEC-BSS', 'AEC-DR-BSS', 'DRAEC-BSS', 'BSS'});
ylabel('SDR Improvement');

figure(2);
boxplot([siri_wpe_nlms_bss, siri_nlms_wpe_bss, siri_dr_aec_bss, siri_aec_dr_bss, siri_draec_bss, siri_bss],...
    'labels', {'WPE-NLMS-BSS', 'NLMS-WPE-BSS', 'DR-AEC-BSS', 'AEC-DR-BSS', 'DRAEC-BSS', 'BSS'});
ylabel('SIR Improvement');

figure(3);
boxplot([sisdri_wpe_nlms_bss, sisdri_nlms_wpe_bss, sisdri_dr_aec_bss, sisdri_aec_dr_bss, sisdri_draec_bss, sisdri_bss],...
    'labels', {'WPE-NLMS-BSS', 'NLMS-WPE-BSS', 'DR-AEC-BSS', 'AEC-DR-BSS', 'DRAEC-BSS', 'BSS'});
ylabel('SI-SDR Improvement');
