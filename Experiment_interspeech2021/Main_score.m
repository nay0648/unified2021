% clear all
% numrpts=30
% savedir=

sdri_wpe_nlms_bss=zeros(numrpts, 1);
sdri_nlms_wpe_bss=zeros(numrpts, 1);
sdri_dr_aec_bss=zeros(numrpts, 1);
sdri_aec_dr_bss=zeros(numrpts, 1);
sdri_draec_bss=zeros(numrpts, 1);
sdri_joint_ss=zeros(numrpts, 1);

siri_wpe_nlms_bss=zeros(numrpts, 1);
siri_nlms_wpe_bss=zeros(numrpts, 1);
siri_dr_aec_bss=zeros(numrpts, 1);
siri_aec_dr_bss=zeros(numrpts, 1);
siri_draec_bss=zeros(numrpts, 1);
siri_joint_ss=zeros(numrpts, 1);

sisdri_wpe_nlms_bss=zeros(numrpts, 1);
sisdri_nlms_wpe_bss=zeros(numrpts, 1);
sisdri_dr_aec_bss=zeros(numrpts, 1);
sisdri_aec_dr_bss=zeros(numrpts, 1);
sisdri_draec_bss=zeros(numrpts, 1);
sisdri_joint_ss=zeros(numrpts, 1);

% signal+interference to interference ratio
siiri_wpe_nlms_bss=zeros(numrpts, 1);
siiri_nlms_wpe_bss=zeros(numrpts, 1);
siiri_dr_aec_bss=zeros(numrpts, 1);
siiri_aec_dr_bss=zeros(numrpts, 1);
siiri_draec_bss=zeros(numrpts, 1);
siiri_joint_ss=zeros(numrpts, 1);

% signal+interference+echo to interference+echo ratio 
sieri_wpe_nlms_bss=zeros(numrpts, 1);
sieri_nlms_wpe_bss=zeros(numrpts, 1);
sieri_dr_aec_bss=zeros(numrpts, 1);
sieri_aec_dr_bss=zeros(numrpts, 1);
sieri_draec_bss=zeros(numrpts, 1);
sieri_joint_ss=zeros(numrpts, 1);

% signal+echo to echo ratio
seeri_wpe_nlms_bss=zeros(numrpts, 1);
seeri_nlms_wpe_bss=zeros(numrpts, 1);
seeri_dr_aec_bss=zeros(numrpts, 1);
seeri_aec_dr_bss=zeros(numrpts, 1);
seeri_draec_bss=zeros(numrpts, 1);
seeri_joint_ss=zeros(numrpts, 1);

for rptcount=1:numrpts
    testdata=audioread([savedir '/input_data_', num2str(rptcount), '.wav']);
    target_interf_echo=audioread([savedir '/target_interf_echo', num2str(rptcount), '.wav']);
    tblock_sdr=size(testdata, 1)/4*2+1:size(testdata, 1)/4*3;
    tblock_siir_x=1:size(testdata, 1)/4;
    tblock_siir_y=size(testdata, 1)/4+1:size(testdata, 1)/4*2;
    tblock_sier_x=size(testdata, 1)/4*3+1:size(testdata, 1)/4*4;
    tblock_sier_y=tblock_sdr;

    % performance indices for input signals
    [sdr, sir, sar, perm]=bss_eval_sources(testdata(tblock_sdr, 1:2)', target_interf_echo(tblock_sdr, 1:2)');
    sdr_input=sdr(1);
    sir_input=sir(1);
    sisdr_input=sisdr(target_interf_echo(tblock_sdr, 1), testdata(tblock_sdr, 1:2));
    siir_input=10 * log10 (sum(testdata(tblock_siir_y, 1).^2) / sum(testdata(tblock_siir_x, 1).^2));
    sier_input=10 * log10 (sum(testdata(tblock_sier_y, 1).^2) / sum(testdata(tblock_sier_x, 1).^2));
    
    out=audioread([savedir '/01_out_wpe_nlms_bss_', num2str(rptcount), '.wav']);
    [sdr, sir, sar, perm]=bss_eval_sources(out(tblock_sdr, :)', target_interf_echo(tblock_sdr, 1:2)');
    sdri_wpe_nlms_bss(rptcount)=sdr(1)-sdr_input;
    siri_wpe_nlms_bss(rptcount)=sir(1)-sir_input;
    sisdri_wpe_nlms_bss(rptcount)=sisdr(target_interf_echo(tblock_sdr, 1), out(tblock_sdr, :))-sisdr_input;
    siiri_wpe_nlms_bss(rptcount)=10 * log10 (sum(out(tblock_siir_y, perm(1)).^2) / sum(out(tblock_siir_x, perm(1)).^2))-siir_input;
    sieri_wpe_nlms_bss(rptcount)=10 * log10 (sum(out(tblock_sier_y, perm(1)).^2) / sum(out(tblock_sier_x, perm(1)).^2))-sier_input;

    out=audioread([savedir '/02_out_nlms_wpe_bss_', num2str(rptcount), '.wav']);
    [sdr, sir, sar, perm]=bss_eval_sources(out(tblock_sdr, :)', target_interf_echo(tblock_sdr, 1:2)');
    sdri_nlms_wpe_bss(rptcount)=sdr(1)-sdr_input;
    siri_nlms_wpe_bss(rptcount)=sir(1)-sir_input;
    sisdri_nlms_wpe_bss(rptcount)=sisdr(target_interf_echo(tblock_sdr, 1), out(tblock_sdr, :))-sisdr_input;
    siiri_nlms_wpe_bss(rptcount)=10 * log10 (sum(out(tblock_siir_y, perm(1)).^2) / sum(out(tblock_siir_x, perm(1)).^2))-siir_input;
    sieri_nlms_wpe_bss(rptcount)=10 * log10 (sum(out(tblock_sier_y, perm(1)).^2) / sum(out(tblock_sier_x, perm(1)).^2))-sier_input;
    
    
    out=audioread([savedir '/03_out_dr_aec_bss_', num2str(rptcount), '.wav']);
    [sdr, sir, sar, perm]=bss_eval_sources(out(tblock_sdr, :)', target_interf_echo(tblock_sdr, 1:2)');
    sdri_dr_aec_bss(rptcount)=sdr(1)-sdr_input;
    siri_dr_aec_bss(rptcount)=sir(1)-sir_input;
    sisdri_dr_aec_bss(rptcount)=sisdr(target_interf_echo(tblock_sdr, 1), out(tblock_sdr, :))-sisdr_input;
    siiri_dr_aec_bss(rptcount)=10 * log10 (sum(out(tblock_siir_y, perm(1)).^2) / sum(out(tblock_siir_x, perm(1)).^2))-siir_input;
    sieri_dr_aec_bss(rptcount)=10 * log10 (sum(out(tblock_sier_y, perm(1)).^2) / sum(out(tblock_sier_x, perm(1)).^2))-sier_input;
    
    out=audioread([savedir '/04_out_aec_dr_bss_', num2str(rptcount), '.wav']);
    [sdr, sir, sar, perm]=bss_eval_sources(out(tblock_sdr, :)', target_interf_echo(tblock_sdr, 1:2)');
    sdri_aec_dr_bss(rptcount)=sdr(1)-sdr_input;
    siri_aec_dr_bss(rptcount)=sir(1)-sir_input;
    sisdri_aec_dr_bss(rptcount)=sisdr(target_interf_echo(tblock_sdr, 1), out(tblock_sdr, :))-sisdr_input;
    siiri_aec_dr_bss(rptcount)=10 * log10 (sum(out(tblock_siir_y, perm(1)).^2) / sum(out(tblock_siir_x, perm(1)).^2))-siir_input;
    sieri_aec_dr_bss(rptcount)=10 * log10 (sum(out(tblock_sier_y, perm(1)).^2) / sum(out(tblock_sier_x, perm(1)).^2))-sier_input;
    
    out=audioread([savedir '/05_out_draec_bss_', num2str(rptcount), '.wav']);
    [sdr, sir, sar, perm]=bss_eval_sources(out(tblock_sdr, :)', target_interf_echo(tblock_sdr, 1:2)');
    sdri_draec_bss(rptcount)=sdr(1)-sdr_input;
    siri_draec_bss(rptcount)=sir(1)-sir_input;
    sisdri_draec_bss(rptcount)=sisdr(target_interf_echo(tblock_sdr, 1), out(tblock_sdr, :))-sisdr_input;
    siiri_draec_bss(rptcount)=10 * log10 (sum(out(tblock_siir_y, perm(1)).^2) / sum(out(tblock_siir_x, perm(1)).^2))-siir_input;
    sieri_draec_bss(rptcount)=10 * log10 (sum(out(tblock_sier_y, perm(1)).^2) / sum(out(tblock_sier_x, perm(1)).^2))-sier_input;
    
    out=audioread([savedir '/06_out_joint_ss_', num2str(rptcount), '.wav']);
    [sdr, sir, sar, perm]=bss_eval_sources(out(tblock_sdr, :)', target_interf_echo(tblock_sdr, 1:2)');
    sdri_joint_ss(rptcount)=sdr(1)-sdr_input;
    siri_joint_ss(rptcount)=sir(1)-sir_input;
    sisdri_joint_ss(rptcount)=sisdr(target_interf_echo(tblock_sdr, 1), out(tblock_sdr, :))-sisdr_input;
    siiri_joint_ss(rptcount)=10 * log10 (sum(out(tblock_siir_y, perm(1)).^2) / sum(out(tblock_siir_x, perm(1)).^2))-siir_input;
    sieri_joint_ss(rptcount)=10 * log10 (sum(out(tblock_sier_y, perm(1)).^2) / sum(out(tblock_sier_x, perm(1)).^2))-sier_input;
end

save([savedir '_sdr.mat'], 'sdri*');
save([savedir '_sir.mat'], 'sdri*');
save([savedir '_sisdr.mat'], 'sisdri*');
save([savedir '_siiri.mat'], 'siiri*');
save([savedir '_sieri.mat'], 'sieri*');

if 0
figure;
boxplot([sdri_wpe_nlms_bss, sdri_nlms_wpe_bss, sdri_dr_aec_bss, sdri_aec_dr_bss, sdri_draec_bss, sdri_joint_ss],...
    'labels', {'WPE-NLMS-BSS', 'NLMS-WPE-BSS', 'DR-AEC-BSS', 'AEC-DR-BSS', 'DRAEC-BSS', 'Joint-SS'});
ylabel('SDR Improvement');
title(savedir,'Interpreter','none')

figure;
boxplot([siiri_wpe_nlms_bss, siiri_nlms_wpe_bss, siiri_dr_aec_bss, siiri_aec_dr_bss, siiri_draec_bss, siiri_joint_ss],...
    'labels', {'WPE-NLMS-BSS', 'NLMS-WPE-BSS', 'DR-AEC-BSS', 'AEC-DR-BSS', 'DRAEC-BSS', 'Joint-SS'});
ylabel('SIIR Improvement');
title(savedir,'Interpreter','none')

figure;
boxplot([sieri_wpe_nlms_bss, sieri_nlms_wpe_bss, sieri_dr_aec_bss, sieri_aec_dr_bss, sieri_draec_bss, sieri_joint_ss],...
    'labels', {'WPE-NLMS-BSS', 'NLMS-WPE-BSS', 'DR-AEC-BSS', 'AEC-DR-BSS', 'DRAEC-BSS', 'Joint-SS'});
ylabel('SIER Improvement');
title(savedir,'Interpreter','none')
end
