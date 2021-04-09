%
% Test how to use wpe.
%
clear;

cfgs='settings/local.m';
x=audioread('../output/input_data_1.wav');
x(:, 3)=[];
y = wpe(x, cfgs); 

audiowrite('test_out.wav', y, 16000);
