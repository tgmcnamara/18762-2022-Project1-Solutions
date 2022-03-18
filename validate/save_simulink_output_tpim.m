clear
clc

im_data = load('IM1.mat');
t_out = im_data.motor_out.vds.Time';
out_signals = zeros(8,length(t_out));
out_signals(1,:) = im_data.motor_out.ids.data';
out_signals(2,:) = im_data.motor_out.iqs.data';
out_signals(3,:) = im_data.motor_out.idr.data';
out_signals(4,:) = im_data.motor_out.iqr.data';
out_signals(5,:) = im_data.motor_out.w.data';
out_signals(6,:) = im_data.motor_out.Te.data';
out_signals(7,:) = im_data.motor_out.vds.data';
out_signals(8,:) = im_data.motor_out.vqs.data';

save('im_simple_circuit_signals.mat', 't_out', 'out_signals');

im_data = load('IM2.mat');
t_out = im_data.motor_out.vds.Time';
out_signals = zeros(6,length(t_out));
out_signals(1,:) = im_data.motor_out.ids.data';
out_signals(2,:) = im_data.motor_out.iqs.data';
out_signals(3,:) = im_data.motor_out.idr.data';
out_signals(4,:) = im_data.motor_out.iqr.data';
out_signals(5,:) = im_data.motor_out.w.data';
out_signals(6,:) = im_data.motor_out.Te.data';
out_signals(7,:) = im_data.motor_out.vds.data';
out_signals(8,:) = im_data.motor_out.vqs.data';
save('im_switch_circuit_signals.mat', 't_out', 'out_signals');

