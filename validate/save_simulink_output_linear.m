clear
clc

name = 'project1_bonus_rlc_circuit';
sim_file = strcat(name, '.slx');
out = sim(sim_file, 0.2);

t_out = out.tout';
n = length(t_out);
voltages_out = zeros(3,n);
for i = 1:3
    voltages_out(i,:) = out.yout{i}.Values.Data;
end

currents_out = zeros(3,n);
for i = 1:3
    currents_out(i,:) = out.yout{i+3}.Values.Data;
end

out_name = strcat(name, '_simulink_output.mat');

save(out_name, 't_out', 'voltages_out', 'currents_out');