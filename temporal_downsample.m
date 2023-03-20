clear all
close all
clc

for  i=1:7
cell_no = i;
cell_filename = append('outputs/long_bicarb_sim_cell_', int2str(cell_no),'_VPLC0.0015.mat');
cell_filename1 = append('downsampled_outputs/long_bicarb_sim_cell_', int2str(cell_no),'_VPLC0.0015.mat');
load(cell_filename, 'sol', 'tim');
time_samples = tim(1:2:end);
ca_solutions = sol(:,1:2:end);

save(cell_filename1,'time_samples','ca_solutions')
end
%  load(cell_filename1)
%  plot(time_samples, ca_solutions(1,:))
%  xlim([0,50])