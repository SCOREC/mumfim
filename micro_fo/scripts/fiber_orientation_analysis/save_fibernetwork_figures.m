clear all; close all; clc;


num_ldstp = 40;

filename = 'FiberNetwork_rank';
dir = '/Users/vicchan/SCOREC/N1P144/vor_10seed_2sup/stress0_07_linear_E4_32/fiber_networks_per_ldstp/';
max_stress = 0.07;
stress = linspace(max_stress/num_ldstp,max_stress,num_ldstp);

for rank_idx = 1:6
    close all
for ldstp_idx = 1:num_ldstp
    plot_supported_networks([dir,filename,num2str(rank_idx-1)],ldstp_idx,ldstp_idx,{[0,0,1],[1,0,0]})
    title(['$$\sigma_{11} =$$ ',num2str(stress(ldstp_idx))],'interpreter','latex') 
    name = ['rank',num2str(rank_idx),'/fiber_network_',num2str(rank_idx),'_',num2str(ldstp_idx),'.jpg'];
    saveas(gca,name)
end
end
