clear all; close all; clc;

num_files = 40; num_rank = 6;
f1_store = zeros(num_rank,num_files);
f2_store = f1_store;
max_stress = 0.07;

filename = 'FiberNetwork_rank';
dir = '/Users/vicchan/SCOREC/Examine_Length_Ratio/vor_10seed_2sup/stress0_07_linear_E4_32/fiber_networks_per_ldstp/';

for rank_idx = 1:num_rank
    for file_idx = 1:num_files
        [f1,f2] = calculate_orientation([dir,filename,num2str(rank_idx-1),'_',num2str(file_idx),'.txt']);
        f1_store(rank_idx,file_idx) = f1;
        f2_store(rank_idx,file_idx) = f2;
    end
end
stress = linspace(max_stress/num_files,max_stress,num_files);

figure
plot(stress,f1_store(1,:),'-o')
hold on
plot(stress,f1_store(2,:),'-xr')
plot(stress,f1_store(3,:),'-sqg')
plot(stress,f1_store(4,:),'-^m')
plot(stress,f1_store(5,:),'-*y')
plot(stress,f1_store(6,:),'-+k')
legend('rank 1','rank 2','rank 3','rank 4','rank 5','rank 6','location','best')
xlabel('$$\sigma_{11}$$','interpreter','latex')
ylabel('$$f$$','interpreter','latex')
title('main fibers')
name = 'main_fiber_orientation_ldstp.eps';
saveas(gca,name,'epsc2')
eps2pdf(name)

figure
plot(stress,f2_store(1,:),'-o')
hold on
plot(stress,f2_store(2,:),'-xr')
plot(stress,f2_store(3,:),'-sqg')
plot(stress,f2_store(4,:),'-^m')
plot(stress,f2_store(5,:),'-*y')
plot(stress,f2_store(6,:),'-+k')
legend('rank 1','rank 2','rank 3','rank 4','rank 5','rank 6','location','best')
xlabel('$$\sigma_{11}$$','interpreter','latex')
ylabel('$$f$$','interpreter','latex')
title('support fibers')
name = 'support_fiber_orientation_ldstp.eps';
saveas(gca,name,'epsc2')
eps2pdf(name)