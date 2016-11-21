clear all; close all; clc;

ldstp = 36;
micro_iterations = 22;
num_rank = 6;

f1_store = zeros(num_rank,micro_iterations);
f2_store = f1_store;

filename = 'MicroFiberNetwork_rank';
dir = '/Users/vicchan/SCOREC/Examine_Length_Ratio/vor_10seed_2sup/stress0_07_linear_E4_32/fiber_networks_per_macro_iter/';

for rank_idx = 1:num_rank
    for file_idx = 1:micro_iterations
        [f1,f2] = calculate_orientation([dir,filename,num2str(rank_idx-1),'_ldstp_',num2str(ldstp),'_',num2str(file_idx),'.txt']);
        f1_store(rank_idx,file_idx) = f1;
        f2_store(rank_idx,file_idx) = f2;
    end
end
iterations = linspace(1,micro_iterations,micro_iterations);

figure
plot(iterations,f1_store(1,:),'-o')
hold on
plot(iterations,f1_store(2,:),'-xr')
plot(iterations,f1_store(3,:),'-sqg')
plot(iterations,f1_store(4,:),'-^m')
plot(iterations,f1_store(5,:),'-*y')
plot(iterations,f1_store(6,:),'-+k')
legend('rank 1','rank 2','rank 3','rank 4','rank 5','rank 6','location','eastoutside')
xlabel('macro iteration number')
ylabel('$$f$$','interpreter','latex')
title(['main fibers, load step = ',num2str(ldstp)])
name = ['main_fiber_orientation_micro_iter_ldstp',num2str(ldstp),'.eps'];
saveas(gca,name,'epsc2')
eps2pdf(name)

figure
plot(iterations,f2_store(1,:),'-o')
hold on
plot(iterations,f2_store(2,:),'-xr')
plot(iterations,f2_store(3,:),'-sqg')
plot(iterations,f2_store(4,:),'-^m')
plot(iterations,f2_store(5,:),'-*y')
plot(iterations,f2_store(6,:),'-+k')
legend('rank 1','rank 2','rank 3','rank 4','rank 5','rank 6','location','eastoutside')
xlabel('macro iteration number')
ylabel('$$f$$','interpreter','latex')
title(['support fibers, load step = ',num2str(ldstp)])
name = ['support_fiber_orientation_macro_iter_ldstp',num2str(ldstp),'.eps'];
saveas(gca,name,'epsc2')
eps2pdf(name)