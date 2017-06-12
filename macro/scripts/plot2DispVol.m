clear; close all; clc;

%% 
macro=8;
directory = ['neuron_emb_',num2str(macro),'_del100_rot0_FT_singlescale_b_n1'];
% figure out number of lines ====
fid = fopen([directory,'/disps.log'],'r');
allText = textscan(fid,'%d%d%f%f%f','HeaderLines',1,'Delimiter',',');
stp = length(unique(allText{1}))-1;
disp = allText{3};
strain = (disp(1:2:end) - disp(2:2:end))/300;
fclose(fid);
% ===============================
num_tags = 15;
tag_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14];
% stp+1: ldstp counts from 0.
% stp+2: additional step due to initial volume written out at ldstp 0.
voldata=dlmread([directory,'/vols.log'],',',[1 0 num_tags*(stp+2) 2]);
count = 1;
volume_init = 0;
for line = 1:num_tags
    if (any(voldata(line,2) == tag_list))
        volume_init = volume_init + voldata(line,3);
    end
end

volume = zeros(1,stp+1);
for line = num_tags+1:num_tags*(stp+2)
    count = voldata(line,1) + 1;
    if (any(voldata(line,2) == tag_list))
        volume(count) = volume(count) + voldata(line,3);
    end
end

plot(strain,(volume)/volume_init,'o')
ylabel('$$\V/V_0$$','interpreter','latex')
xlabel('$$(l-l_0)/l_0$$ (bulk strain)','interpreter','latex')