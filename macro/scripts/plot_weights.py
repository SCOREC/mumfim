import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import math

"""plot_weights.py is a python script for plotting micro_fo_weight.#.log files.
1st argument is the name of the file (without .#.log).
2nd argument is the number of microprocessors.
3rd argument is the number of load steps.
4th argument is the output name.
"""
dirname = sys.argv[1]
num_files = int(sys.argv[2]) #(equal to number of processors: 0,...,N-1)
load_steps = int(sys.argv[3])
output_name = (sys.argv[4])
 
# obtain filename from commandline argument:sys.argv[1].
#filename=os.path.abspath('')+'/'+dirname
filename=dirname

# Generate list of filenames to be read in.
fnames = []
for ii in range(0,num_files):
    fnames.append(filename+"."+str(ii)+".log")

data = []
for files in fnames:
    dataIn = np.genfromtxt(files,delimiter=',',skip_header=1)
    data.append(dataIn)

data=np.array(data)

## Bar Plot
if (load_steps<4):
    plt_cols = math.ceil(math.sqrt(float(load_steps)))+1
    plt_rows = math.floor(math.sqrt(float(load_steps)))
else:
    plt_cols = math.ceil(math.sqrt(float(load_steps)))
    plt_rows = math.floor(math.sqrt(float(load_steps)))

print plt_cols
print plt_rows

# Specify two figures
fig=plt.figure(1)
fig2=plt.figure(2)

plt_height=1.0
plt_width=0.30
fig.subplots_adjust(hspace=plt_height,wspace=plt_width)
fig2.subplots_adjust(hspace=plt_height,wspace=plt_width)


avg_iter=[]
for ll in range(0,load_steps):
    ax = fig.add_subplot(plt_cols,plt_rows,ll+1)
    ax2 = fig2.add_subplot(plt_cols,plt_rows,ll+1)
    proc_axis=[]
    iter_axis=[]
    gt_avg=[]
    lt_avg=[]
    gt_iter=0
    lt_iter=0
    total_iter=0
    avg_iter.append(data[:,ll,1].sum()/num_files)

    for ii in range(0,num_files): # Loop through number of processors
        proc_axis.append(ii)
        iter_axis.append(data[ii,ll,1])
        #total_iter=total_iter+iter_axis[ii]

        if (iter_axis[ii]>avg_iter[ll]):
            gt_avg.append(iter_axis[ii]-avg_iter[ll])
            lt_avg.append(0)
            gt_iter = gt_iter + gt_avg[ii]
        else:
            gt_avg.append(0)
            lt_avg.append(iter_axis[ii]-avg_iter[ll])
            lt_iter = lt_iter + abs(lt_avg[ii])

    # calculate average deviation
    avg_deviation = (gt_iter+lt_iter)/num_files/avg_iter[ll]*100
    max_deviation = max([max(gt_avg),max(lt_avg)])/avg_iter[ll]*100

# plot bar plot in figure 1 for each load step.
    ax.bar(proc_axis,iter_axis,color='k')
    #plot average
    ax.axhline(y=avg_iter[ll],color='r')
    #ax.set_ylabel('iterations')
    #ax.set_xlabel('proc #')
    ax.set_title('load step '+str(ll+1)+'\n'\
                 +'avg iter:'+str(avg_iter[ll])[0:6]+'\n'\
                 +'avg dev:'+str(avg_deviation)[0:4]+'%\n'\
                 +'max dev:'+str(max_deviation)[0:4]+'%',fontsize=10)
    
    ax.set_xlim([0,num_files])

# label figure 1 of each load step with % iterations greater than and less than the average, which is stored in avg_iter.
#    ax.annotate('gt % ='+str(gt_iter/total_iter*100)[0:2]+'%',xy=(0.05,0.85),xycoords='axes fraction',size=10)
#    ax.annotate('lt % ='+str(lt_iter/total_iter*100)[0:2]+'%',xy=(0.05,0.70),xycoords='axes fraction',size=10)

# plot bar plot in figure 2 for each load step.
    ax2.bar(proc_axis,gt_avg,color='b')
    ax2.bar(proc_axis,lt_avg,color='r')
    ax2.set_title('load step '+str(ll+1)+'\n'\
                 +'avg iter:'+str(avg_iter[ll])[0:6]+'\n'\
                 +'avg dev:'+str(avg_deviation)[0:4]+'%\n'\
                 +'max dev:'+str(max_deviation)[0:4]+'%',fontsize=10)
    ax2.set_xlim([0,num_files])

# label figure 2 of each load step with % iterations greater than and less than the average, which is stored in avg_iter.
#    ax2.annotate('gt % ='+str(gt_iter/total_iter*100)[0:2]+'%',xy=(0.05,0.85),xycoords='axes fraction',size=10)
#    ax2.annotate('lt % ='+str(lt_iter/total_iter*100)[0:2]+'%',xy=(0.05,0.70),xycoords='axes fraction',size=10)

# label figure 1
fig.text(0.04,0.5,"Newton Iterations",rotation="vertical",va="center")
fig.text(0.5,0.04,"Processor Number",ha="center")

# label figure 2
fig2.text(0.04,0.5,"Newton Iterations",rotation="vertical",va="center")
fig2.text(0.5,0.04,"Processor Number",ha="center")

plt.figure(1)
plt.savefig(output_name+'_1.pdf',bbox_inches='tight')

plt.figure(2)
plt.savefig(output_name+'_2.pdf',bbox_inches='tight')
#plt.show()




