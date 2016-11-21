import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import math
import gantt_util as util

## Read files
dirname_macro = sys.argv[1]
num_files_macro = int(sys.argv[2])
dirname_micro = sys.argv[3]
num_files_micro = int(sys.argv[4])
method = sys.argv[5]
output_file = sys.argv[6]

data_macro, proc_macro, f_macro, f_macro_no = util.fileIO(dirname_macro,num_files_macro)
data_micro, proc_micro, f_micro, f_micro_no = util.fileIO(dirname_micro,num_files_micro)
## calculate % macro/micro files used
percent_micro = float(f_micro)/(float(f_micro)+float(f_micro_no))*100
percent_macro = float(f_macro)/(float(f_macro)+float(f_macro_no))*100

# calculate average time of micro-scale POST_SEND
micro_PS_avg = util.calculate_micro_PS_avg(data_micro,proc_micro)

# replace macro-scale time before #,# with proper entries of micro_PS_avg
data_macro = util.correct_macro_PS_time(data_macro,proc_macro,micro_PS_avg)

# set transtion between iterations to be IDLE in macro-scale.
data_macro = util.set_macro_it_transiton_IDLE(data_macro,proc_macro)

# reformat data_micro and data_macro to 'start', 'end', and 'status' format.
micro_start, micro_end, micro_stat = util.format_data(data_micro,proc_micro)
macro_start, macro_end, macro_stat = util.format_data(data_macro,proc_macro)

# determine average or max_min times.
stime_micro_IDLE, etime_micro_IDLE = util.extract_times(micro_start, micro_end, micro_stat, proc_micro, ' IDLE', method)
stime_macro_IDLE, etime_macro_IDLE = util.extract_times(macro_start, macro_end, macro_stat, proc_macro, ' IDLE', method)

stime_micro_LB, etime_micro_LB = util.extract_times(micro_start, micro_end, micro_stat, proc_micro, ' LB', method)
stime_macro_LB, etime_macro_LB = util.extract_times(macro_start, macro_end, macro_stat, proc_macro, ' LB', method)

max_macro_time = np.amax(macro_end[:][:])
max_micro_time = np.amax(micro_end[:][:])

macro_raw_IDLE = 0
micro_raw_IDLE = 0
macro_raw_LB = 0
micro_raw_LB = 0
## macro_IDLE_raw percentage
for ii in range(0,len(stime_macro_IDLE)):
    macro_raw_IDLE = macro_raw_IDLE + etime_macro_IDLE[ii] - stime_macro_IDLE[ii]

## micro_IDLE_raw pecentage
for ii in range(0,len(stime_micro_IDLE)):
    micro_raw_IDLE = micro_raw_IDLE + etime_micro_IDLE[ii] - stime_micro_IDLE[ii]

macro_raw_IDLE_pct = macro_raw_IDLE/max_macro_time*100
micro_raw_IDLE_pct = micro_raw_IDLE/max_micro_time*100

## macro_LB_raw percentage
for ii in range(0,len(stime_macro_LB)):
    macro_raw_LB = macro_raw_LB + etime_macro_LB[ii] - stime_macro_LB[ii]

## micro_LB_raw pecentage
for ii in range(0,len(stime_micro_LB)):
    micro_raw_LB = micro_raw_LB + etime_micro_LB[ii] - stime_micro_LB[ii]

macro_raw_LB_pct = macro_raw_LB/max_macro_time*100
micro_raw_LB_pct = micro_raw_LB/max_micro_time*100

weighted_IDLE_pct = (num_files_macro*macro_raw_IDLE+num_files_micro*micro_raw_IDLE)/(num_files_macro*max_macro_time+num_files_micro*max_micro_time)*100
 
## =========
# Plot data
## =========
fig=plt.figure()
ax=fig.add_subplot(111)

## Plot ACTIVE
ax.barh(1,max_macro_time,color='g',alpha=0.7,align='center')
ax.barh(0,max_micro_time,color='g',alpha=0.7,align='center')

## Plot IDLE
for ii in range(0,len(stime_macro_IDLE)):
    ax.barh(1,etime_macro_IDLE[ii]-stime_macro_IDLE[ii],left=stime_macro_IDLE[ii],color='r',align='center')

for ii in range(0,len(stime_micro_IDLE)):
    ax.barh(0,etime_micro_IDLE[ii]-stime_micro_IDLE[ii],left=stime_micro_IDLE[ii],color='r',align='center')

## Plot LB
for ii in range(0,len(stime_macro_LB)):
    ax.barh(1,etime_macro_LB[ii]-stime_macro_LB[ii],left=stime_macro_LB[ii],color='b',align='center')
for ii in range(0,len(stime_micro_LB)):
    ax.barh(0,etime_micro_LB[ii]-stime_micro_LB[ii],left=stime_micro_LB[ii],color='b',align='center')

## label x and y labels.
ax.set(yticks=[0,1],yticklabels=['micro','macro'])
ax.set_xlabel('Time (sec)')

## label title.
# macro_micro size and method use
title_method = sys.argv[6]+'_'+sys.argv[4]+', method:'+method 
title_fidelity = str(percent_micro)[0:4]+'%micro files and '+str(percent_macro)[0:4]+'%macro files'
title_raw_stats_IDLE = '%micro IDLE (raw):'+str(micro_raw_IDLE_pct)[0:4]+', %macro IDLE (raw):'+str(macro_raw_IDLE_pct)[0:4]
title_raw_stats_LB = ' %micro LB (raw):'+str(micro_raw_LB_pct)[0:4]+', %macro LB (raw):'+str(macro_raw_LB_pct)[0:4]
title_weighted_IDLE = ' %IDLE (weighted):'+str(weighted_IDLE_pct)[0:4]

ax.set_title(title_method+'\n'+title_fidelity+'\n'\
            +title_raw_stats_IDLE+'\n'+title_raw_stats_LB+'\n'\
            +title_weighted_IDLE)
fig.savefig(output_file+'.pdf',bbox_inches='tight')
#plt.show()



