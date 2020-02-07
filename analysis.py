import matplotlib.pyplot as plt
import pandas
import copy
import pdb
import pandas as pd
from utils import *

# plt.plot([7813, 7813], [0., 0.9], "r:")         # Not shown
# plt.plot([-50000, 7813], [0.9, 0.9], "r:")      # Not shown
# plt.plot([-50000, 7813], [0.4368, 0.4368], "r:")# Not shown
# plt.plot([7813], [0.9], "ro")                   # Not shown
# plt.plot([7813], [0.4368], "ro")                # Not shown

def save_fig(fig_id, tight_layout=True, fig_extension="png", resolution=300):
    path = os.path.join(PLOTS_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)

    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)

"""
analysis ideas:

comparison family graphs - which kind of graphs which is better

variance ? with a,b,c

comparison how objective changes in each iter in nr / bs - what kind of path tehy follow ,  why does it start 
very close to the answer / does it change from family to family where it starts?

are the performance of the algorithm change with lambda

base vs reliable solution - gains?

"""

# get the subplot from handsonml notebooks"
def graph_family(family_name, cplex_times, nr_times):
    node_num = ['2^10', '2^11', '2^12', '2^13', '2^14']
    plt.plot(node_num, cplex_times, 'b:o', linewidth=2, label='CPLEX')
    plt.plot(node_num, nr_times, 'g-+', linewidth=2, label='NR')

    plt.ylabel('Running time', fontsize=16)
    plt.xlabel('Number of nodes', fontsize=16)
    # plt.axis([xs,xe,ys,ye])
    plt.grid(True)
    plt.title(family_name, fontsize=16)
    plt.legend(loc='lower right', fontsize=16)
    save_fig(family_name)



experiment = 'graph_families'

cplex_times = []
nr_times = []

bases = ['netgen', 'goto']
tail = 'a'
exponents = (np.arange(10, 15)).astype(str)
types = ['_lo_8_', '_sr_', '_8_', ]

for base in bases:
    for atype in types:
        cplex_times = []
        nr_times = []
        for exponent in exponents:
            filename = experiment + '/' + base + '/' + base + atype + exponent + tail
            data_dic = load_run(filename)
            cplex_times.append(data_dic['solver_elapsed'])
            nr_times.append(data_dic['nr_elapsed'])
        family_type = atype.strip('_')
        family_name = base.upper() + '-' + family_type + ' family'
        graph_family(family_name, cplex_times, nr_times)

pdb.set_trace()











def get_level_time(gap_list, time_list):
    level_time = {}
    levels = [1e3,1e0,1e-3,1e-6]
    for level in levels:
        level_time[level] = 0

        for j in range(10):
            try:
                ind = np.where(np.array(gap_list[j])<level)[0][0]
                level_time[level] += time_list[j][ind]
            except:
                pdb.set_trace()
        level_time[level] = level_time[level]/10.0

    sorted_d = sorted(level_time.items(), key=lambda x: x[0])
    gap, time = zip(*sorted_d)
    return gap, time

def read_log(logfilename):
    f = open(logfilename + '_processed.txt', 'w')
    ff = open(logfilename + '.txt', 'r+')
    w_to_processed = False

    for line in ff.readlines():
        if line.find('ITE PFEAS') >= 0:
            w_to_processed = True
        if line.find('Optimizer terminated') >=0:
            w_to_processed = False
        if w_to_processed:
            f.write(line)
    f.close()

    log = pandas.read_table(logfilename + '_processed.txt', delim_whitespace=True, sep='\t', lineterminator='\n')
    return log

def find_eq_obj_ind(alg, best_obj, descent_to):
    len_alg = len(alg)
    for i in range(len_alg):
        if 100*((alg[i]/benchmark)-1) <= descent_to:
            return i+1




lams = np.linspace(0, 5, num=50)
n = 250
arcs = 1000
mu=50
var=50
d_top = 1
mucosts=[]
varcosts=[]
seed = 9
lam_costs = {}
iters = 0 
vallams = []
for lam in lams:
    if (iters%1==0 and iters<=10) or (iters%5==0 and iters>10) or lam==lams[-1] :
        lamstr = str(round(lam, 6))
        lamstr = lamstr.replace('.','')

        save_extension = str(mu) + str(var) + '_' + str(arcs) + '_' + str(seed) + str(d_top)

        experiment_name = 'varying_lams'
        mucost = load( 'mu_cost_' + lamstr + '_' + save_extension , n, experiment_name)
        varcost = load( 'var_cost_' + lamstr + '_' + save_extension, n, experiment_name)
        mucosts.append(mucost)
        varcosts.append(varcost)
        vallams.append(lam)
        lam_costs[lam] = (mucost, varcost)
    iters += 1

x=[]
y=[]
iters = 0

for lam in lams:
    if (iters%1==0 and iters<=10) or (iters%5==0 and iters>10) or lam==lams[-1] :
        y.append(lam_costs[lam])
        x.append(lam)
    iters +=1

fig, ax1 = plt.subplots()

# for xe, ye in zip(x, y):
#     plt.scatter([xe] * len(ye), ye)
ax2 = ax1.twinx()
pdb.set_trace()
ax1.scatter(vallams, mucosts)
ax2.scatter(vallams, varcosts)


# plt.xticks(list(np.arange(lams)))
def major_formatter(x, pos):
    return "[%.2f]" % x
import matplotlib
ax1.set_xticks([lams[0], lams[5], lams[10], lams[15], lams[20], lams[25], lams[30], lams[35], lams[40], lams[45], lams[-1]])
ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
# ax1.xaxis.set_ticks

lns1 = ax1.plot(vallams, mucosts, label='Mean')
lns2 = ax2.plot(vallams, varcosts, '--', label='Standard Deviation')


# added these three lines
lns = lns1+lns2
labs = [l.get_label() for l in lns]

box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width, box.height * 0.8])
ax2.set_position([box.x0, box.y0, box.width , box.height* 0.8])


ax1.legend(lns, labs, loc='center top', bbox_to_anchor=(1, 0.5))




ax1.set_xlabel('Reliability Parameter')
ax1.set_ylabel('Mean Cost')
ax2.set_ylabel('Standard Deviation Cost')


# plt.plot(lams, mucosts, label='Mean')
# plt.plot(lams, varcosts, '-', label='Standard Deviation')
# plt.legend(loc='upper right')
# plt.xlabel('Value of Reliability')
# plt.ylabel('Mean Cost')
# plt.ylabel('Mean Cost')

# ax2.set_ylabel('Y2 data', color='b')
plt.savefig('plots/' + 'meanvsstd')
# plt.show()
#


lams = np.linspace(0, 5, num=50)
n = 250
arcs = 2000
mu=50
var=50
d_top = 1
mucosts=[]
varcosts=[]
seed = 0
lam_costs = {}
iters = 0 
vallams = []
for lam in lams:
    if (iters%1==0 and iters<=10) or (iters%5==0 and iters>10) or lam==lams[-1] :
        lamstr = str(round(lam, 6))
        lamstr = lamstr.replace('.','')

        save_extension = str(mu) + str(var) + '_' + str(arcs) + '_' + str(seed) + str(d_top)

        experiment_name = 'varying_lams'
        mucost = load( 'mu_cost_' + lamstr + '_' + save_extension , n, experiment_name)
        varcost = load( 'var_cost_' + lamstr + '_' + save_extension, n, experiment_name)
        mucosts.append(mucost)
        varcosts.append(varcost)
        vallams.append(lam)
        lam_costs[lam] = (mucost, varcost)
    iters += 1

x=[]
y=[]
iters = 0

for lam in lams:
    if (iters%1==0 and iters<=10) or (iters%5==0 and iters>10) or lam==lams[-1] :
        y.append(lam_costs[lam])
        x.append(lam)
    iters +=1

fig, ax1 = plt.subplots()

# for xe, ye in zip(x, y):
#     plt.scatter([xe] * len(ye), ye)
ax2 = ax1.twinx()
pdb.set_trace()
ax1.scatter(vallams, mucosts)
ax2.scatter(vallams, varcosts)


def major_formatter(x, pos):
    return "[%.2f]" % x
import matplotlib
ax1.set_xticks([lams[0], lams[5], lams[10], lams[15], lams[20], lams[25], lams[30], lams[35], lams[40], lams[45], lams[-1]])
ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
# ax1.xaxis.set_ticks

lns1 = ax1.plot(vallams, mucosts, label='Mean')
lns2 = ax2.plot(vallams, varcosts, '--', label='Standard Deviation')


# added these three lines
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=1)


ax1.set_xlabel('Value of Reliability')
ax1.set_ylabel('Mean Cost')
ax2.set_ylabel('Standard Deviation Cost')


# plt.plot(lams, mucosts, label='Mean')
# plt.plot(lams, varcosts, '-', label='Standard Deviation')
# plt.legend(loc='upper right')
# plt.xlabel('Value of Reliability')
# plt.ylabel('Mean Cost')
# plt.ylabel('Mean Cost')

# ax2.set_ylabel('Y2 data', color='b')
plt.savefig('plots/' + 'meanvsstd_8degree')
# plt.show()
#







#
nodes = [5000, 10000, 25000, 50000]
### Experiment Analysis ###
### Experiment 1 ###
experiment_name = 'performance'
lam = 0.7
mu = 50
d_top = 1
variances = [20]

seeds=8 * np.arange(10)

# plt.figure()

bx_nrtd = {}
bx_btd = {}
bx_nrbtd = {}
bx_msktd = {}

nr_time_dict = {}
bs_time_dict = {}
nrbs_time_dict = {}
msk_time_dict = {}

nr_gap_dict = {}
bs_gap_dict = {}
nrbs_gap_dict = {}
msk_gap_dict = {}

for n in nodes:
    i = 0
    arcs = n*4

    nr_time_dict[n] = {}
    bs_time_dict[n] = {}
    nrbs_time_dict[n] = {}
    msk_time_dict[n] = {}

    nr_gap_dict[n] = {}
    bs_gap_dict[n] = {}
    nrbs_gap_dict[n] = {}
    msk_gap_dict[n] = {}

    bx_nrtd[n] = {}
    bx_btd[n] = {}
    bx_nrbtd[n] = {}
    bx_msktd[n] = {}

    for var in variances:
        print(n, var)
        i += 1

        bx_nrtd[n][var] = []
        bx_btd[n][var] = []
        bx_nrbtd[n][var] = []
        bx_msktd[n][var] = []


        nr_gap_dict[n][var] = {}
        bs_gap_dict[n][var] = {}
        nrbs_gap_dict[n][var] = {}
        msk_gap_dict[n][var] = {}

        nr_time_dict[n][var] = {}
        bs_time_dict[n][var] = {}
        nrbs_time_dict[n][var] = {}
        msk_time_dict[n][var] = {}

        setup_times_cvx_l = []
        solve_times_nr_l = []
        solve_times_bs_l = []
        solve_times_alg3_l = []

        for seed in seeds:

            if n == 50000:
                if seed == 8:
                    seed = 16
                if seed == 9:
                    continue
            lamstr=str(lam)
            lamstr = lamstr.replace(".","-")
            if lam == 1:
                lamstr_c = '1-0'
            elif lam == 0.7:
                lamstr_c = '0-7'
            else:
                lamstr_c = '0-1'

            save_extension = lamstr + '_' + str(mu) + str(var) + '_' + str(arcs) + '_' + str(seed) + str(d_top)

            times_alg3 = np.array(load('times_alg3' + save_extension, n, experiment_name))
            times_nr = np.array(load('times_nr' + save_extension, n, experiment_name))
            times_bs = np.array(load('times_bs' + save_extension, n, experiment_name))

            logfilename = 'saved_runs/' + str(n) + '/' + experiment_name + lamstr_c + '_' + str(mu) + str(var) +  '_' + str(arcs) + '_' + str(seed)

            log = read_log(logfilename)
            overhead_time = np.array(load('mosek_overhead_time' + save_extension, n, experiment_name))
            log = log[log['PRSTATUS']>=0.997]
            log = log[log['PRSTATUS']<=1.003]
            cvx_objs = log['POBJ'].tolist()
            times_cvx = log['TIME'].tolist()
            times_cvx = np.array(times_cvx) + overhead_time

            bs_objs = np.array(load('bs_objs' + save_extension, n, experiment_name))
            nr_objs = np.array(load('nr_objs' + save_extension, n, experiment_name))
            alg3_objs = np.array(load('alg3_objs' + save_extension, n, experiment_name))

            bs_min = min(bs_objs)
            nr_min = min(nr_objs)
            alg3_min = min(alg3_objs)
            cvx_min = cvx_objs[-1]

            benchmark = min(bs_min, nr_min, alg3_min, cvx_min)
            descent_to = max(100*((nr_min/benchmark)-1), 100*((bs_min/benchmark)-1), 100*((alg3_min/benchmark)-1), 100*((cvx_min/benchmark)-1))

            print('min percentage diff: %', descent_to)
            # if descent_to > 1e-5:
            #     pdb.set_trace()
            # descent_to = min(descent_to, 1e-4)
            # descent_to=1e-4

            # bs_i = find_eq_obj_ind(bs_objs, benchmark, descent_to)
            # nr_i = find_eq_obj_ind(nr_objs, benchmark, descent_to)
            # alg3_i = find_eq_obj_ind(alg3_objs, benchmark, descent_to)

            # nr_objs = nr_objs[:nr_i]
            # bs_objs = bs_objs[:bs_i]
            # alg3_objs = alg3_objs[:alg3_i]

            # times_alg3 = times_alg3[:alg3_i]
            # times_nr = times_nr[:nr_i]
            # times_bs = times_bs[:bs_i]

            gap_bs = abs(100*(np.array(bs_objs)/benchmark -1))
            gap_nr = abs(100*(np.array(nr_objs)/benchmark -1))
            gap_alg3 = abs(100*(np.array(alg3_objs)/benchmark -1))
            gap_cvx = abs(100*(np.array(cvx_objs)/benchmark -1))
            setup_times_cvx_l.append(times_cvx[-1])
            solve_times_nr_l.append(times_nr[-1])
            solve_times_bs_l.append(times_bs[-1])
            solve_times_alg3_l.append(times_alg3[-1])

            nr_time_dict[n][var][seed] = times_nr
            nrbs_time_dict[n][var][seed] = times_alg3
            bs_time_dict[n][var][seed] = times_bs
            msk_time_dict[n][var][seed] = times_cvx

            nr_gap_dict[n][var][seed] = gap_nr
            nrbs_gap_dict[n][var][seed] = gap_alg3
            bs_gap_dict[n][var][seed] = gap_bs
            msk_gap_dict[n][var][seed] = gap_cvx

            nr_time_dict[n][var][seed] = times_nr[1:]
            nrbs_time_dict[n][var][seed] = times_alg3[1:]
            bs_time_dict[n][var][seed] = times_bs[1:]

            nr_gap_dict[n][var][seed] = gap_nr[1:]
            nrbs_gap_dict[n][var][seed] = gap_alg3[1:]
            bs_gap_dict[n][var][seed] = gap_bs[1:]

        bx_nrtd[n][var] = solve_times_nr_l
        bx_btd[n][var] = solve_times_bs_l
        bx_nrbtd[n][var] = solve_times_alg3_l
        bx_msktd[n][var] = setup_times_cvx_l

n = 5000
n2 = 10000  
n3 = 25000
n4 = 50000

print(tabulate([['NR', np.mean(bx_nrtd[n][var]), np.mean(bx_nrtd[n2][var]), np.mean(bx_nrtd[n3][var]), np.mean(bx_nrtd[n4][var])],['NRBS', np.mean(bx_nrbtd[n][var]), np.mean(bx_nrbtd[n2][var]), np.mean(bx_nrbtd[n3][var]), np.mean(bx_nrbtd[n4][var])],['BS', np.mean(bx_btd[n][var]), np.mean(bx_btd[n2][var]), np.mean(bx_btd[n3][var]), np.mean(bx_btd[n4][var])],['MSK', np.mean(bx_msktd[n][var]), np.mean(bx_msktd[n2][var]), np.mean(bx_msktd[n3][var]), np.mean(bx_msktd[n4][var])]], headers=['n=5000', 'n=100000', 'n=25000', 'n=50000']))


from scipy import stats
sample_size = 10
t_critical = stats.t.ppf(q = 0.9, df=sample_size)

yer_nr = []
yer_nrbs = []
yer_bs = []
yer_msk = []

nr_bars = []
nrbs_bars = []
bs_bars = []
msk_bars = []

for n in [10000, 25000, 50000]:
    nr_bars.append(np.mean(bx_nrtd[n][var]))
    sample_stdev = np.std(bx_nrtd[n][var])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_nr.append(t_critical * sigma)

    nrbs_bars.append(np.mean(bx_nrbtd[n][var]))
    sample_stdev = np.std(bx_nrbtd[n][var])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_nrbs.append(t_critical * sigma)

    bs_bars.append(np.mean(bx_btd[n][var]))
    sample_stdev = np.std(bx_btd[n][var])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_bs.append(t_critical * sigma)

    msk_bars.append(np.mean(bx_msktd[n][var]))
    sample_stdev = np.std(bx_msktd[n][var])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_msk.append(t_critical * sigma)


plt.figure()

barWidth = 0.2
rnr = np.arange(len(nr_bars))
rnrbs = [x + barWidth for x in rnr]
rbs = [x + 2*barWidth for x in rnr]
rmsk = [x + 3*barWidth for x in rnr]

plt.bar(rnr, nr_bars, width = barWidth, hatch='///', edgecolor = 'black', color='w', yerr=yer_nr, capsize=7, label='NR')
plt.bar(rnrbs, nrbs_bars, width = barWidth, hatch='\\\\\\', edgecolor = 'black', color='w',yerr=yer_nrbs, capsize=7, label='NR-BSC')
plt.bar(rbs, bs_bars, width = barWidth, hatch='xxx', edgecolor = 'black',color='w', yerr=yer_bs, capsize=7, label='BSC')
plt.bar(rmsk, msk_bars, width = barWidth, hatch='...', edgecolor = 'black', color='w',yerr=yer_msk, capsize=7, label='MOSEK')
plt.ylabel('Time(s)')
plt.xticks([(r + 2*barWidth) for r in range(len(nr_bars))], ['10000', '25000', '50000'])
plt.xlabel('# of Nodes')
plt.legend()
plt.tight_layout()
plt.savefig('plots/' + 'comparison')








nodes = [25000]
### Experiment Analysis ###
### Experiment 1 ###
experiment_name = 'varying_lams_2'
lams = [0.1, 0.5, 1, 5]
mu = 50
d_top = 1
variances = [50]

seeds=8 * np.arange(10)

# plt.figure()

bx_nrtd = {}
bx_btd = {}
bx_nrbtd = {}
bx_msktd = {}

nr_time_dict = {}
bs_time_dict = {}
nrbs_time_dict = {}
msk_time_dict = {}

nr_gap_dict = {}
bs_gap_dict = {}
nrbs_gap_dict = {}
msk_gap_dict = {}
var=50
for n in nodes:
    i = 0
    arcs = n*4

    nr_time_dict[n] = {}
    bs_time_dict[n] = {}
    nrbs_time_dict[n] = {}
    msk_time_dict[n] = {}

    nr_gap_dict[n] = {}
    bs_gap_dict[n] = {}
    nrbs_gap_dict[n] = {}
    msk_gap_dict[n] = {}

    bx_nrtd[n] = {}
    bx_btd[n] = {}
    bx_nrbtd[n] = {}
    bx_msktd[n] = {}

    for lam in lams:
        i += 1

        bx_nrtd[n][lam] = []
        bx_btd[n][lam] = []
        bx_nrbtd[n][lam] = []
        bx_msktd[n][lam] = []


        nr_gap_dict[n][lam] = {}
        bs_gap_dict[n][lam] = {}
        nrbs_gap_dict[n][lam] = {}
        msk_gap_dict[n][lam] = {}

        nr_time_dict[n][lam] = {}
        bs_time_dict[n][lam] = {}
        nrbs_time_dict[n][lam] = {}
        msk_time_dict[n][lam] = {}

        setup_times_cvx_l = []
        solve_times_nr_l = []
        solve_times_bs_l = []
        solve_times_alg3_l = []

        for seed in seeds:

            lamstr=str(lam)
            lamstr = lamstr.replace(".","-")

            if lam == 1:
                lamstr_c = '1-0'
            elif lam == 5:
                lamstr_c = '5-0'
            elif lam == 0.5:
                lamstr_c = '0-5'
            else:
                lamstr_c = '0-1'

            print(lam, lamstr_c)

            save_extension = lamstr_c + '_' + str(mu) + str(var) + '_' + str(arcs) + '_' + str(seed) + str(d_top)

            times_alg3 = np.array(load('times_alg3' + save_extension, n, experiment_name))
            times_nr = np.array(load('times_nr' + save_extension, n, experiment_name))
            times_bs = np.array(load('times_bs' + save_extension, n, experiment_name))
            if lam==5:
                lamstr_c = str(5)

                    
            logfilename = 'saved_runs/' + str(n) + '/' + experiment_name + lamstr_c + '_' + str(mu) + str(var) +  '_' + str(arcs) + '_' + str(seed)

            log = read_log(logfilename)
            overhead_time = np.array(load('mosek_overhead_time' + save_extension, n, experiment_name))
            log = log[log['PRSTATUS']>=0.997]
            log = log[log['PRSTATUS']<=1.003]
            cvx_objs = log['POBJ'].tolist()
            times_cvx = log['TIME'].tolist()
            times_cvx = np.array(times_cvx) + overhead_time

            bs_objs = np.array(load('bs_objs' + save_extension, n, experiment_name))
            nr_objs = np.array(load('nr_objs' + save_extension, n, experiment_name))
            alg3_objs = np.array(load('alg3_objs' + save_extension, n, experiment_name))

            bs_min = min(bs_objs)
            nr_min = min(nr_objs)
            alg3_min = min(alg3_objs)
            cvx_min = cvx_objs[-1]

            print(cvx_min, nr_min, bs_min, alg3_min)

            benchmark = min(bs_min, nr_min, alg3_min, cvx_min)
            descent_to = max(100*((nr_min/benchmark)-1), 100*((bs_min/benchmark)-1), 100*((alg3_min/benchmark)-1), 100*((cvx_min/benchmark)-1))

            print('min percentage diff: %', descent_to)
            # if descent_to > 1e-5:
            #     pdb.set_trace()
            # descent_to = min(descent_to, 1e-4)
            # descent_to=1e-4

            # bs_i = find_eq_obj_ind(bs_objs, benchmark, descent_to)
            # nr_i = find_eq_obj_ind(nr_objs, benchmark, descent_to)
            # alg3_i = find_eq_obj_ind(alg3_objs, benchmark, descent_to)

            # nr_objs = nr_objs[:nr_i]
            # bs_objs = bs_objs[:bs_i]
            # alg3_objs = alg3_objs[:alg3_i]

            # times_alg3 = times_alg3[:alg3_i]
            # times_nr = times_nr[:nr_i]
            # times_bs = times_bs[:bs_i]

            gap_bs = abs(100*(np.array(bs_objs)/benchmark -1))
            gap_nr = abs(100*(np.array(nr_objs)/benchmark -1))
            gap_alg3 = abs(100*(np.array(alg3_objs)/benchmark -1))
            gap_cvx = abs(100*(np.array(cvx_objs)/benchmark -1))

            setup_times_cvx_l.append(times_cvx[-1])
            solve_times_nr_l.append(times_nr[-1])
            solve_times_bs_l.append(times_bs[-1])
            solve_times_alg3_l.append(times_alg3[-1])

            nr_time_dict[n][lam][seed] = times_nr
            nrbs_time_dict[n][lam][seed] = times_alg3
            bs_time_dict[n][lam][seed] = times_bs
            msk_time_dict[n][lam][seed] = times_cvx

            nr_gap_dict[n][lam][seed] = gap_nr
            nrbs_gap_dict[n][lam][seed] = gap_alg3
            bs_gap_dict[n][lam][seed] = gap_bs
            msk_gap_dict[n][lam][seed] = gap_cvx

            nr_time_dict[n][lam][seed] = times_nr[1:]
            nrbs_time_dict[n][lam][seed] = times_alg3[1:]
            bs_time_dict[n][lam][seed] = times_bs[1:]

            nr_gap_dict[n][lam][seed] = gap_nr[1:]
            nrbs_gap_dict[n][lam][seed] = gap_alg3[1:]
            bs_gap_dict[n][lam][seed] = gap_bs[1:]

        bx_nrtd[n][lam] = solve_times_nr_l
        bx_btd[n][lam] = solve_times_bs_l
        bx_nrbtd[n][lam] = solve_times_alg3_l
        bx_msktd[n][lam] = setup_times_cvx_l

lam = 0.1
lam2 = 0.5  
lam3 = 1
lam4 = 5
pdb.set_trace()
# print(tabulate([['NR', np.mean(bx_nrtd[n][var]), np.mean(bx_nrtd[n2][var]), np.mean(bx_nrtd[n3][var]), np.mean(bx_nrtd[n4][var])],['NRBS', np.mean(bx_nrbtd[n][var]), np.mean(bx_nrbtd[n2][var]), np.mean(bx_nrbtd[n3][var]), np.mean(bx_nrbtd[n4][var])],['BS', np.mean(bx_btd[n][var]), np.mean(bx_btd[n2][var]), np.mean(bx_btd[n3][var]), np.mean(bx_btd[n4][var])],['MSK', np.mean(bx_msktd[n][var]), np.mean(bx_msktd[n2][var]), np.mean(bx_msktd[n3][var]), np.mean(bx_msktd[n4][var])]], headers=['n=5000', 'n=100000', 'n=25000', 'n=50000']))

print(tabulate([['NR', np.mean(bx_nrtd[n][lam]), np.mean(bx_nrtd[n][lam2]), np.mean(bx_nrtd[n][lam3]), np.mean(bx_nrtd[n][lam4])],['NRBS', np.mean(bx_nrbtd[n][lam]), np.mean(bx_nrbtd[n][lam2]), np.mean(bx_nrbtd[n][lam3]), np.mean(bx_nrbtd[n][lam4])],['BS', np.mean(bx_btd[n][lam]), np.mean(bx_btd[n][lam2]), np.mean(bx_btd[n][lam3]), np.mean(bx_btd[n][lam4])],['MSK', np.mean(bx_msktd[n][lam]), np.mean(bx_msktd[n][lam2]), np.mean(bx_msktd[n][lam3]), np.mean(bx_msktd[n][lam4])]], headers=['lam1=0.1', 'lam2=0.5', 'lam3=1', 'lam4=5']))


from scipy import stats
sample_size = 10
t_critical = stats.t.ppf(q = 0.9, df=sample_size)

yer_nr = []
yer_nrbs = []
yer_bs = []
yer_msk = []

nr_bars = []
nrbs_bars = []
bs_bars = []
msk_bars = []

for lam in lams:
    nr_bars.append(np.mean(bx_nrtd[n][lam]))
    sample_stdev = np.std(bx_nrtd[n][lam])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_nr.append(t_critical * sigma)

    nrbs_bars.append(np.mean(bx_nrbtd[n][lam]))
    sample_stdev = np.std(bx_nrbtd[n][lam])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_nrbs.append(t_critical * sigma)

    bs_bars.append(np.mean(bx_btd[n][lam]))
    sample_stdev = np.std(bx_btd[n][lam])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_bs.append(t_critical * sigma)

    msk_bars.append(np.mean(bx_msktd[n][lam]))
    sample_stdev = np.std(bx_msktd[n][lam])    # Get the sample standard deviation
    sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
    yer_msk.append(t_critical * sigma)


plt.figure()

barWidth = 0.2
rnr = np.arange(len(nr_bars))
rnrbs = [x + barWidth for x in rnr]
rbs = [x + 2*barWidth for x in rnr]
rmsk = [x + 3*barWidth for x in rnr]

plt.gca().set_title('25000 Nodes')
plt.bar(rnr, nr_bars, width = barWidth, hatch='///', edgecolor = 'black', color='w', yerr=yer_nr, capsize=7, label='NR')
plt.bar(rnrbs, nrbs_bars, width = barWidth, hatch='\\\\\\', edgecolor = 'black', color='w',yerr=yer_nrbs, capsize=7, label='NR-BSC')
plt.bar(rbs, bs_bars, width = barWidth, hatch='xxx', edgecolor = 'black',color='w', yerr=yer_bs, capsize=7, label='BSC')
plt.bar(rmsk, msk_bars, width = barWidth, hatch='...', edgecolor = 'black', color='w',yerr=yer_msk, capsize=7, label='MOSEK')
plt.ylabel('Time(s)')
plt.xticks([(r + 2*barWidth) for r in range(len(nr_bars))], ['0.1', '0.5', '1', '5'])
plt.xlabel('Lambda')
plt.legend()
plt.tight_layout()
plt.savefig('plots/' + 'something')











nodes = [25000]
var = 20
lam = 0.7
for n in nodes:
    i = 0
    arcs = [n*2, n*4, n*8]


    nr_time_dict[n] = {}
    bs_time_dict[n] = {}
    nrbs_time_dict[n] = {}
    msk_time_dict[n] = {}

    nr_gap_dict[n] = {}
    bs_gap_dict[n] = {}
    nrbs_gap_dict[n] = {}
    msk_gap_dict[n] = {}

    bx_nrtd[n] = {}
    bx_btd[n] = {}
    bx_nrbtd[n] = {}
    bx_msktd[n] = {}

    for arc in arcs:
        experiment_name = 'network_density'

        if arc == n*4:
            experiment_name = 'performance'
        
        print(n, arc)
        i += 1

        bx_nrtd[n][arc] = []
        bx_btd[n][arc] = []
        bx_nrbtd[n][arc] = []
        bx_msktd[n][arc] = []


        nr_gap_dict[n][arc] = {}
        bs_gap_dict[n][arc] = {}
        nrbs_gap_dict[n][arc] = {}
        msk_gap_dict[n][arc] = {}

        nr_time_dict[n][arc] = {}
        bs_time_dict[n][arc] = {}
        nrbs_time_dict[n][arc] = {}
        msk_time_dict[n][arc] = {}

        setup_times_cvx_l = []
        solve_times_nr_l = []
        solve_times_bs_l = []
        solve_times_alg3_l = []

        for seed in seeds:

            lamstr=str(lam)
            lamstr = lamstr.replace(".","-")
            if lam == 1:
                lamstr_c = '1-0'
            elif lam == 0.7:
                lamstr_c = '0-7'
            else:
                lamstr_c = '0-1'

            save_extension = lamstr + '_' + str(mu) + str(var) + '_' + str(arc) + '_' + str(seed) + str(d_top)

            times_alg3 = np.array(load('times_alg3' + save_extension, n, experiment_name))
            times_nr = np.array(load('times_nr' + save_extension, n, experiment_name))
            times_bs = np.array(load('times_bs' + save_extension, n, experiment_name))

            logfilename = 'saved_runs/' + str(n) + '/' + experiment_name + lamstr_c + '_' + str(mu) + str(var) +  '_' + str(arc) + '_' + str(seed)
            log = read_log(logfilename)
            overhead_time = np.array(load('mosek_overhead_time' + save_extension, n, experiment_name))
            log = log[log['PRSTATUS']>=0.997]
            log = log[log['PRSTATUS']<=1.003]
            cvx_objs = log['POBJ'].tolist()
            times_cvx = log['TIME'].tolist()
            times_cvx = np.array(times_cvx) + overhead_time

            bs_objs = np.array(load('bs_objs' + save_extension, n, experiment_name))
            nr_objs = np.array(load('nr_objs' + save_extension, n, experiment_name))
            alg3_objs = np.array(load('alg3_objs' + save_extension, n, experiment_name))

            bs_min = min(bs_objs)
            nr_min = min(nr_objs)
            alg3_min = min(alg3_objs)
            cvx_min = cvx_objs[-1]

            benchmark = min(bs_min, nr_min, alg3_min, cvx_min)
            descent_to = max(100*((nr_min/benchmark)-1), 100*((bs_min/benchmark)-1), 100*((alg3_min/benchmark)-1), 100*((cvx_min/benchmark)-1))

            print('min percentage diff: %', descent_to)

            bs_i = -1 #find_eq_obj_ind(bs_objs, benchmark, descent_to)
            nr_i = -1 #find_eq_obj_ind(nr_objs, benchmark, descent_to)
            alg3_i = -1 #find_eq_obj_ind(alg3_objs, benchmark, descent_to)

            nr_objs = nr_objs[:nr_i]
            bs_objs = bs_objs[:bs_i]
            alg3_objs = alg3_objs[:alg3_i]

            times_alg3 = times_alg3[:alg3_i]
            times_nr = times_nr[:nr_i]
            times_bs = times_bs[:bs_i]

            gap_bs = abs(100*(np.array(bs_objs)/benchmark -1))
            gap_nr = abs(100*(np.array(nr_objs)/benchmark -1))
            gap_alg3 = abs(100*(np.array(alg3_objs)/benchmark -1))
            gap_cvx = abs(100*(np.array(cvx_objs)/benchmark -1))
            setup_times_cvx_l.append(times_cvx[-1])
            solve_times_nr_l.append(times_nr[-1])
            solve_times_bs_l.append(times_bs[-1])
            solve_times_alg3_l.append(times_alg3[-1])

            nr_time_dict[n][arc][seed] = times_nr
            nrbs_time_dict[n][arc][seed] = times_alg3
            bs_time_dict[n][arc][seed] = times_bs
            msk_time_dict[n][arc][seed] = times_cvx

            nr_gap_dict[n][arc][seed] = gap_nr
            nrbs_gap_dict[n][arc][seed] = gap_alg3
            bs_gap_dict[n][arc][seed] = gap_bs
            msk_gap_dict[n][arc][seed] = gap_cvx

        bx_nrtd[n][arc] = solve_times_nr_l
        bx_btd[n][arc] = solve_times_bs_l
        bx_nrbtd[n][arc] = solve_times_alg3_l
        bx_msktd[n][arc] = setup_times_cvx_l



plt.figure()

# var = 0.7

# two_bars= []
# yer_2brs = []
# four_bars =[]
# yer_4brs =[]
# eight_bars =[]
# yer_8brs =[]

# yer_nr = []
# yer_nrbs = []
# yer_bs = []
# yer_msk = []

# nr_bars = []
# nrbs_bars = []
# bs_bars = []
# msk_bars = []


# for n in [10000]:
#     arcs = [n*2, n*4, n*8]
#     for arc in arcs:
#         # if arc == n*2:
#         nr_bars.append(np.mean(bx_nrtd[n][arc]))
#         sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
#         sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         yer_nr.append(t_critical * sigma)

#         bs_bars.append(np.mean(bx_btd[n][arc]))
#         sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
#         sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         yer_bs.append(t_critical * sigma)

#         nrbs_bars.append(np.mean(bx_nrbtd[n][arc]))
#         sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
#         sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         yer_nrbs.append(t_critical * sigma)

#         msk_bars.append(np.mean(bx_msktd[n][arc]))
#         sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
#         sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         yer_msk.append(t_critical * sigma)
#         # elif arc == n*4:
#         #     four_bars.append(np.mean(bx_nrtd[n][arc]))
#         #     sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_4brs.append(t_critical * sigma)

#         #     four_bars.append(np.mean(bx_btd[n][arc]))
#         #     sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_4brs.append(t_critical * sigma)

#         #     four_bars.append(np.mean(bx_nrbtd[n][arc]))
#         #     sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_4brs.append(t_critical * sigma)



#         #     four_bars.append(np.mean(bx_msktd[n][arc]))
#         #     sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_4brs.append(t_critical * sigma)
#         # else:
#         #     eight_bars.append(np.mean(bx_nrtd[n][arc]))
#         #     sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_8brs.append(t_critical * sigma)

#         #     eight_bars.append(np.mean(bx_btd[n][arc]))
#         #     sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_8brs.append(t_critical * sigma)

#         #     eight_bars.append(np.mean(bx_nrbtd[n][arc]))
#         #     sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_8brs.append(t_critical * sigma)

#         #     eight_bars.append(np.mean(bx_msktd[n][arc]))
#         #     sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
#         #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
#         #     yer_8brs.append(t_critical * sigma)


# r2 = np.arange(len(two_bars))
# # r4 = [x + barWidth for x in r2]
# r8 = [x + 2*barWidth for x in r2]


barWidth = 0.2
rnr = np.arange(len(nr_bars))
rnrbs = [x + barWidth for x in rnr]
rbs = [x + 2*barWidth for x in rnr]
rmsk = [x + 3*barWidth for x in rnr]
#
 # color = 'red',
# plt.bar(rnr, nr_bars, width = barWidth,  edgecolor = 'black', yerr=yer_nr, capsize=7, label='NR')
# plt.bar(rbs, bs_bars, width = barWidth,  edgecolor = 'black', yerr=yer_bs, capsize=7, label='BSC')
# plt.bar(rnrbs, nrbs_bars, width = barWidth,  edgecolor = 'black', yerr=yer_nrbs, capsize=7, label='NR-BSC')
# plt.bar(rmsk, msk_bars, width = barWidth, edgecolor = 'black', yerr=yer_msk, capsize=7, label='MOSEK')
# plt.ylabel('Time(s)')
# plt.xticks([(r + 3*barWidth) for r in range(len(nr_bars))], ['5000', '10000', '25000'])
# plt.xlabel('# of Nodes')
# plt.legend()
# plt.tight_layout()
# plt.savefig('plots/' + 'comparison')
# ax = plt.subplot(2, 1, 1)
# plt.gca().set_title('10000 Nodes')


# plt.bar(rnr, nr_bars, width = barWidth, hatch='///', edgecolor = 'black', color='w', yerr=yer_nr, capsize=7, label='NR')
# plt.bar(rnrbs, nrbs_bars, width = barWidth, hatch='\\\\\\', edgecolor = 'black', color='w',yerr=yer_nrbs, capsize=7, label='NR-BSC')
# plt.bar(rbs, bs_bars, width = barWidth, hatch='xxx', edgecolor = 'black',color='w', yerr=yer_bs, capsize=7, label='BSC')
# plt.bar(rmsk, msk_bars, width = barWidth, hatch='...', edgecolor = 'black', color='w',yerr=yer_msk, capsize=7, label='MOSEK')

# plt.bar(rnr, two_bars, width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer_2brs[0], capsize=3, label='NR')
# plt.bar(rbs, two_bars[1], width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer_2brs[1], capsize=3, label='BSC')
# plt.bar(rnrbs, two_bars[2] , width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer_2brs[2], capsize=3, label='NR-BSC')
# plt.bar(rmsk, two_bars[3], width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer_2brs[3], capsize=3, label='MSK')

# plt.ylabel('Time(s)')
# plt.xticks([(r + 2*barWidth) for r in range(len(nr_bars))], ['2', '4', '8'])
# plt.xlabel('Average Degree of Nodes')
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots/' + str(n) + 'densityfig10kn')



# plt.bar(r2, two_bars, width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer_2brs, capsize=3)
# plt.bar(r4, four_bars, width = barWidth, color = 'orange', edgecolor = 'black', yerr=yer_4brs, capsize=3)
# plt.bar(r8, eight_bars, width = barWidth, color = 'purple', edgecolor = 'black', yerr=yer_8brs, capsize=3)
# plt.ylabel('Time(s)')
# plt.xticks([r for r in range(len(r4))], ['NR', 'BSC', 'NR-BSC', 'MSK'])
# plt.xlabel('Algorithm')
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots/' + str(n) + 'densityfig')


# print(tabulate([['NR', (100*eight_bars[0]/two_bars[0]) - 100], ['BS', (100*eight_bars[1]/two_bars[1]) - 100], ['NRBS', (100*eight_bars[2]/two_bars[2]) - 100],['MSK', (100*eight_bars[3]/two_bars[3]) - 100]], headers=['%inc']))
# print(tabulate([['NR', (100*four_bars[0]/two_bars[0]) - 100], ['BS', (100*four_bars[1]/two_bars[1]) - 100], ['NRBS', (100*four_bars[2]/two_bars[2]) - 100],['MSK', (100*four_bars[3]/two_bars[3]) - 100]], headers=['%inc']))


yer_nr = []
yer_nrbs = []
yer_bs = []
yer_msk = []

nr_bars = []
nrbs_bars = []
bs_bars = []
msk_bars = []


for n in [25000]:
    arcs = [n*2, n*4, n*8]
    for arc in arcs:
        # if arc == n*2:
        nr_bars.append(np.mean(bx_nrtd[n][arc]))
        sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
        sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        yer_nr.append(t_critical * sigma)

        bs_bars.append(np.mean(bx_btd[n][arc]))
        sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
        sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        yer_bs.append(t_critical * sigma)

        nrbs_bars.append(np.mean(bx_nrbtd[n][arc]))
        sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
        sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        yer_nrbs.append(t_critical * sigma)

        msk_bars.append(np.mean(bx_msktd[n][arc]))
        sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
        sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        yer_msk.append(t_critical * sigma)
            # two_bars.append(np.mean(bx_nrtd[n][arc]))
            # sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
            # sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
            # yer_2brs.append(t_critical * sigma)

            # two_bars.append(np.mean(bx_btd[n][arc]))
            # sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
            # sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
            # yer_2brs.append(t_critical * sigma)

            # two_bars.append(np.mean(bx_nrbtd[n][arc]))
            # sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
            # sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
            # yer_2brs.append(t_critical * sigma)

            # two_bars.append(np.mean(bx_msktd[n][arc]))
            # sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
            # sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
            # yer_2brs.append(t_critical * sigma)
        # elif arc == n*4:
        #     four_bars.append(np.mean(bx_nrtd[n][arc]))
        #     sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_4brs.append(t_critical * sigma)

        #     four_bars.append(np.mean(bx_btd[n][arc]))
        #     sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_4brs.append(t_critical * sigma)

        #     four_bars.append(np.mean(bx_nrbtd[n][arc]))
        #     sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_4brs.append(t_critical * sigma)



        #     four_bars.append(np.mean(bx_msktd[n][arc]))
        #     sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_4brs.append(t_critical * sigma)
        # else:
        #     eight_bars.append(np.mean(bx_nrtd[n][arc]))
        #     sample_stdev = np.std(bx_nrtd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_8brs.append(t_critical * sigma)

        #     eight_bars.append(np.mean(bx_btd[n][arc]))
        #     sample_stdev = np.std(bx_btd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_8brs.append(t_critical * sigma)

        #     eight_bars.append(np.mean(bx_nrbtd[n][arc]))
        #     sample_stdev = np.std(bx_nrbtd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_8brs.append(t_critical * sigma)

        #     eight_bars.append(np.mean(bx_msktd[n][arc]))
        #     sample_stdev = np.std(bx_msktd[n][arc])    # Get the sample standard deviation
        #     sigma = sample_stdev/np.sqrt(sample_size)  # Standard deviation estimate
        #     yer_8brs.append(t_critical * sigma)

# r2 = np.arange(len(two_bars))
# r4 = [x + barWidth for x in r2]
# r8 = [x + 2*barWidth for x in r2]

barWidth = 0.2
rnr = np.arange(len(nr_bars))
rnrbs = [x + barWidth for x in rnr]
rbs = [x + 2*barWidth for x in rnr]
rmsk = [x + 3*barWidth for x in rnr]

# plt.subplot(2, 2, 3)
# plt.gca().set_title('25000 Nodes')

# ax = plt.subplot(2, 1, 2)
pdb.set_trace()
plt.gca().set_title('25000 Nodes')

plt.bar(rnr, nr_bars, width = barWidth, hatch='///', edgecolor = 'black', color='w', yerr=yer_nr, capsize=7, label='NR')
plt.bar(rnrbs, nrbs_bars, width = barWidth, hatch='\\\\\\', edgecolor = 'black', color='w',yerr=yer_nrbs, capsize=7, label='NR-BSC')
plt.bar(rbs, bs_bars, width = barWidth, hatch='xxx', edgecolor = 'black',color='w', yerr=yer_bs, capsize=7, label='BSC')
plt.bar(rmsk, msk_bars, width = barWidth, hatch='...', edgecolor = 'black', color='w',yerr=yer_msk, capsize=7, label='MOSEK')

plt.ylabel('Time(s)')
plt.xticks([(r + 2*barWidth) for r in range(len(nr_bars))], ['2', '4', '8'])
plt.xlabel('Average Degree of Nodes')
plt.tight_layout()
plt.legend()
plt.savefig('plots/' + str(n) + 'densityfig25kn')

# plt.bar(r2, two_bars, width = barWidth, color = 'cyan', edgecolor = 'black', yerr=yer_2brs, capsize=3)
# plt.bar(r4, four_bars, width = barWidth, color = 'orange', edgecolor = 'black', yerr=yer_4brs, capsize=3)
# plt.bar(r8, eight_bars, width = barWidth, color = 'purple', edgecolor = 'black', yerr=yer_8brs, capsize=3)
# plt.ylabel('Time(s)')
# plt.xticks([r for r in range(len(r4))], ['NR', 'BSC', 'NR-BSC', 'MSK'])
# plt.xlabel('Algorithm')
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots/' + str(n) + 'densityfig')

# plt.show()


# print(tabulate([['NR', (100*eight_bars[0]/two_bars[0]) - 100], ['BS', (100*eight_bars[1]/two_bars[1]) - 100], ['NRBS', (100*eight_bars[2]/two_bars[2]) - 100],['MSK', (100*eight_bars[3]/two_bars[3]) - 100]], headers=['%inc']))
# print(tabulate([['NR', (100*four_bars[0]/two_bars[0]) - 100], ['BS', (100*four_bars[1]/two_bars[1]) - 100], ['NRBS', (100*four_bars[2]/two_bars[2]) - 100],['MSK', (100*four_bars[3]/two_bars[3]) - 100]], headers=['%inc']))



# pdb.set_trace()





    # try:
    #     cvx.append(np.mean(setup_times_cvx_l))
    #     bs.append(np.mean(solve_times_bs_l))
    #     nr.append(np.mean(solve_times_nr_l))
    #     alg3.append(np.mean(solve_times_alg3_l))
    # except:
    #     pdb.set_trace()


# fig = ax.get_figure()

df = pd.DataFrame({'MOSEK': cvx, 'BS': bs, 'NR': nr, 'NR-BS': alg3}, index=lams)

ax = df.plot.bar(rot=0)
ax.set_ylabel('Time')
ax.set_xlabel('Lambda')
ax.set_title('Lambda Effect')
fig = ax.get_figure()
fig.show()
fig.savefig("varying_lambda node_no: " + str(n))
pdb.set_trace()

    # plt.figure()
    # #pdb.set_trace()
    # plt.plot(gap_bs_l[0], label='BSC', marker='o')
    # plt.plot(gap_nr_l[0], label='NR', marker='o')
    # plt.plot(gap_alg3_l[0], label='NR-BSC', marker='o')
    #
    # # plt.yscale('log')
    # plt.title('Number of nodes: ' + str(n))
    # plt.legend(loc='upper right')
    # plt.xlabel('Iteration number')
    # plt.ylabel('Objective Gap')
    # plt.savefig('plots/' + str(n))

### Experiment Analysis ###
### Experiment 3 ###
experiment_name = 'network_density'
lam = 0.1
mu = 10
variances = 10
num_arcs = [n*2, n*4, n*8]
seeds=8 * np.arange(10)

for arcs in num_arcs:
    for seed in seeds:
        pass

### Experiment Analysis ###
### Experiment 4 ###
experiment_name = 'performance'
lam = 0.1
mu = 10
var = 10
arcs = n*2
seeds=99 * np.arange(100)
data = []

for seed in seeds:

    lamstr=str(lam)
    lamstr = lamstr.replace(".","-")

    save_extension = lamstr + '_' + str(mu) + str(var) + '_' + str(arcs) + '_' + str(seed)

    times_alg3 = np.array(load('times_alg3' + save_extension, n, experiment_name))
    times_nr = np.array(load('times_nr' + save_extension, n, experiment_name))
    times_bs = np.array(load('times_bs' + save_extension, n, experiment_name))
    solve_time_cvx = np.array(load('solve_time_cvx' + save_extension, n, experiment_name))

    best_alg = str(np.array(load('best_soln_algo' + save_extension, n, experiment_name)))
    best_obj = float(np.array(load('benchmark_soln' + save_extension, n, experiment_name)))

    logfilename = 'saved_runs/' + str(n) + '/' + lamstr + '_' + str(mu) + str(var) +  '_' + str(arcs) + '_' + str(seed)
    log = read_log(logfilename)

    bs_objs = np.array(load('bs_objs' + save_extension, n, experiment_name))
    nr_objs = np.array(load('nr_objs' + save_extension, n, experiment_name))
    alg3_objs = np.array(load('alg3_objs' + save_extension, n, experiment_name))

    bs_i = find_eq_obj_ind(bs_objs, best_obj)
    nr_i = find_eq_obj_ind(nr_objs, best_obj)
    alg3_i = find_eq_obj_ind(alg3_objs, best_obj)

    print(times_bs[bs_i], times_alg3[alg3_i], times_nr[nr_i], float(solve_time_cvx))
    data.append([times_bs[bs_i], times_alg3[alg3_i], times_nr[nr_i], float(solve_time_cvx)])

df = pd.DataFrame(data,columns=['BSC','NR-BSC','NR','CVX'])
ax = df.plot(kind='box', figsize=[16, 8])
fig = ax.get_figure()
fig.show()
pdb.set_trace()


cvx_dict= {}
cvx = []
bs = []
nr = []
alg3 = []

for n in [5000]:

    lams = [0.1, 0.5, 5, 10, 100]
    mus = [1, 10, 100]
    variances = [1, 10, 50, 100]
    num_arcs = [n*2, n*4, n*8]
    for mu in mus:
        for var in variances:
            for seed in seeds:
                for lam in lams:
                    print('mu - var - lam - seed', mu, var, lam, seed, '============')
                    # lam = lams[0]
                    # mu = mus[0]
                    # var = variances[0]
                    arcno = num_arcs[0]
                    # seed = seeds[0]
                    # i = iss[0]
                    save_extension = str(lam) + '_' + str(mu) + str(var) + '_' + str(arcno) + '_' + str(seed)

                    times_alg3 = np.array(load('times_alg3' + save_extension, n))
                    times_nr = np.array(load('times_nr' + save_extension, n))
                    times_bs = np.array(load('times_bs' + save_extension, n))
                    solve_time_cvx = np.array(load('solve_time_cvx' + save_extension, n))

                    print('=======')
                    print(arcno)
                    print(times_bs[-1], times_alg3[-1], times_nr[-1], solve_time_cvx)



    times_cvx_l = []
    gap_cvx_l= []
    gap_bs_l= []
    times_bs_l= []
    bs_iters_l= []
    gap_nr_l= []
    times_nr_l= []
    nr_iters_l= []
    gap_alg3_l= []
    times_alg3_l= []
    alg3_iters_l= []
    solve_times_cvx_l = []
    setup_times_cvx_l = []
    solve_times_nr_l = []
    solve_times_bs_l = []
    solve_times_alg3_l = []
    
    
    myrange =10
    if n==75000:
        myrange = 5
    for i in range(myrange):
        solve_times_cvx_l.append(load('solve_time_cvx' + str(i), n))
        # setup_times_cvx_l.append(load('setup_time_cvx' + str(i), n))
    
        gap_bs = np.array(load('gap_bs' + str(i), n))
        gap_bs[gap_bs<0.01]=0
        gap_bs_l.append(gap_bs)
        times_bs_l.append(load('times_bs' + str(i), n))
        solve_times_bs_l.append(load('times_bs' + str(i), n)[-1])
        # bs_subproblem_times=load(bs_subproblem_times, 'bs_subproblem_times' + str(i), n)
        bs_iters = np.array(load('bs_iters' + str(i), n)) + 1
        bs_iters_l.append(bs_iters)
    
        gap_nr = np.array(load('gap_nr' + str(i), n))
        gap_nr[gap_nr<0.01]=0
        gap_nr_l.append(gap_nr)
        times_nr_l.append(load('times_nr' + str(i), n))
        solve_times_nr_l.append(load('times_nr' + str(i), n)[-1])
        nr_iters_l.append(load('nr_iters' + str(i), n))
        # nr_subproblem_times=load(nr_subproblem_times, 'nr_subproblem_times' + str(i), n)
    
        gap_alg3 = np.array(load('gap_alg3' + str(i), n))
        gap_alg3[gap_alg3<0.01]=0
        gap_alg3_l.append(gap_alg3)
        times_alg3_l.append(load('times_alg3' + str(i), n))
        solve_times_alg3_l.append(load('times_alg3' + str(i), n)[-1])
        alg3_iters_l.append(load('alg3_iters' + str(i), n))
        # alg3_subproblem_times=load(alg3_subproblem_times, 'alg3_subproblem_times' + str(i), n)
    
    
    cvx.append(np.mean(solve_times_cvx_l))
    bs.append(np.mean(solve_times_bs_l))
    nr.append(np.mean(solve_times_nr_l))
    alg3.append(np.mean(solve_times_alg3_l))
    
    print('cvx ' + str(n) + ' ' + str(np.mean(solve_times_cvx_l)) + '' )
    print('bs ' + str(n) + ' ' + str(np.mean(solve_times_bs_l))+ ' ' + str(np.mean(bs_iters_l)))
    print('nr ' + str(n) + ' ' + str(np.mean(solve_times_nr_l)) + ' ' + str(np.mean(nr_iters_l)))
    print('alg3 ' + str(n) +  ' ' + str(np.mean(solve_times_alg3_l))+ ' ' +str(np.mean(alg3_iters_l)))
    
    print('======================/=============================')
    

    pdb.set_trace()

    bs_dict[n] = {}
    bs_dict[n]['mean solve time'] = np.mean(solve_times_bs_l)

    nr_dict[n] = {}
    nr_dict[n]['mean solve time'] = np.mean(solve_times_nr_l)

    alg3[n] = {}
    alg3[n]['mean solve time'] = np.mean(solve_times_alg3_l)


    cvx_dict[n]['mean setup time'] = np.mean(setup_times_cvx_l)



    print(n)
    # gap_cvx, times_cvx = get_level_time(gap_cvx_l, times_cvx_l)
    gap_bs, times_bs = get_level_time(gap_bs_l, times_bs_l)
    gap_nr, times_nr = get_level_time(gap_nr_l, times_nr_l)
    gap_alg3, times_alg3 = get_level_time(gap_alg3_l, times_alg3_l)


    df = pd.DataFrame({'BS': times_bs[::-1], 'NR': times_nr[::-1], 'NR-BS': times_alg3[::-1]}, index=gap_nr[::-1])
    ax = df.plot.bar(rot=0)

    ax.set_ylabel('Time')
    ax.set_xlabel('Gap')
    ax.set_title('Time ')
    fig = ax.get_figure()
    fig.savefig("network_" + str(n) )
pdb.set_trace()


plt.figure()
#pdb.set_trace()
plt.plot(gap_bs_l[0], label='BSC', marker='o')
plt.plot(gap_nr_l[0], label='NR', marker='o')
plt.plot(gap_alg3_l[0], label='NR-BSC', marker='o')

# plt.yscale('log')
plt.title('Number of nodes: ' + str(n))
plt.legend(loc='upper right')
plt.xlabel('Iteration number')
plt.ylabel('Objective Gap')
plt.savefig('plots/' + str(n))

plt.figure()
#pdb.set_trace()
plt.plot(gap_bs_l[2], label='BSC', marker='o')
plt.plot(gap_nr_l[2], label='NR', marker='o')
plt.plot(gap_alg3_l[2], label='NR-BSC', marker='o')

# plt.yscale('log')
plt.title('Number of nodes: ' + str(n))
plt.legend(loc='upper right')
plt.xlabel('Iteration number')
plt.ylabel('Objective Gap')
plt.savefig('plots/' + str(n) + '_different_run')

plt.figure()
#pdb.set_trace()
plt.plot(times_bs_l[2],gap_bs_l[2], label='BSC', marker='o', linestyle='dashed')
plt.plot(times_nr_l[2],gap_nr_l[2], label='NR', marker='D', linestyle='dashed')
#plt.plot(times_alg3_l[2], gap_alg3_l[2],label='NR-BSC', marker='o')

# plt.yscale('log')
plt.title('Number of nodes: ' + str(n))
plt.legend(loc='upper right')
plt.xlabel('Time (s)')
plt.ylabel('Objective Gap')
plt.savefig('plots/' + str(n) + 'timevsiter')
pdb.set_trace()
##do stacked bar plot for cvx over number of nodes
#pdb.set_trace()
#pdb.set_trace()

nodesnumbers = [5000,10000,25000,50000,75000]


df = pd.DataFrame({'Mosek': cvx,'BSC': bs, 'NR': nr, 'NR-BSC': alg3 }, index=nodesnumbers)
ax = df.plot.bar(rot=0)

ax.set_ylabel('Time (seconds)')
ax.set_xlabel('# Nodes')
ax.set_title('')
#for i in ax.patches:
#    # get_width pulls left or right; get_y pushes up or down
#    ax.text(i.get_x()+.04,i.get_height()+12000,"BSC","MOSEK","NR","NR-BSC", fontsize=11, color='dimgrey')
fig = ax.get_figure()
fig.savefig("comparison")





    # TODO:
    # maybe no presolve or cut in early iters for cvx

    # pdb.set_trace()
    # pdb.set_trace()


lams = np.linspace(0, 5, num=50)
n = 5000
arcs = 40000
mu=20
var=60
d_top = 1
mucosts=[]
varcosts=[]
seed = 0
lam_costs = {}
iters = 0 
vallams = []
for lam in lams:
    if (iters%1==0 and iters<=10) or (iters%5==0 and iters>10) or lam==lams[-1] :
        lamstr = str(round(lam, 6))
        lamstr = lamstr.replace('.','')

        save_extension = str(mu) + str(var) + '_' + str(arcs) + '_' + str(seed) + str(d_top)

        experiment_name = 'varying_lams'
        mucost = load( 'mu_cost_' + lamstr + '_' + save_extension , n, experiment_name)
        varcost = load( 'var_cost_' + lamstr + '_' + save_extension, n, experiment_name)
        mucosts.append(mucost)
        varcosts.append(varcost)
        vallams.append(lam)
        lam_costs[lam] = (mucost, varcost)
    iters += 1

x=[]
y=[]
iters = 0

for lam in lams:
    if (iters%1==0 and iters<=10) or (iters%5==0 and iters>10) or lam==lams[-1] :
        y.append(lam_costs[lam])
        x.append(lam)
    iters +=1

fig, ax1 = plt.subplots()

# for xe, ye in zip(x, y):
#     plt.scatter([xe] * len(ye), ye)
ax2 = ax1.twinx()
pdb.set_trace()
ax1.scatter(vallams, mucosts)
ax2.scatter(vallams, varcosts)

# plt.xticks(list(np.arange(lams)))
def major_formatter(x, pos):
    return "[%.2f]" % x
import matplotlib
ax1.set_xticks([lams[0], lams[5], lams[10], lams[15], lams[20], lams[25], lams[30], lams[35], lams[40], lams[45], lams[-1]])
ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%0.1f'))
# ax1.xaxis.set_ticks

lns1 = ax1.plot(vallams, mucosts, label='Mean')
lns2 = ax2.plot(vallams, varcosts, '--', label='Standard Deviation')


# added these three lines
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=1)


ax1.set_xlabel('Value of Reliability')
ax1.set_ylabel('Mean Cost')
ax2.set_ylabel('Standard Deviation Cost')


# plt.plot(lams, mucosts, label='Mean')
# plt.plot(lams, varcosts, '-', label='Standard Deviation')
# plt.legend(loc='upper right')
# plt.xlabel('Value of Reliability')
# plt.ylabel('Mean Cost')
# plt.ylabel('Mean Cost')

# ax2.set_ylabel('Y2 data', color='b')
plt.savefig('plots/' + 'meanvsstd')
# plt.show()


