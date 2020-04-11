import matplotlib.pyplot as plt
import pandas
import copy
import pdb
import pandas as pd
from utils import *
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mtick


# plt.plot([7813, 7813], [0., 0.9], "r:")         # Not shown
# plt.plot([-50000, 7813], [0.9, 0.9], "r:")      # Not shown
# plt.plot([-50000, 7813], [0.4368, 0.4368], "r:")# Not shown
# plt.plot([7813], [0.9], "ro")                   # Not shown
# plt.plot([7813], [0.4368], "ro")                # Not shown

def am_i_close(a, b, perc_rel_tol):
    if 100*abs(a - b) / min(a, b) <= perc_rel_tol:
        return True
    else:
        return False

def save_fig(fig_id, tight_layout=False, fig_extension="png", resolution=300):
    path = os.path.join(PLOTS_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)

    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)

def mean_std(array):

    # mean = array.nanmean(axis=1)
    mean = np.nanmean(np.where(array != 0, array, np.nan), 1)
    std = np.nanstd(np.where(array != 0, array, np.nan), 1)
    # std = array.nanstd(axis=1)
    print(mean)
    return mean, std

def graph_family(family_name, cplex_times, nr_times, bsc_times, nr_iters_, bs_iters_, cplex_infeas_):
    node_num = [4096, 8192, 8192*2, 8192*4]
    node_num = [r'$2^{12}$', r'$2^{13}$', r'$2^{14}$', r'$2^{15}$']
    print(family_name)

    print('cplex')
    cplex_mean, cplex_std = mean_std(cplex_times)
    print('nr')
    nr_mean, nr_std = mean_std(nr_times)
    print('bsc')
    bsc_mean, bsc_std = mean_std(bsc_times)

    print(cplex_infeas_, nr_iters_, bs_iters_)
    plt.plot(node_num, cplex_mean, 'b--o',
             linewidth=1, label='CPLEX', markersize=4)
    plt.plot(node_num, bsc_mean, 'r--^',
             linewidth=1, label='BSC', markersize=4)
    plt.plot(node_num, nr_mean, 'g--s', linewidth=1, label='NR', markersize=4)

    plt.text(node_num[0], cplex_mean[0], 'CPLEX')
    plt.text(node_num[0], nr_mean[0], 'NR')
    plt.text(node_num[0], bsc_mean[0], 'BSC')

    # plt.fill_between(node_num, cplex_mean - cplex_std, cplex_mean + cplex_std, facecolor='blue', alpha=0.1);
    # plt.fill_between(node_num, nr_mean - nr_std, nr_mean + nr_std, facecolor='green', alpha=0.1);
    # plt.fill_between(node_num, bsc_mean - bsc_std, bsc_mean + bsc_std, facecolor='red', alpha=0.1);

    plt.ylabel('Running time (logscale)', fontsize=12)
    plt.xlabel('Number of nodes', fontsize=12)
    plt.yscale('log')
    # plt.xticks([])
    # plt.tick_params(bottom='off')
    # plt.xscale('log')
    # plt.xticks(np.arange(min(np.log(node_num)), max(np.log(node_num)), 1.0))
    # labels_ = [r'$2^{12}$', r'$2^{13}$', r'$2^{14}$']
    # plt.xticks(node_num, labels_, rotation=45)  # Set text labels.

    if family_name.lower().find('sr') >= 0:
        plt.ylim(top=1e4)
        plt.ylim(bottom=1e0)

        # plt.xlim(xmin=4000)

    else:
        plt.ylim(top=1e4)
        plt.ylim(bottom=1e0)
        
        # plt.xlim(xmin=4000)
    # plt.axis([xs,xe,ys,ye])

    plt.grid(True)
    plt.title(family_name, fontsize=12)
    plt.legend(loc='upper left', fontsize=10)

def graph_base_vs_reliable():
    # graph base vs reliable
    experiment = 'base_vs_reliable'
    # lams = [0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7.5, 10, 25, 50, 75, 100]#, 250, 500, 750, 1000]
    lams = np.logspace(-1, 3, 200)
    lams = [lams[i] for i in range(len(lams)) if i%4==0]
    exponent = '10'
    tail = 'e'
    base = 'netgen'
    atype = '_8_'

    rel_mean = []
    rel_std = []
    rel_obj = []
    base_obj = []

    for lam in lams:

        lam_dir = str(lam).replace('.', '_')

        filename = experiment + '/' + base + '/' + \
                    base + atype + exponent + tail + '_' + lam_dir
        data_dic = load_run(filename)

        rel_mean.append(data_dic['rel_mean'])
        rel_std.append(data_dic['rel_std'])
        rel_obj.append(data_dic['solver_obj'])
        base_obj.append(data_dic['base_obj'])

    fig, ax1 = plt.subplots()

    rel_mean = np.array(rel_mean)/1e6
    rel_std = np.array(rel_std)/1e6

    ax1.set_xlabel(r'$\bar{\lambda}$')
    ax1.set_ylabel('Mean (' + r'$10^6$' + ')')
    ax1.plot(lams, rel_mean, 'b--+', linewidth=1, label='mean', markersize=4)
    ax1.tick_params(axis='y')


    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Std. Dev. (' + r'$10^6$' + ')')  # we already handled the x-label with ax1
    ax2.plot(lams, rel_std, 'r-+', linewidth=1, label='std', markersize=4)

    ax2.tick_params(axis='y')


    plt.grid(True)
    ax1.legend(loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2,
            borderaxespad=0, frameon=False)
    ax2.legend(loc='lower right', bbox_to_anchor= (0.4, 1.01), ncol=2,
            borderaxespad=0, frameon=False)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    # ax1.set_xscale('log')
    # ax2.set_xscale('log')

    save_fig('mean_std_tradeoff')
    plt.clf()


    rel_gap = []
    for i in range(len(rel_obj)):
        rel_gap.append(100*abs(base_obj[i] - rel_obj[i])/min(base_obj[i], rel_obj[i]))

    plt.plot(lams, rel_gap , 'k--+', linewidth=1, label='gap', markersize=4)
    plt.ylabel('Relative Objective Gap', fontsize=10)
    plt.xlabel(r'$\bar{\lambda}$', fontsize=10)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
    # plt.yscale('log')
    # plt.xscale('log')
    plt.grid(True)

    save_fig('obj_tradeoff')
    plt.clf()

def graph_gap_levels():
    

    experiment = 'graph_families'

    bases = ['netgen'] 
    tails = ['e']
    exponents = (np.arange(14, 15)).astype(str)
    types = ['_lo_sr_']


    for base in bases:

        i = 0
        for atype in types:

            if base == 'goto' and atype.find('lo') > 0:
                continue
            i += 1
            # if base == 'netgen':
            #     perc_rel_tol = 1e-3

            #     figure_grid = 220
            #     plt.subplot(figure_grid + i)
            # else:
            #     perc_rel_tol = 1e-3

            #     figure_grid = 120
            #     plt.subplot(figure_grid + i)


            for exponent in exponents:
                bs_iter_list = []
                nr_iter_list = []
                cplex_infeas_list = []
                j = 0

                for tail in tails:
                    row = np.where(exponents == exponent)[0][0]
                    col = j

                    cut_idx = None
                    weird = False
                    network = base + atype + exponent + tail
                    filename = experiment + '/' + base + '/' + network
                    data_dic = load_run(filename)
                    cplex_obj = data_dic['solver_obj']

                    cplex_iter_objs = data_dic['solver_primals']
                    cplex_iter_times = data_dic['solver_iter_times']
                    cplex_feasible = data_dic['solver_feasible']

                    cplex_iter_objs[-1] = cplex_obj

                    nr_iter_objs = data_dic['nr_iter_objs']
                    nr_iter_elapsed = data_dic['nr_iter_elapsed']

                    bs_iter_objs = data_dic['bs_iter_objs']
                    bs_iter_elapsed = data_dic['bs_iter_elapsed']

                    lb_elapsed = data_dic['elapsed_lower_bound']
                    ub_elapsed = data_dic['elapsed_upper_bound']

                    nr_iter_elapsed = lb_elapsed + np.array(nr_iter_elapsed)
                    bs_iter_elapsed = lb_elapsed + \
                        ub_elapsed + np.array(bs_iter_elapsed)

                    for _ind in range(len(cplex_iter_objs) - 1, 0, -1):
                        if not am_i_close(cplex_iter_objs[_ind], cplex_iter_objs[_ind - 1], perc_rel_tol=1e-7):
                            cplex_iter_objs = cplex_iter_objs[:_ind + 1]
                            cplex_iter_times = cplex_iter_times[:_ind + 1]
                            break

                    plt.figure(figsize=(16,8))
                    
                    best_known = min(min(nr_iter_objs), min(bs_iter_objs), cplex_iter_objs[-1])

                    def to_arr(sequence):
                        return np.array(sequence)

                    nr_iter_objs = to_arr(nr_iter_objs)
                    bs_iter_objs = to_arr(bs_iter_objs)
                    cplex_iter_objs = to_arr(cplex_iter_objs)

                    def prep(seq):
                        seq = abs(seq - best_known)
                        seq *= 100
                        seq /= best_known
                        return seq

                    nr_iter_objs = prep(nr_iter_objs)
                    bs_iter_objs = prep(bs_iter_objs)
                    cplex_iter_objs = prep(cplex_iter_objs)

                    plt.plot(nr_iter_elapsed, nr_iter_objs,'g--s', linewidth=1, label='NR')
                    plt.plot(bs_iter_elapsed, bs_iter_objs,'r--^', linewidth=1, label='BSC')
                    plt.plot(cplex_iter_times, cplex_iter_objs,'b--o', linewidth=1, label='CPLEX')

                    # plt.plot([0, nr_times_e[tail][0]], [nr_objs_e[tail][
                    #          0], nr_objs_e[tail][0]], "k:")  # Not shown
                    # plt.plot([nr_times_e[tail][0], nr_times_e[tail][0]], [
                    #          1e-10, nr_objs_e[tail][0]], "k:")  # Not shown
                    # plt.plot([nr_times_e[tail][0], nr_times_e[tail][0]], [
                    #          1e-10, cplex_objs_e[tail][6]], "k:")  # Not shown
                    # plt.plot([0, nr_times_e[tail][0]], [cplex_objs_e[tail][
                    #          6], cplex_objs_e[tail][6]], "k:")  # Not shown

                    # plt.plot([7813], [0.9], "ko")                   # Not shown

                    plt.ylabel('Relative Objective Gap % (logscale)', fontsize=12)
                    plt.yscale('log')
                    plt.ylim(bottom=1e-5)
                    # plt.xlim(left=0)

                    plt.xlabel('Running time', fontsize=12)
                    plt.grid(True)
                    plt.legend(loc='upper right', fontsize=10)

                    save_fig('gap_level')
                    plt.clf()

def graph_comparison():

    experiment = 'graph_families'

    bases = ['netgen']  # , 'goto']
    tails = ['a','b','c','d','e']
    exponents = (np.arange(12, 16)).astype(str)
    types = ['_lo_8_', '_8_', '_lo_sr_', '_sr_']
    temp_tails = ['a','b','c','d','e']

    for base in bases:
        plt.figure(figsize=(16, 8))

        i = 0
        for atype in types:

            if base == 'goto' and atype.find('lo') > 0:
                continue
            i += 1
            if base == 'netgen':
                perc_rel_tol = 1e-3

                figure_grid = 220
                plt.subplot(figure_grid + i)
            else:
                perc_rel_tol = 1e-3

                figure_grid = 120
                plt.subplot(figure_grid + i)

            cplex_times = np.zeros((len(exponents), len(temp_tails)))
            nr_times = np.zeros((len(exponents), len(temp_tails)))
            bsc_times = np.zeros((len(exponents), len(temp_tails)))

            nr_iters_ = np.zeros(len(exponents))
            bs_iters_ = np.zeros(len(exponents))
            cplex_infeas_ = np.zeros(len(exponents))

            for exponent in exponents:

                bs_iter_list = []
                nr_iter_list = []
                cplex_infeas_list = []
                j = 0

                for tail in tails:
                    row = np.where(exponents == exponent)[0][0]
                    col = j

                    cut_idx = None
                    weird = False

                    if exponent != '15':
                        network = base + atype + exponent + tail
                        filename = experiment + '/' + base + '/' + network
                        data_dic = load_run(filename)
                        cplex_times[row, col] = data_dic['solver_elapsed']
                        cplex_obj = data_dic['solver_obj']

                        cplex_iter_objs = data_dic['solver_primals']
                        cplex_duals = data_dic['solver_duals']

                        cplex_gaps = 100*abs(np.array(cplex_iter_objs[:-1]) - np.array(cplex_duals))/np.minimum(np.array(cplex_iter_objs[:-1]),np.array(cplex_duals))

                        cplex_iter_times = data_dic['solver_iter_times']
                        cplex_feasible = data_dic['solver_feasible']

                        cplex_iter_objs[-1] = cplex_obj

                        nr_iter_objs = data_dic['nr_iter_objs']
                        nr_iter_elapsed = data_dic['nr_iter_elapsed']

                        bs_iter_objs = data_dic['bs_iter_objs']
                        bs_iter_elapsed = data_dic['bs_iter_elapsed']

                        lb_elapsed = data_dic['elapsed_lower_bound']
                        ub_elapsed = data_dic['elapsed_upper_bound']
                    else:
                        network = base + atype + exponent + tail
                        filename = experiment + '/' + base + '/' + network + '_cp_results'
                        cp_dic = load_run(filename)

                        filename = experiment + '/' + base + '/' + network + '_nr_results'
                        nr_dic = load_run(filename)

                        filename = experiment + '/' + base + '/' + network + '_bs_results'
                        bs_dic = load_run(filename)

                        cplex_times[row, col] = cp_dic['solver_elapsed']
                        cplex_obj = cp_dic['solver_obj']

                        cplex_iter_objs = cp_dic['solver_primals']
                        cplex_duals = cp_dic['solver_duals']

                        cplex_gaps = 100*abs(np.array(cplex_iter_objs[:-1]) - np.array(cplex_duals))/np.minimum(np.array(cplex_iter_objs[:-1]),np.array(cplex_duals))

                        cplex_iter_times = cp_dic['solver_iter_times']
                        cplex_feasible = cp_dic['solver_feasible']

                        cplex_iter_objs[-1] = cplex_obj

                        nr_iter_objs = nr_dic['nr_iter_objs']
                        nr_iter_elapsed = nr_dic['nr_iter_elapsed']

                        bs_iter_objs = bs_dic['bs_iter_objs']
                        bs_iter_elapsed = bs_dic['bs_iter_elapsed']

                        lb_elapsed = bs_dic['elapsed_lower_bound']
                        ub_elapsed = bs_dic['elapsed_upper_bound']

                    nr_iter_elapsed = lb_elapsed + np.array(nr_iter_elapsed)
                    bs_iter_elapsed = lb_elapsed + \
                        ub_elapsed + np.array(bs_iter_elapsed)

                    # continue

                    # for nr_ind in range(len(nr_iter_objs) - 1, 0, -1):
                    #     if not am_i_close(nr_iter_objs[nr_ind], nr_iter_objs[nr_ind - 1], perc_rel_tol=1e-8):

                    #         nr_iter_objs = nr_iter_objs[:nr_ind + 1]
                    #         nr_iter_elapsed = nr_iter_elapsed[:nr_ind + 1]
                    #         break

                    # for bs_ind in range(len(bs_iter_objs) - 1, 0, -1):
                    #     if not am_i_close(bs_iter_objs[bs_ind], bs_iter_objs[bs_ind - 1], perc_rel_tol=1e-8):
                    #         bs_iter_objs = bs_iter_objs[:bs_ind + 1]
                    #         bs_iter_elapsed = bs_iter_elapsed[:bs_ind + 1]
                    #         break

                    for _ind in range(len(cplex_iter_objs) - 1, 0, -1):
                        if not am_i_close(cplex_iter_objs[_ind], cplex_iter_objs[_ind - 1], perc_rel_tol=1e-7):
                            cplex_iter_objs = cplex_iter_objs[:_ind + 1]
                            cplex_iter_times = cplex_iter_times[:_ind + 1]
                            break

                    # if atype == '_sr_' and exponent == '13' and (tail == 'a' or tail == 'e') and base == 'goto':
                    #     nr_times_g[tail] = nr_iter_elapsed
                    #     nr_objs_g[tail] = nr_iter_objs

                    #     bs_times_g[tail] = bs_iter_elapsed
                    #     bs_objs_g[tail] = bs_iter_objs

                    #     cplex_times_g[tail] = cplex_iter_times
                    #     cplex_objs_g[tail] = cplex_iter_objs

                    def get_indices(algo, iter_objs, common, common_2):

                        ind = None
                        trial = False
                        if trial:

                            if algo == 'cplex' and (iter_objs[0] < iter_objs[1]):
                                for p in range(len(iter_objs) - 1, 0, -1):
                                    if iter_objs[p] < common:
                                        ind = p
                                        break

                                    elif iter_objs[p] == common:
                                        ind = p
                                        break
                            else:
                                for p in range(len(iter_objs) - 1, 0, -1):
                                    if iter_objs[p] > common:
                                        ind = p + 1
                                        break

                                    elif iter_objs[p] == common:
                                        ind = p
                                        break

                        else:

                            if algo != 'cplex':
                                for p in range(len(iter_objs)):
                                    if am_i_close(iter_objs[p], common, perc_rel_tol=perc_rel_tol):
                                        ind = p
                                        break

                            else:

                                if iter_objs[0] > iter_objs[1]:
                                    for p in range(len(iter_objs)):
                                        if am_i_close(iter_objs[p], common, perc_rel_tol=perc_rel_tol):
                                            ind = p
                                            break

                                else:
                                    for p in range(len(iter_objs) - 1, 0, -1):
                                        if not am_i_close(iter_objs[p], common, perc_rel_tol=perc_rel_tol):
                                            ind = p
                                            break
                                    if p + 1 != len(iter_objs):
                                        p += 1

                        if ind is None:
                            ind = len(iter_objs) - 1
                        if ind == len(iter_objs):
                            ind = len(iter_objs) - 1

                        return ind

                    def find_indices(iter_objs, algo, common, common_2=None):
                        ind = get_indices(algo, iter_objs, common, common_2)

                        return ind

                    # if data_dic['solver_infeas'] < 1e-6:

                    #     order = [min(nr_iter_objs), min(bs_iter_objs), cplex_obj]
                    #     order = np.array(order)

                    #     common = max(order)
                    #     common_ind = np.argmax(order)

                    #     order = np.delete(order, common_ind)
                    #     common_2 = max(order)

                    #     nr_i = find_indices(nr_iter_objs, 'nr', common, common_2)
                    #     bs_i = find_indices(bs_iter_objs, 'bs', common, common_2)
                    #     cp_i = find_indices(
                    #         cplex_iter_objs, 'cplex', common, common_2)

                    # else:
                    #     common = max(min(nr_iter_objs), min(bs_iter_objs))

                    #     nr_i = find_indices(nr_iter_objs, 'nr', common, common_2)
                    #     bs_i = find_indices(bs_iter_objs, 'bs', common, common_2)
                    #     cp_i = find_indices(
                    #         cplex_iter_objs, 'cplex', common, common_2)

                    if cplex_iter_objs[0] > cplex_iter_objs[1]:
                        common = min(min(nr_iter_objs), min(bs_iter_objs))
                    else:
                        common = min(min(nr_iter_objs), min(bs_iter_objs), min(cplex_iter_objs))


                    nr_i = find_indices(nr_iter_objs, 'nr', common)
                    bs_i = find_indices(bs_iter_objs, 'bs', common)
                    cp_i = find_indices(cplex_iter_objs, 'cplex', common)

                    if nr_i is None:
                        nr_i = len(nr_iter_elapsed) - 1

                    if bs_i is None:
                        bs_i = len(bs_iter_elapsed) - 1

                    # print('*********')
                    # print('common: ', common)
                    # print('network: ', network)
                    # print('cplex iters: ', cplex_iter_objs)
                    # print('nr_iter_objs: ', nr_iter_objs)
                    # print('bs_iter_objs: ', bs_iter_objs)
                    # print('---------')
                    
                    print('nr_obj: ', nr_iter_objs[nr_i])
                    print('bs_obj: ', bs_iter_objs[bs_i])
                    print('cp_obj: ', cplex_iter_objs[cp_i])
                    # try:

                    #     print('cplex end obj: ', cplex_iter_objs[cp_i])
                    # except:
                    #     pdb.set_trace()
                    # print('cplex_gap: ', cplex_gaps[cp_i])
                    
                    # print('nr_tm: ', nr_iter_elapsed[nr_i])
                    # print('bs_tm: ', bs_iter_elapsed[bs_i])
                    # print('cplex end tm: ', cplex_iter_times[cp_i])
                    # print('&&&&&&&&&&&&')

                    cplex_times[row, col] = cplex_iter_times[cp_i]

                    nr_elapsed = nr_iter_elapsed[nr_i]
                    bs_elapsed = bs_iter_elapsed[bs_i]

                    bs_iter_list.append(bs_i)
                    nr_iter_list.append(nr_i)
                    cplex_infeas_list.append(data_dic['solver_infeas'])

                    nr_times[row, col] = nr_elapsed
                    bsc_times[row, col] = bs_elapsed

                    j += 1

                nr_iters_[row] = np.mean(nr_iter_list)
                bs_iters_[row] = np.mean(bs_iter_list)
                cplex_infeas_[row] = np.mean(cplex_infeas_list)

            family_type = atype.strip('_')
            if atype == '_8_':
                num_arcs = '8n'

            elif atype == '_lo_8_':
                num_arcs = '8n'

            elif atype == '_sr_':
                num_arcs = 'n' + r'$\sqrt{n}$'

            else:
                num_arcs = 'n' + r'$\sqrt{n}$'

            family_name = base.upper() + '-' + family_type.upper() + '  (m=' + num_arcs + ')'
            graph_family(family_name, cplex_times, nr_times,
                         bsc_times, nr_iters_, bs_iters_, cplex_infeas_)
        plt.tight_layout(pad=1.5)
        save_fig('comparison_' + base.lower())
        plt.clf()

    plt.figure(figsize=(16, 8))
    i = 0
    # GOTO INFEASIBILITY
    infeasiblity_tails = ['e']
    for tail in infeasiblity_tails:

        i += 1
        figure_grid = 120
        plt.subplot(figure_grid + i)

        # best_known = min(nr_objs_g[tail][nrmin], bs_objs_g[tail][bsmin])
        best_known = min(min(nr_objs_g[tail]), min(bs_objs_g[tail]))

        nrmin = np.argmin(nr_objs_g[tail])
        bsmin = np.argmin(bs_objs_g[tail])

        nr_objs_g[tail] = nr_objs_g[tail][:nrmin + 1]
        nr_times_g[tail] = nr_times_g[tail][:nrmin + 1]

        bs_objs_g[tail] = bs_objs_g[tail][:bsmin + 1]
        bs_times_g[tail] = bs_times_g[tail][:bsmin + 1]

        nr_objs_g[tail] = np.array(nr_objs_g[tail])
        bs_objs_g[tail] = np.array(bs_objs_g[tail])
        cplex_objs_g[tail] = np.array(cplex_objs_g[tail])

        nr_objs_g[tail] -= best_known
        bs_objs_g[tail] -= best_known
        cplex_objs_g[tail] -= best_known

        nr_objs_g[tail] /= best_known
        bs_objs_g[tail] /= best_known
        cplex_objs_g[tail] /= best_known

        plt.plot(nr_times_g[tail], nr_objs_g[tail], 'g:s', linewidth=1, label='NR')
        plt.plot(bs_times_g[tail], bs_objs_g[tail],
                 'r:^', linewidth=1, label='BSC')
        # plt.plot(np.cumsum(cplex_times_g[tail]), abs(cplex_objs_g[tail]), 'b:o', linewidth=1, label='CPLEX')
        plt.plot(cplex_times_g[tail], abs(cplex_objs_g[tail]),
                 'b:o', linewidth=1, label='CPLEX')

        plt.ylabel('Relative objective gap (logscale)', fontsize=12)
        plt.yscale('log')
        plt.ylim(bottom=1e-8)
        plt.xlabel('Running time', fontsize=12)
        # plt.axis([xs,xe,ys,ye])
        family_name = 'GOTO' + '-' + 'SR' + '-13' + \
            tail + '  (m=' + 'n' + r'$\sqrt{n}$' + ')'
        plt.grid(True)
        plt.title(family_name, fontsize=12)
        plt.legend(loc='upper right', fontsize=10)

    save_fig('infeasibility_' + 'sr')
    plt.clf()

def graph_lambar_experiments():


    experiment = 'varying_lambar'
    bases = ['netgen']
    tails = ['a','b','c','d','e','f','g','h','j','k']

    exponents = np.array(['12'])
    types = ['_lo_8_', '_8_', '_lo_sr_', '_sr_']

    lams = np.array([0.01, 10, 1000])


    # plt.figure(figsize=(16, 8))
    for base in bases:
        i = 0
        for atype in types:
            plt.figure(figsize=(16, 8))

            if base == 'goto' and atype.find('lo') >= 0:
                continue

            i += 1
            if base == 'netgen':
                elapsed = np.zeros((3, len(lams)))

                perc_rel_tol = 1e-2

                # figure_grid = 220
                # plt.subplot(figure_grid + i)
            else:

                perc_rel_tol = 1e-2

                figure_grid = 120
                plt.subplot(figure_grid + i)

            for lam in [0.01, 10, 1000]:

                lam_dir = str(lam).replace('.', '_')
                if lam_dir == '1_0':
                    lam_dir = '1'
                print(lam_dir)

                cplex_times = np.zeros((len(exponents), len(tails)))
                nr_times = np.zeros((len(exponents), len(tails)))
                bsc_times = np.zeros((len(exponents), len(tails)))

                j = 0
                for exponent in exponents:
                    for tail in tails:

                        filename = experiment + '/' + base + '/' + \
                            base + atype + exponent + tail + '_' + lam_dir
                        data_dic = load_run(filename)

                        lam_col = np.where(lams == lam)[0][0]

                        row = np.where(exponents == exponent)[0][0]
                        col = j

                        cplex_times[row, col] = data_dic['solver_elapsed']
                        cplex_obj = data_dic['solver_obj']

                        cplex_iter_objs = data_dic['solver_primals']
                        cplex_iter_times = data_dic['solver_iter_times']
                        cplex_feasible = data_dic['solver_feasible']
                        # cplex_iter_times_bk = data_dic['solver_iter_objs']
                        # cplex_iter_times = [data_dic['solver_elapsed'] - sum(cplex_iter_times_bk)] + cplex_iter_times

                        cplex_iter_objs[-1] = cplex_obj

                        nr_iter_objs = data_dic['nr_iter_objs']
                        nr_iter_elapsed = data_dic['nr_iter_elapsed']

                        bs_iter_objs = data_dic['bs_iter_objs']
                        bs_iter_elapsed = data_dic['bs_iter_elapsed']

                        lb_elapsed = data_dic['elapsed_lower_bound']
                        ub_elapsed = data_dic['elapsed_upper_bound']

                        nr_iter_elapsed = lb_elapsed + np.array(nr_iter_elapsed)
                        bs_iter_elapsed = lb_elapsed + \
                            ub_elapsed + np.array(bs_iter_elapsed)



                        def get_indices(algo, iter_objs, common, common_2):

                            ind = None
                            trial = False
                            if trial:

                                if algo == 'cplex' and (iter_objs[0] < iter_objs[1]):
                                    for p in range(len(iter_objs) - 1, 0, -1):
                                        if iter_objs[p] < common:
                                            ind = p
                                            break

                                        elif iter_objs[p] == common:
                                            ind = p
                                            break
                                else:
                                    for p in range(len(iter_objs) - 1, 0, -1):
                                        if iter_objs[p] > common:
                                            ind = p + 1
                                            break

                                        elif iter_objs[p] == common:
                                            ind = p
                                            break

                            else:

                                if algo != 'cplex':
                                    for p in range(len(iter_objs)):
                                        if am_i_close(iter_objs[p], common, perc_rel_tol=perc_rel_tol):
                                            ind = p
                                            break

                                else:

                                    if iter_objs[0] > iter_objs[1]:
                                        for p in range(len(iter_objs)):
                                            if am_i_close(iter_objs[p], common, perc_rel_tol=perc_rel_tol):
                                                ind = p
                                                break

                                    else:
                                        for p in range(len(iter_objs) - 1, 0, -1):
                                            if not am_i_close(iter_objs[p], common, perc_rel_tol=perc_rel_tol):
                                                ind = p
                                                break
                                        if p + 1 != len(iter_objs):
                                            p += 1

                            if ind is None:
                                ind = len(iter_objs) - 1
                            if ind == len(iter_objs):
                                ind = len(iter_objs) - 1

                            return ind

                        def find_indices(iter_objs, algo, common, common_2):
                            ind = get_indices(algo, iter_objs, common, common_2)

                            return ind

                        if data_dic['solver_infeas'] < 1e-6:

                            order = [min(nr_iter_objs), min(bs_iter_objs), cplex_obj]
                            order = np.array(order)

                            common = max(order)
                            common_ind = np.argmax(order)

                            order = np.delete(order, common_ind)
                            common_2 = max(order)

                            nr_i = find_indices(nr_iter_objs, 'nr', common, common_2)
                            bs_i = find_indices(bs_iter_objs, 'bs', common, common_2)
                            cp_i = find_indices(
                                cplex_iter_objs, 'cplex', common, common_2)

                        else:
                            common = max(min(nr_iter_objs), min(bs_iter_objs))

                            nr_i = find_indices(nr_iter_objs, 'nr', common, common_2)
                            bs_i = find_indices(bs_iter_objs, 'bs', common, common_2)
                            cp_i = find_indices(
                                cplex_iter_objs, 'cplex', common, common_2)


                        cplex_times[row, col] = cplex_iter_times[cp_i]

                        nr_elapsed = nr_iter_elapsed[nr_i]
                        bs_elapsed = bs_iter_elapsed[bs_i]

                        nr_times[row, col] = nr_elapsed
                        try:
                            bsc_times[row, col] = bs_elapsed
                        except:
                            pdb.set_trace()
                        j += 1

                elapsed[0, lam_col] = nr_times.mean(axis=1)
                elapsed[1, lam_col] = bsc_times.mean(axis=1)
                elapsed[2, lam_col] = cplex_times.mean(axis=1)

            barWidth = 0.2
            r_nr = np.arange(3)
            r_bsc = [x + barWidth for x in r_nr]
            r_cplex = [x + 2 * barWidth for x in r_nr]

            # plt.rcParams["font.size"] = 14

            plt.bar(r_cplex, elapsed[2, :], width=barWidth, hatch='xxx', edgecolor='#087efe', color='#087efe',
                    ecolor='#c6ccce', alpha=0.8,  capsize=5, label='CPLEX')

            plt.bar(r_bsc, elapsed[1, :], width=barWidth, hatch='\\\\\\', edgecolor='#FF4500', color='#FF4500',
                    ecolor='#c6ccce', alpha=0.8, capsize=5, label='BSC')

            plt.bar(r_nr, elapsed[0, :], width=barWidth, hatch='///', edgecolor='#b7fe00',
                    color='#b7fe00', ecolor='#c6ccce', alpha=0.8, capsize=5, label='NR')
            
            

            # tap_means_scaled[i]/2
            for pind in r_nr:
                
                
                plt.annotate('{0:1.1f}'.format(elapsed[2, pind]), (pind + 2 * barWidth, 0), textcoords='offset points', xytext=(
                    0, 20), ha='center', va='bottom', rotation=70,  size=20)
                plt.annotate('{0:1.1f}'.format(elapsed[1, pind]), (pind + barWidth, 0), textcoords='offset points', xytext=(
                    0, 20), ha='center', va='bottom', rotation=70,  size=20)
                plt.annotate('{0:1.1f}'.format(elapsed[0, pind]), (pind, 0), textcoords='offset points', xytext=(
                    0, 20), ha='center', va='bottom', rotation=70, size=20)

            plt.ylabel('Running time', fontsize=24)

            plt.xticks([(r + barWidth) for r in range(3)],
                       [r'$\bar\lambda$' + '=0.01', r'$\bar\lambda$' + '=10', r'$\bar\lambda$' + '=1000'], fontsize=20)
            plt.yticks(fontsize=20)
            family_type = atype.strip('_')
            if atype == '_8_':
                num_arcs = '8n'

            elif atype == '_lo_8_':
                num_arcs = '8n'

            elif atype == '_sr_':
                num_arcs = 'n' + r'$\sqrt{n}$'

            else:
                num_arcs = 'n' + r'$\sqrt{n}$'

            family_name = base.upper() + '-' + family_type.upper() + '  (m=' + num_arcs + ')'

            plt.title(family_name, fontsize=24)
            plt.tight_layout()
            # plt.figtext(0.5, 0.01, txt, wrap=True,
            #             ha='center', va="bottom", fontsize=7)
            plt.legend(fontsize=24)
        # plt.suptitle('Performance sensitivity to ' + r'$\lambda$')
            save_fig('varying_lambar_' + family_type.lower().replace('_',''), tight_layout=True)
            plt.clf()

# graph_gap_levels()
# graph_lambar_experiments()
# graph_base_vs_reliable()
graph_comparison()



