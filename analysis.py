import matplotlib.pyplot as plt
import pandas
import copy
import pdb
import pandas as pd
from utils import *
import numpy as np
# plt.plot([7813, 7813], [0., 0.9], "r:")         # Not shown
# plt.plot([-50000, 7813], [0.9, 0.9], "r:")      # Not shown
# plt.plot([-50000, 7813], [0.4368, 0.4368], "r:")# Not shown
# plt.plot([7813], [0.9], "ro")                   # Not shown
# plt.plot([7813], [0.4368], "ro")                # Not shown

def save_fig(fig_id, tight_layout=False, fig_extension="png", resolution=300):
    path = os.path.join(PLOTS_PATH, fig_id + "." + fig_extension)
    print("Saving figure", fig_id)

    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format=fig_extension, dpi=resolution)



def mean_std(array):
    mean = array.mean(axis=1)
    std = array.std(axis=1)
    print(mean)
    return np.log(mean), np.log(std)


def graph_family(family_name, cplex_times, nr_times, bsc_times, nr_iters_, bs_iters_, cplex_infeas_):
    node_num = [r'$2^{10}$', r'$2^{11}$', r'$2^{12}$', r'$2^{13}$']

    print(family_name)

    print('cplex')
    cplex_mean, cplex_std = mean_std(cplex_times)
    print('nr')
    nr_mean, nr_std = mean_std(nr_times)
    print('bsc')
    bsc_mean, bsc_std = mean_std(bsc_times)
    
    print(cplex_infeas_, nr_iters_, bs_iters_)

    plt.plot(node_num, cplex_mean, 'b:o', linewidth=2, label='CPLEX')
    plt.plot(node_num, nr_mean, 'g-s', linewidth=2, label='NR')
    plt.plot(node_num, bsc_mean, 'r--v', linewidth=2, label='BSC')

    # plt.fill_between(node_num, cplex_mean - cplex_std, cplex_mean + cplex_std, facecolor='blue', alpha=0.1);
    # plt.fill_between(node_num, nr_mean - nr_std, nr_mean + nr_std, facecolor='green', alpha=0.1);
    # plt.fill_between(node_num, bsc_mean - bsc_std, bsc_mean + bsc_std, facecolor='red', alpha=0.1);

    plt.ylabel('Running time (logscale)', fontsize=12)
    plt.xlabel('Number of nodes', fontsize=12)
    # plt.axis([xs,xe,ys,ye])

    plt.grid(True)
    plt.title(family_name, fontsize=12)
    plt.legend(loc='upper left', fontsize=10)



experiment = 'graph_families'

cplex_times = []
nr_times = []

bases = ['netgen', 'goto']
tails = ['a', 'b', 'c', 'd', 'e']
exponents = (np.arange(10, 14)).astype(str)
types = ['_8_', '_sr_', '_lo_8_', '_lo_sr_']



nr_times_g = {}
nr_objs_g = {}

bs_times_g = {}
bs_objs_g = {}

cplex_times_g = {}
cplex_objs_g = {}

nr_times_g_reg = {}
nr_objs_g_reg = {}

bs_times_g_reg = {}
bs_objs_g_reg = {}

cplex_times_g_reg = {}
cplex_objs_g_reg = {}

nr_times_e = {}
nr_objs_e = {}

bs_times_e = {}
bs_objs_e = {}

cplex_times_e = {}
cplex_objs_e = {}

for base in bases:
    plt.figure(figsize=(16,8))

    i = 0
    for atype in types:
        if base == 'goto' and atype.find('lo')>0:
            continue
        i += 1
        if base == 'netgen':
            rel_tol = 1e-8

            figure_grid = 220
            plt.subplot(figure_grid + i)
        else:
            rel_tol = 1e-8

            figure_grid = 120
            plt.subplot(figure_grid + i)

        cplex_times = np.zeros((len(exponents), len(tails)))
        nr_times = np.zeros((len(exponents), len(tails)))
        bsc_times = np.zeros((len(exponents), len(tails)))

        nr_iters_ = np.zeros(len(exponents))
        bs_iters_ = np.zeros(len(exponents))
        cplex_infeas_ = np.zeros(len(exponents))

        for exponent in exponents:
            bs_iter_list = []
            nr_iter_list = []
            cplex_infeas_list = []
            j = 0

            for tail in tails:
                cut_idx = None
                weird = False
                filename = experiment + '/' + base + '/' + base + atype + exponent + tail
                data_dic = load_run(filename)
                row = np.where(exponents==exponent)[0][0]
                col = j
                cplex_times[row, col] = data_dic['solver_elapsed']
                cplex_obj = data_dic['solver_obj']

                cplex_iter_objs = data_dic['solver_iter_elapsed']
                cplex_iter_times = data_dic['solver_iter_objs']

                cplex_iter_objs[-1] = cplex_obj

                nr_iter_objs = data_dic['nr_iter_objs']
                nr_iter_elapsed = data_dic['nr_iter_elapsed']

                bs_iter_objs = data_dic['bs_iter_objs']
                bs_iter_elapsed = data_dic['bs_iter_elapsed']

                lb_elapsed = data_dic['elapsed_lower_bound']
                ub_elapsed = data_dic['elapsed_upper_bound']

                nr_iter_elapsed = lb_elapsed + np.array(nr_iter_elapsed)
                bs_iter_elapsed = lb_elapsed + ub_elapsed + np.array(bs_iter_elapsed)

                if atype=='_sr_' and exponent=='12' and (tail=='a' or tail=='e') and base=='goto':
                    
                    nr_times_g[tail] = nr_iter_elapsed
                    nr_objs_g[tail] = nr_iter_objs

                    bs_times_g[tail] = bs_iter_elapsed
                    bs_objs_g[tail] = bs_iter_objs

                    cplex_times_g[tail] = [data_dic['solver_elapsed'] - sum(cplex_iter_times)]
                    cplex_times_g[tail] = cplex_times_g[tail] + cplex_iter_times
                    cplex_objs_g[tail] = cplex_iter_objs

                if atype=='_8_' and exponent=='12' and (tail=='a' or tail=='e') and base=='goto':
                    nr_times_g_reg[tail] = nr_iter_elapsed
                    nr_objs_g_reg[tail] = nr_iter_objs

                    bs_times_g_reg[tail] = bs_iter_elapsed
                    bs_objs_g_reg[tail] = bs_iter_objs

                    cplex_times_g_reg[tail] = [data_dic['solver_elapsed'] - sum(cplex_iter_times)]
                    cplex_times_g_reg[tail] = cplex_times_g_reg[tail] + cplex_iter_times
                    cplex_objs_g_reg[tail] = cplex_iter_objs

                def get_indices(rel_tol):
                    nr_i, bs_i = None, None

                    for p in range(len(nr_iter_objs)):
                        if math.isclose(cplex_obj, nr_iter_objs[p], rel_tol=rel_tol):
                            nr_i = p
                            break

                        if p==len(nr_iter_objs)-1:
                            weird = True


                    for p in range(len(bs_iter_objs)):
                        if math.isclose(cplex_obj, bs_iter_objs[p], rel_tol=rel_tol):
                            bs_i = p
                            break


                        if p==len(bs_iter_objs)-1:
                            weird = True

                    return nr_i, bs_i


                nr_i, bs_i = get_indices(rel_tol)


                if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol):
                    nr_i, bs_i = get_indices(rel_tol*10)

                    if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol*10):
                        print(base)
                        print('weird case')
                        print('solver infeas: ', print(data_dic['solver_infeas']))

                        nr_i, bs_i = get_indices(rel_tol*100)

                        if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol*100):

                            nr_i, bs_i = get_indices(rel_tol*10000)
                        

                # print(data_dic['solver_infeas'])
                # print('cplex objs:')
                # print(cplex_iter_objs)
                # print('cplex times:')
                # print(cplex_iter_times)
                # print('nr objs:')
                # print(nr_iter_objs)
                # print('nr times:')
                # print(nr_iter_elapsed)
                # print('bs objs:')
                # print(bs_iter_objs)
                # print('bs times:')
                # print(bs_iter_elapsed)

                # print('nr_i', nr_i)
                # print('bs_i', bs_i)

                # for idx in range(len(cplex_iter_objs)):
                #     if math.isclose(cplex_iter_objs[idx], nr_iter_objs[nr_i], rel_tol=rel_tol):
                #         cut_idx = idx
                #         if len(cplex_iter_objs) - cut_idx > 5:
                #                 pdb.set_trace()
                #                 print('cut idx found')
                #         break

                # if cut_idx == None:
                #     for idx in range(len(cplex_iter_objs)):
                #         if math.isclose(cplex_iter_objs[idx], nr_iter_objs[nr_i], rel_tol=rel_tol*10):
                #             cut_idx = idx
                #             if len(cplex_iter_objs) - cut_idx > 5:
                #                 print('cut idx found')
                #             break

                # if cut_idx == None:
                #     cut_idx = len(cplex_iter_objs)-1

                cplex_times[row, col] = data_dic['solver_elapsed'] #- sum(cplex_iter_times[cut_idx:])


                if atype=='_sr_' and exponent=='12' and (tail=='a' or tail=='e') and base=='netgen':
                    print('convergence graph work')

                    nr_times_e[tail] = nr_iter_elapsed[:nr_i+1]
                    nr_objs_e[tail] = nr_iter_objs[:nr_i+1]

                    bs_times_e[tail] = bs_iter_elapsed[:bs_i+1]
                    bs_objs_e[tail] = bs_iter_objs[:bs_i+1]

                    cplex_times_e[tail] = [data_dic['solver_elapsed'] - sum(cplex_iter_times)]
                    cplex_times_e[tail] = cplex_times_e[tail] + cplex_iter_times
                    cplex_objs_e[tail] = cplex_iter_objs
                    # cplex_times_e = cplex_times_e[:cut_idx+1]
                    # cplex_objs_e = cplex_iter_objs[:cut_idx+1]

                nr_elapsed = nr_iter_elapsed[nr_i]
                bs_elapsed = bs_iter_elapsed[bs_i]

                bs_iter_list.append(bs_i)
                nr_iter_list.append(nr_i)
                cplex_infeas_list.append(data_dic['solver_infeas'])

                nr_times[row, col] = nr_elapsed
                bsc_times[row, col] = bs_elapsed

                # if bs_elapsed < nr_elapsed:
                #     pdb.set_trace()

                # if weird:
                #     print('instance: ', filename)
                #     print('nr_iterations: \n', nr_iter_objs)
                #     print('bs_iterations: \n', bs_iter_objs)
                #     print('cplex_iterations: \n', cplex_iter_objs)
                #     print('================')
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
        graph_family(family_name, cplex_times, nr_times, bsc_times, nr_iters_, bs_iters_, cplex_infeas_)
    # plt.subplots.adjust(hspace = 0.2)
    plt.tight_layout(pad=1.5)

    if base == 'netgen':
        plt.plot([0.5, 0.5], [0, 1], marker=',', color='lightgreen', lw=1 ,transform=plt.gcf().transFigure, clip_on=False)
        plt.plot([0, 1], [0.5, 0.5], marker=',', color='lightgreen', lw=1 ,transform=plt.gcf().transFigure, clip_on=False)
    
    save_fig('comparison_' + base.lower())
    plt.clf()


plt.figure(figsize=(16,8))
### NETGEN convergence
convergence_tails = ['a', 'e']
i=0
for tail in convergence_tails:
    i +=1
    figure_grid = 120
    plt.subplot(figure_grid + i)

    best_known = min(min(nr_objs_e[tail]), min(bs_objs_e[tail]))
    nr_objs_e[tail][-1] = best_known
    bs_objs_e[tail][-1] = best_known
    cplex_objs_e[tail][-1] = best_known


    ref = cplex_objs_e[tail][-1] 
    round_amount = 8 - len(str(cplex_objs_e[tail][-1])[:str(cplex_objs_e[tail][-1]).find('.')])
    cind_it = 0
    for cind in range(len(cplex_objs_e[tail]))[::-1]:
        cind_it += 1
        if round(ref, round_amount) == cplex_objs_e[tail][cind]:
            cplex_objs_e[tail][cind] = ref
        if cind_it > 5:
            break

    nr_objs_e[tail] -= best_known
    bs_objs_e[tail] -= best_known
    cplex_objs_e[tail] -= best_known


    nr_objs_e[tail] = np.array(nr_objs_e[tail])
    indice = np.where(nr_objs_e[tail]==0.0)[0][0]
    nr_objs_e[tail][indice] = 1

    bs_objs_e[tail] = np.array(bs_objs_e[tail])
    indice = np.where(bs_objs_e[tail]==0.0)[0][0]
    bs_objs_e[tail][indice] = 1

    cplex_objs_e[tail] = np.array(cplex_objs_e[tail])
    indice = np.where(cplex_objs_e[tail]==0.0)[0][0]
    cplex_objs_e[tail][indice] = 1

    plt.plot(nr_times_e[tail], np.log(nr_objs_e[tail]), 'g-s', linewidth=2, label='NR')
    plt.plot(bs_times_e[tail], np.log(bs_objs_e[tail]), 'r--v', linewidth=2, label='BSC')
    plt.plot(np.cumsum(cplex_times_e[tail]), np.log(cplex_objs_e[tail]), 'b:o', linewidth=2, label='CPLEX')

    plt.ylabel('Objective gap (logscale)', fontsize=12)
    plt.xlabel('Running time', fontsize=12)
    # plt.axis([xs,xe,ys,ye])
    family_name = 'NETGEN' + '-' + 'SR' + '-12' + tail + '  (m=' + 'n' + r'$\sqrt{n}$' + ')' 
    plt.grid(True)
    plt.title(family_name, fontsize=12)
    plt.legend(loc='upper right', fontsize=10)

save_fig('convergence_netgen')
plt.clf()


plt.figure(figsize=(16,8))
i=0
#### GOTO INFEASIBILITY
infeasiblity_tails = ['a', 'e']
for tail in infeasiblity_tails:

    i += 1
    figure_grid = 120
    plt.subplot(figure_grid + i)

    nrmin = np.argmin(nr_objs_g[tail])
    bsmin = np.argmin(bs_objs_g[tail])
    best_known = min(nr_objs_g[tail][nrmin], bs_objs_g[tail][bsmin])

    nr_objs_g[tail] = nr_objs_g[tail][:nrmin+1]
    nr_times_g[tail] = nr_times_g[tail][:nrmin+1]
    nr_objs_g[tail][-1] = best_known

    bs_objs_g[tail] = bs_objs_g[tail][:bsmin+1]
    bs_times_g[tail] = bs_times_g[tail][:bsmin+1]
    bs_objs_g[tail][-1] = best_known

    nr_objs_g[tail] -= best_known
    bs_objs_g[tail] -= best_known

    nr_objs_g[tail] = np.array(nr_objs_g[tail])
    indice = np.where(nr_objs_g[tail]==0.0)[0][0]
    nr_objs_g[tail][indice] = 1

    bs_objs_g[tail] = np.array(bs_objs_g[tail])
    indice = np.where(bs_objs_g[tail]==0.0)[0][0]
    bs_objs_g[tail][indice] = 1


    cplex_objs_g[tail] -= best_known
    cplex_objs_g[tail] = np.array(cplex_objs_g[tail])
    cplex_scale = np.log(abs(cplex_objs_g[tail])) * np.sign(cplex_objs_g[tail])

    plt.plot(nr_times_g[tail], np.log(nr_objs_g[tail]), 'g-s', linewidth=2, label='NR')
    plt.plot(bs_times_g[tail], np.log(bs_objs_g[tail]), 'r--v', linewidth=2, label='BSC')
    plt.plot(np.cumsum(cplex_times_g[tail]), cplex_scale, 'b:o', linewidth=2, label='CPLEX')

    plt.ylabel('Objective gap (logscale)', fontsize=12)
    plt.xlabel('Running time', fontsize=12)
    # plt.axis([xs,xe,ys,ye])
    family_name = 'GOTO' + '-' + 'SR' + '-12' + tail + '  (m=' + 'n' + r'$\sqrt{n}$' + ')' 
    plt.grid(True)
    plt.title(family_name, fontsize=12)
    plt.legend(loc='upper right', fontsize=10)
    
save_fig('infeasibility_' + 'sr')
plt.clf()

plt.figure(figsize=(16,8))
i=0
#### GOTO INFEASIBILITY
infeasiblity_tails = ['a', 'e']
for tail in infeasiblity_tails:

    i += 1
    figure_grid = 120
    plt.subplot(figure_grid + i)

    nrmin = np.argmin(nr_objs_g_reg[tail])
    bsmin = np.argmin(bs_objs_g_reg[tail])
    best_known = min(nr_objs_g_reg[tail][nrmin], bs_objs_g_reg[tail][bsmin])

    nr_objs_g_reg[tail] = nr_objs_g_reg[tail][:nrmin+1]
    nr_times_g_reg[tail] = nr_times_g_reg[tail][:nrmin+1]
    nr_objs_g_reg[tail][-1] = best_known

    bs_objs_g_reg[tail] = bs_objs_g_reg[tail][:bsmin+1]
    bs_times_g_reg[tail] = bs_times_g_reg[tail][:bsmin+1]
    bs_objs_g_reg[tail][-1] = best_known

    nr_objs_g_reg[tail] -= best_known
    bs_objs_g_reg[tail] -= best_known

    nr_objs_g_reg[tail] = np.array(nr_objs_g_reg[tail])
    indice = np.where(nr_objs_g_reg[tail]==0.0)[0][0]
    nr_objs_g_reg[tail][indice] = 1

    bs_objs_g_reg[tail] = np.array(bs_objs_g_reg[tail])
    indice = np.where(bs_objs_g_reg[tail]==0.0)[0][0]
    bs_objs_g_reg[tail][indice] = 1


    cplex_objs_g_reg[tail] -= best_known
    cplex_objs_g_reg[tail] = np.array(cplex_objs_g_reg[tail])
    cplex_scale = np.log(abs(cplex_objs_g_reg[tail])) * np.sign(cplex_objs_g_reg[tail])


    plt.plot(nr_times_g_reg[tail], np.log(nr_objs_g_reg[tail]), 'g-s', linewidth=2, label='NR')
    plt.plot(bs_times_g_reg[tail], np.log(bs_objs_g_reg[tail]), 'r--v', linewidth=2, label='BSC')
    plt.plot(np.cumsum(cplex_times_g_reg[tail]), cplex_scale, 'b:o', linewidth=2, label='CPLEX')

    plt.ylabel('Objective gap (logscale)', fontsize=12)
    plt.xlabel('Running time', fontsize=12)
    # plt.axis([xs,xe,ys,ye])
    family_name = 'GOTO' + '-' + '8' + '-12' + tail + '  (m=' + 'n' + r'$\sqrt{n}$' + ')' 
    plt.grid(True)
    plt.title(family_name, fontsize=12)
    plt.legend(loc='upper right', fontsize=10)

save_fig('infeasibility_' + '8')
plt.clf()



pdb.set_trace()

experiment = 'varying_lambar'
bases = ['netgen']
tails = ['a', 'b', 'c', 'd', 'e']
exponents = np.array(['11'])
types = ['_8_', '_sr_', '_lo_8_', '_lo_sr_']

lams = np.array([0.001, 0.1, 10, 1000])


plt.figure(figsize=(16,8))
for base in bases:
    i = 0
    for atype in types:
        if base == 'goto' and atype.find('lo')>=0:
            continue

        i+=1
        if base == 'netgen':
            elapsed = np.zeros((3,4))

            rel_tol = 1e-8

            figure_grid = 220
            plt.subplot(figure_grid + i)
        else:

            rel_tol = 1e-8

            figure_grid = 120
            plt.subplot(figure_grid + i)

        for lam in lams:

            if lam>=10:
                lam=int(lam)
            lam_dir = str(lam).replace('.','')

            print(lam_dir)

            cplex_times = np.zeros((len(exponents), len(tails)))
            nr_times = np.zeros((len(exponents), len(tails)))
            bsc_times = np.zeros((len(exponents), len(tails)))


            j = 0
            for exponent in exponents:
                for tail in tails:

                    filename = experiment + '/' + base + '/' + base + atype + exponent + tail +'_'+lam_dir
                    data_dic = load_run(filename)
                    
                    lam_col = np.where(lams==lam)[0][0]
                    
                    row = np.where(exponents==exponent)[0][0]
                    col = j
                    print(filename)

                    cplex_times[row, col] = data_dic['solver_elapsed']
                    cplex_obj = data_dic['solver_obj']

                    cplex_iter_objs = data_dic['solver_iter_elapsed']
                    cplex_iter_times = data_dic['solver_iter_objs']

                    cplex_iter_objs[-1] = cplex_obj

                    nr_iter_objs = data_dic['nr_iter_objs']
                    nr_iter_elapsed = data_dic['nr_iter_elapsed']

                    bs_iter_objs = data_dic['bs_iter_objs']
                    bs_iter_elapsed = data_dic['bs_iter_elapsed']

                    lb_elapsed = data_dic['elapsed_lower_bound']
                    ub_elapsed = data_dic['elapsed_upper_bound']

                    nr_iter_elapsed = lb_elapsed + np.array(nr_iter_elapsed)
                    bs_iter_elapsed = lb_elapsed + ub_elapsed + np.array(bs_iter_elapsed)



                    def get_indices(rel_tol):
                        nr_i, bs_i = None, None

                        for p in range(len(nr_iter_objs)):
                            if math.isclose(cplex_obj, nr_iter_objs[p], rel_tol=rel_tol):
                                nr_i = p
                                break

                            if p==len(nr_iter_objs)-1:
                                weird = True


                        for p in range(len(bs_iter_objs)):
                            if math.isclose(cplex_obj, bs_iter_objs[p], rel_tol=rel_tol):
                                bs_i = p
                                break


                            if p==len(bs_iter_objs)-1:
                                weird = True

                        return nr_i, bs_i


                    nr_i, bs_i = get_indices(rel_tol)

                    if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol):
                        nr_i, bs_i = get_indices(rel_tol*10)

                        if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol*10):
                            print(base)
                            print('weird case')
                            print('solver infeas: ', print(data_dic['solver_infeas']))

                            nr_i, bs_i = get_indices(rel_tol*100)

                            if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol*100):

                                nr_i, bs_i = get_indices(rel_tol*1000)

                                if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol*1000):

                                    nr_i, bs_i = get_indices(rel_tol*10000)
                                    
                                    if not math.isclose(cplex_obj, min(nr_iter_objs), rel_tol=rel_tol*10000):
                                    
                                        nr_i, bs_i = get_indices(rel_tol*100000)
                                        pdb.set_trace()


                    # for idx in range(len(cplex_iter_objs)):
                    #     if math.isclose(cplex_iter_objs[idx], nr_iter_objs[nr_i], rel_tol=rel_tol):
                    #         cut_idx = idx
                    #         break

                    # if cut_idx == None:
                    #     for idx in range(len(cplex_iter_objs)):
                    #         if math.isclose(cplex_iter_objs[idx], nr_iter_objs[nr_i], rel_tol=sec_rel_tol):
                    #             cut_idx = idx
                    #             break

                    # if cut_idx == None:
                    #     cut_idx = len(cplex_iter_objs)-1

                    cplex_time = data_dic['solver_elapsed'] #- sum(cplex_iter_times[cut_idx:])
                    cplex_times[row, col] = cplex_time
                    
                    nr_elapsed = nr_iter_elapsed[nr_i]
                    bs_elapsed = bs_iter_elapsed[bs_i]

                    nr_times[row, col] = nr_elapsed
                    try:
                        bsc_times[row, col] = bs_elapsed
                    except:
                        pdb.set_trace()
                    j+=1
            elapsed[0, lam_col] = nr_times.mean(axis=1)
            elapsed[1, lam_col] = bsc_times.mean(axis=1)
            elapsed[2, lam_col] = cplex_times.mean(axis=1) 


        barWidth = 0.2
        r_nr = np.arange(4)
        r_bsc = [x + barWidth for x in r_nr]
        r_cplex = [x + 2 * barWidth for x in r_nr]


        plt.rcParams["font.size"] = 8
        plt.bar(r_nr, elapsed[0,:], width=barWidth, hatch='///', edgecolor='#b7fe00', color='#b7fe00',ecolor='#c6ccce', alpha=0.8, capsize=5, label='NR')
        plt.bar(r_bsc, elapsed[1,:], width=barWidth, hatch='\\\\\\', edgecolor='#FF4500', color='#FF4500',
                ecolor='#c6ccce', alpha=0.8, capsize=5, label='BSC')
        plt.bar(r_cplex, elapsed[2,:], width=barWidth, hatch='xxx', edgecolor='#087efe', color='#087efe',
                ecolor='#c6ccce', alpha=0.8,  capsize=5, label='CPLEX')


        # tap_means_scaled[i]/2
        for pind in r_nr:
            plt.annotate('{0:1.1f}'.format(elapsed[0,pind]), (pind, 0), textcoords='offset points', xytext=(
                0, 20), ha='center', va='bottom', rotation=70, size=8)
            plt.annotate('{0:1.1f}'.format(elapsed[1,pind]), (pind + barWidth, 0), textcoords='offset points', xytext=(
                0, 20), ha='center', va='bottom', rotation=70,  size='smaller')
            plt.annotate('{0:1.1f}'.format(elapsed[2,pind]), (pind + 2 * barWidth, 0), textcoords='offset points', xytext=(
                0, 20), ha='center', va='bottom', rotation=70,  size='smaller')


        plt.ylabel('Running time (logscale)')

        plt.xticks([(r + barWidth) for r in range(4)],
                   [r'$\lambda$' + '=0.001', r'$\lambda$' + '=0.1', r'$\lambda$' + '=10', r'$\lambda$' + '=1000'])
       

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

        plt.title(family_name, fontsize=7)

        # plt.figtext(0.5, 0.01, txt, wrap=True,
        #             ha='center', va="bottom", fontsize=7)
        plt.legend(fontsize=8)
    # plt.suptitle('Performance sensitivity to ' + r'$\lambda$')
    save_fig('varying_lambar', tight_layout=True)

















