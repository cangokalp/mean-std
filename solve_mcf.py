from utils import *
from MCF_DiGraph import *

# import matplotlib.pyplot as plt
import pandas
import random
import itertools
import copy
# from tabulate import tabulate
import cProfile
import pstats
# sys.path.append('/usr/local/Cellar/graph-tool/2.27_1/lib/python3.7/site-packages/')
# import graph_tool as gt

import cplex
from prettytable import PrettyTable
import caffeine


def cvxpy_any(G, solver='ECOS'):

    start = time.time()
    x = cp.Variable(G.m)
    constraints = [0 <= x, x <= G.cap, G.A@x == G.b]
    objective = cp.Minimize(G.mu.T * x + G.lambar*cp.norm(cp.multiply(np.sqrt(G.var), x), 2))
    prob = cp.Problem(objective, constraints)

    # result = prob.solve(solver=solver, verbose=True)  # gurobi mosek compare

    result = prob.solve(solver='MOSEK', verbose=True, mosek_params={'MSK_DPAR_INTPNT_CO_TOL_REL_GAP':1e-12,'MSK_DPAR_INTPNT_CO_TOL_INFEAS':1e-12,'MSK_IPAR_INTPNT_MAX_ITERATIONS':1000})
    print(objective.value)
    
    cvx_soln = x.value
    cvx_obj = objective.value

    # print(prob.solver_stats.solve_time, prob.solver_stats.setup_time)
    cvx_elapsed = prob.solver_stats.solve_time
    elapsed = time.time() - start

    return cvx_obj, elapsed, cvx_soln


def cvxpy_solve(G):
    start = time.time()

    prob = cplex.Cplex()

    x_names = ['x' + str(i) for i in range(G.m)]
    lin_obj = G.mu
    prob.variables.add(obj=lin_obj, lb=np.zeros(G.m), ub=G.cap, names=x_names)

    prob.linear_constraints.add(rhs=G.b, senses='E' * G.n)
    prob.linear_constraints.set_coefficients(zip(G.rows, G.cols, G.values))

    prob.objective.set_sense(prob.objective.sense.minimize)

    decoy_name = ['decoy']
    decoyvar = prob.variables.add(obj=[G.lambar], lb=[0], names=decoy_name)

    qc_names_prepped = decoy_name + x_names
    qc_rhs_prepped = [-1.0]
    qc_rhs_prepped.extend(G.var)

    prob.quadratic_constraints.add(quad_expr=[qc_names_prepped,
                                              qc_names_prepped,
                                              qc_rhs_prepped],
                                   sense='L',
                                   rhs=0,
                                   name='q1')
    cplex_log = "cplex_log.txt"
    prob.set_log_stream(cplex_log)
    prob.set_error_stream(cplex_log)
    prob.set_warning_stream(cplex_log)
    prob.set_results_stream(cplex_log)
    prob.parameters.barrier.display.set(2)
    # prob.parameters.barrier.qcpconvergetol.set(1e-9)
    # prob.parameters.barrier.algorithm.set(1)
    # prob.parameters.barrier.limits.corrections.set(50)

    prob.solve()

    obj = prob.solution.get_objective_value()
    soln = np.array(prob.solution.get_values())

    m = prob.solution.quality_metric
    max_x, max_infeas = prob.solution.get_float_quality(
        [m.max_x, m.max_primal_infeasibility])

    elapsed = time.time() - start

    return obj, elapsed, soln, max_infeas


def cvxpy_solve_additive(G, lam, prob=None, warm_start=False, lp=False, solver='MOSEK', bound_lam=False):

    start = time.time()

    if not warm_start:
        prob = cplex.Cplex()
        if lp:
            if bound_lam:
                obj = np.sqrt(G.var)
            else:
                obj = (G.mu + lam)
            prob.variables.add(obj=obj, lb=np.zeros(G.m), ub=G.cap)
        else:
            prob.variables.add(lb=np.zeros(G.m), ub=G.cap)

            if not bound_lam:
                prob.objective.set_linear([(int(i), j)
                                           for i, j in zip(np.arange(G.m), G.mu)])
                var = lam * np.array(G.var) * 2.0

            else:
                var = np.array(G.var)

            prob.objective.set_quadratic(var)

        prob.linear_constraints.add(rhs=G.b, senses='E' * G.n)
        prob.linear_constraints.set_coefficients(zip(G.rows, G.cols, G.values))

        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.set_log_stream(None)
        prob.set_error_stream(None)
        prob.set_warning_stream(None)
        prob.set_results_stream(None)
        prob.parameters.barrier.display.set(0)

    else:
        if lp:
            if bound_lam:
                coef = G.var
            else:
                coef = (G.mu + lam)
            prob.objective.set_linear([(int(i), j)
                                       for i, j in zip(np.arange(G.m), coef)])
        else:
            var = lam * np.array(G.var) * 2.0
            prob.objective.set_quadratic(var)

    prob.solve()

    obj = prob.solution.get_objective_value()
    soln = np.array(prob.solution.get_values())

    elapsed = time.time() - start

    m = prob.solution.quality_metric
    max_x, max_infeas = prob.solution.get_float_quality(
        [m.max_x, m.max_primal_infeasibility])

    return obj, elapsed, soln, prob, max_infeas


def cvxpy_solve_xi(G, soln, lam, prob=None, warm_start=False, lp=False, solver='MOSEK', iters=0):

    start = time.time()

    diff = soln - G.cap
    x_ = soln
    x_zero = np.argwhere(abs(soln) < 1e-6).ravel().astype(int)
    x_nonzero = np.argwhere(soln > 1e-6).ravel().astype(int)
    x_u = np.argwhere(abs(diff) < 1e-6).ravel().astype(int)
    x_btw = np.array(list(set(x_nonzero).difference(set(x_u)))).astype(int)

    diff_ = G.cap[x_btw] - x_[x_btw]
    b = np.r_[np.zeros(len(x_zero)), np.zeros(
        len(x_u)), -x_[x_u], diff_, -x_[x_btw]]
    names = [str(i) + 'S' for i in range(len(b))]

    cols = np.r_[x_zero, x_u, x_u, x_btw, x_btw]
    senses = np.r_[['E'] * len(x_zero), ['L'] * len(x_u), ['G']
                   * len(x_u), ['L'] * len(x_btw), ['G'] * len(x_btw)]
    rows = np.arange(len(cols))
    values = np.ones(len(rows))

    if not warm_start:

        prob = cplex.Cplex()
        if lp:
            obj = np.ones(G.m)
            prob.variables.add(obj=obj, lb=-np.inf * np.ones(G.m), ub=G.cap)
        else:
            prob.variables.add(lb=-np.inf * np.ones(G.m), ub=G.cap)
            lin_coeffs = 2 * G.var * x_
            quad_coeffs = lam * G.var * 2
            prob.objective.set_linear(
                [(int(i), j) for i, j in zip(np.arange(G.m), lin_coeffs)])
            prob.objective.set_quadratic(quad_coeffs)

        prob.linear_constraints.add(rhs=np.zeros(G.n), senses='E' * G.n)
        prob.linear_constraints.set_coefficients(zip(G.rows, G.cols, G.values))

        prob.linear_constraints.add(rhs=b, names=names, senses=senses)
        prob.linear_constraints.set_coefficients(
            zip([int(i) for i in rows], [int(i) for i in cols], values))

        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.set_log_stream(None)
        prob.set_error_stream(None)
        prob.set_warning_stream(None)
        prob.set_results_stream(None)
        prob.parameters.barrier.display.set(0)
    else:

        prob.linear_constraints.delete(
            G.n, prob.linear_constraints.get_num() - 1)
        prob.linear_constraints.add(rhs=b, names=names, senses=senses)
        prob.linear_constraints.set_coefficients(
            zip([int(i) for i in rows], [int(i) for i in cols], values))

        if not lp:
            lin_coeffs = 2 * G.var * x_
            quad_coeffs = lam * G.var * 2
            prob.objective.set_linear(
                [(int(i), j) for i, j in zip(np.arange(G.m), lin_coeffs)])
            prob.objective.set_quadratic(quad_coeffs)

    prob.solve()

    obj = prob.solution.get_objective_value()
    soln = np.array(prob.solution.get_values())

    elapsed = time.time() - start

    m = prob.solution.quality_metric
    max_x, max_infeas = prob.solution.get_float_quality(
        [m.max_x, m.max_primal_infeasibility])

    return elapsed, soln, prob, max_infeas


def bs_cvxpy(G, low=0, high=1000, prob=None, lp=False, test=False, solver_weird=False):
    print('bsc')

    # test
    test = test
    if test:
        warm_start = False
        for i in np.linspace(low, high / 2, 300):
            lam = i
            obj_, elapsed, x, prob = cvxpy_solve_additive(
                G, lam, prob=prob, warm_start=warm_start, lp=lp, bound_lam=False)
            var_cost = np.multiply(G.var, x).dot(x)
            obj = G.mu.dot(x) + G.lambar * (np.sqrt(var_cost))

            f1 = lam - G.lambar / (2.0 * np.sqrt(var_cost))
            f2 = lam - G.lambar / np.sqrt(var_cost)

            t = PrettyTable(['lam_obj', 'lambar_obj', 'time',
                             'derivative_eq', 'func_eq', 'lambda'])
            t.add_row([obj_, obj, time.time() - start, f1, f2, lam])
            print(t)

            if abs(f1) < 1e-1 or abs(f2) < 1e-1:
                print('now debug')
                print(f1, f2, obj_, obj)
                pdb.set_trace()

        pdb.set_trace()
    start = time.time()

    stop_tol = 1e-6

    f = 100
    iters = 0
    found = False

    if lp:
        if prob is None:
            high = np.ones(G.m) * G.lambar
            low = np.ones(G.m) * low

    if prob is not None:
        warm_start = True
    else:
        warm_start = False

    mid = 10000
    iter_objs = []
    iter_elapsed = []

    infeas = []
    lams = []
    fs = []
    while not found:
        iters += 1
        mid_prev = mid
        mid = (high + low) / 2.0

        if iters > 1 and warm_start == False:
            warm_start = True

        lams.append(mid)
        obj, elapsed, x, prob, max_infeas = cvxpy_solve_additive(
            G, mid, prob=prob, warm_start=warm_start, lp=lp)

        infeas.append(max_infeas)
        var_cost = np.multiply(G.var, x).dot(x)
        obj = G.mu.dot(x) + G.lambar * (np.sqrt(var_cost))

        # soln_diff = np.linalg.norm(solver_soln - x)
        print(obj, time.time() - start, f)

        iter_objs.append(obj)
        iter_elapsed.append(time.time() - start)

        if lp:
            f = mid - G.lambar * np.multiply(G.var, x) / np.sqrt(var_cost)
        else:
            f = mid - G.lambar / (2.0 * np.sqrt(var_cost))

        # if lp:
        #     if np.all(f) < stop_tol or iters > 9:
        #         found = True
        #         break
        # else:
        #     fs.append(f)
        #     if solver_weird:
        #         if abs(f) <= stop_tol:
        #             pdb.set_trace()
        #             found = True
        #             break
        #     else:
        #         if abs(f) < stop_tol:
        #             found = True
        #             break

        if lp:
            pos = np.argwhere(np.sign(f) == np.sign(1))
            neg = np.argwhere(np.sign(f) == -np.sign(1))

            high[pos] = mid[pos]
            low[neg] = mid[neg]
        else:
            if np.sign(f) == np.sign(1):
                high = mid
            else:
                low = mid
        if iters > 9:
            break

    elapsed = time.time() - start
    infeas = np.array(infeas)
    return obj, elapsed, x, iter_objs, iter_elapsed, infeas.mean(), lams, fs


def nr_cvxpy(G, low=0, high=1000, prob=None, lp=False, lam_init=None, solver_weird=False):
    print('nr')
    start = time.time()
    stop_tol = 1e-6

    if lp:
        if prob is None:
            high = np.ones(G.m) * G.lambar
            low = np.ones(G.m) * low

    lam = (high + low) / 2.0

    if lam_init is not None:
        lam = lam_init

    found = False
    f = 100
    iters = 0
    warm_start = False
    prob = None
    warm_start_xi = False
    prob_xi = None

    iter_elapsed = []
    iter_objs = []

    iter_var_times = []
    iter_xi_times = []

    lam_prev = 1000
    mean = []
    var = []

    infeas = []
    fs = []
    lams = []
    while not found:
        iters += 1

        if iters > 1:
            warm_start = True
            warm_start_xi = True

        lams.append(lam)

        obj, elapsed, x, prob, max_infeas = cvxpy_solve_additive(
            G, lam, prob=prob, warm_start=warm_start, lp=lp)
        iter_var_times.append(elapsed)
        infeas.append(max_infeas)
        var_cost = np.multiply(G.var, x).dot(x)

        cur_mean = G.mu.dot(x)
        obj = cur_mean + G.lambar * (np.sqrt(var_cost))
        iter_objs.append(obj)
        iter_elapsed.append(time.time() - start)

        mean.append(cur_mean)
        var.append(var_cost)

        # soln_diff = np.linalg.norm(x - solver_soln)
        print(obj, time.time() - start, f)  # , soln_diff)

        if lp:
            f = lam - G.lambar * np.multiply(G.var, x) / np.sqrt(var_cost)
        else:
            f = lam - G.lambar / (2.0 * np.sqrt(var_cost))

        # if lp:
        #     if np.all(f) < stop_tol or iters > 5:
        #         found = True
        #         break
        # else:
        #     fs.append(f)
        #     if solver_weird:
        #         if abs(f) <= stop_tol:
        #             found = True
        #             break
        #     else:
        #         if abs(f) < stop_tol:
        #             found = True
        #             break

        elapsed_xi, xi, prob_xi, max_infeas = cvxpy_solve_xi(
            G, x, lam, prob=prob_xi, warm_start=warm_start_xi, lp=lp, iters=iters)
        iter_xi_times.append(elapsed_xi)
        infeas.append(max_infeas)

        if lp:
            var_x = np.multiply(G.var, x)
            f_lam_der = 1.0 - ((G.lambar * np.multiply(G.var, xi)) / (np.sqrt(var_cost)) - (
                (G.lambar * np.multiply(var_x.dot(var_x), xi)) / (var_cost**(3.0 / 2.0))))
        else:
            xi_cost = np.multiply(G.var, x).dot(xi)
            f_lam_der = 1.0 + (G.lambar * xi_cost) / \
                (2 * var_cost**(3.0 / 2.0))

        lam_prev = lam
        lam = lam - f / f_lam_der

        if lp:
            lam = np.maximum(lam, np.zeros(G.m))

        if iters > 5:
            break

    infeas = np.array(infeas)

    elapsed = time.time() - start
    return obj, elapsed, x, iter_objs, iter_elapsed, iter_xi_times, iter_var_times, mean, var,  infeas.mean(), lams, fs


def plot_flam(G, answer):
    import matplotlib.pyplot as plt
    lambar_l = [1e1, 1e3, 1e5, 1e7]
    for lambar in lambar_l:
        sigma_cost = get_sigma_cost(G, 0)
        lower_lam = lambar / (2 * sigma_cost)
        sigma_cost = get_sigma_cost(G, lower_lam)
        print('f_low: ', lower_lam - lambar / (2 * sigma_cost))
        sigma_cost = get_sigma_cost(G, 1e9)
        upper_lam = lambar / (2 * sigma_cost)
        sigma_cost = get_sigma_cost(G, upper_lam)

        print('f_high: ', upper_lam - lambar / (2 * sigma_cost))

        print('lambar: ', lambar)
        print('lower: ', lower_lam)
        print('upper: ', upper_lam)
        print('answer: ', answer)

        lamlist = []
        flamlist = []
        for lam in np.linspace(lower_lam, upper_lam, num=100, endpoint=False):
            sigma_cost = get_sigma_cost(G, lam)
            flam = lam - lambar / (2 * sigma_cost)
            lamlist.append(lam)
            flamlist.append(flam)

        print(lamlist[0], flamlist[0])
        plt.plot(lamlist, flamlist, 'o')
        plt.show()
        pdb.set_trace()


def alg3(**kwargs):
    print('-----ALG3------')

    if kwargs['G'] != None:
        G = copy.deepcopy(kwargs['G'])
    else:
        G = create_random_graph(**kwargs)
    start = time.time()
    if kwargs['R'] != None:
        R = copy.deepcopy(kwargs['R'])
    else:
        R = G.build_res_network()

    lambar = kwargs['lambar']
    stop_tol = kwargs['stop_tol']
    subproblem_tol = kwargs['subproblem_tol']
    precision_start_tol = kwargs['precision_start_tol']
    muscalar = kwargs['mu_scalar']
    disc_me_tot = 0

    disc_tot = 0
    nmcc_time_tot = 0
    scc_time_tot = 0
    lams = []
    # m=G.nxg.number_of_edges()
    # decoyG = copy.deepcopy(G)
    # decoyG.set_lambda(lam=0.01)
    # decoyR = copy.deepcopy(R)
    # nx.set_edge_attributes(decoyG.nxg, 0, 'mu')
    # varcost = None
    # discount, nmcc_time, scc_time, varcost=solve_mcf_sa(
    #     decoyG, decoyR, varcost, fully=False, tol=subproblem_tol)
    # disc_tot += discount
    # nmcc_time_tot += nmcc_time
    # scc_time_tot += scc_time
    # lam_high = lambar/(2*np.sqrt(varcost))
    # lam = lam_high/2
    varcost = kwargs['start_varcost']
    nmcc_time_start = kwargs['add_time']
    lam_high = kwargs['lam_high']
    lam = kwargs['startlam']
    # lams.append(lam)

    # lam_high, f_lam, lam_bound_time,varcost = get_upper_bound(G, R, lambar, subproblem_tol,muscalar)
    # uppertime = time.time()-start_uptime
    # lam_bound_time = 0
    # uppertime = 0
    # lam = lam_high

    lam_low = kwargs['lam_low']
    iters = 0
    found = False
    # gap=[]
    # gap_perc=[]
    subproblem_times = []
    times = []
    objs = []
    sub_start = time.time()
    discount = time.time()
    G.set_muscalar(kwargs['mu_scalar'])
    # sigmacost = np.sqrt(varcost)
    # algo_obj = (lambar * sigmacost + G.mu_cost() *
    #             kwargs['mu_scalar']) * kwargs['obj_scalar']
    # algo_obj= (lambar * sigmacost + G.mu_cost()) * kwargs['obj_scalar']
    # gap_perc.append(abs(1 - algo_obj / kwargs['cvx_obj']))
    # gap.append(abs(algo_obj - kwargs['cvx_obj']))
    # objs.append(algo_obj)
    # disc_tot += time.time() - discount
    subproblem_times.append(nmcc_time_start)
    times.append(time.time() - start - disc_tot + nmcc_time_start)
    # f = lam - lambar / (2 * sigmacost)
    f = 100
    lam_prev = None
    f_prev = np.inf
    xicost = 1e6
    nfultol = 1e-2
    nfultol2 = 1e-1

    while not found:

        iters += 1
        G.set_lambda(lam=lam)

        # if abs(f_lam) <= precision_start_tol:
        discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
            G, R, varcost, fully=True, nfullytol=nfultol, var_cost_ful=kwargs['varcostful'], difftol=1e-2)
        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time
        if f <= kwargs['precision_start_tol']:
            nfultol = 1e-4
            nfultol2 = 1e-1

        # else:
        #     if iters!=1:
        #         discount, nmcc_time, scc_time, varcost=solve_mcf_sa(
        #             G, R, varcost, fully=False, tol=subproblem_tol)
        #         disc_tot += discount
        #         nmcc_time_tot += nmcc_time
        #         scc_time_tot += scc_time
        #     else:
        #         pass

        subproblem_times.append(time.time() - sub_start - discount + nmcc_time)
        sigmacost = np.sqrt(varcost)
        disc_me = time.time()

        algo_obj = (lambar * sigmacost + G.mu_cost())

        # gap_perc.append(abs(1 - algo_obj / kwargs['cvx_obj']))
        # gap.append(abs(algo_obj - kwargs['cvx_obj']))
        objs.append(algo_obj)
        lams.append(lam)

        # print('algo_obj and f and lam ', algo_obj, f, lam)
        times.append(time.time() - start - disc_tot +
                     nmcc_time_tot - scc_time_tot + nmcc_time_start)

        disc_tot += time.time() - disc_me

        f_prev = f
        f = lam - lambar / (2 * sigmacost)

        if abs(f) < kwargs['stop_tol']:  # or algo_obj<=kwargs['cvx_obj']:
            found = True
            break

        # if np.linalg.norm(f_lam) <= np.linalg.norm(f_prev):

        if abs(f) <= abs(f_prev):
            if iters == 1:
                G.find_feasible_flow_xi()
                R_xi = G.build_res_network_xi()

            discount, nmcc_time, scc_time, xicost = solve_xi_mmc(
                G, R_xi, xicost, fully=True, nfullytol=nfultol2, tol=kwargs['sensitivity_tol'], difftol=1e0, vartop=kwargs['var_top'])

            disc_tot += discount
            nmcc_time_tot += nmcc_time
            scc_time_tot += scc_time
            f_lam_der = 1 + (lambar / 2) * (varcost**(-3 / 2)) * xicost
            lam_prev = lam
            lam = lam_prev - f / f_lam_der

            lam_high_prev = lam_high
            lam_low_prev = lam_low
            if lam_low <= lam <= lam_high:
                if f > 0:
                    lam_high = lam
                else:
                    lam_low = lam

            else:
                if lam_low > lam:
                    lam_high = lam_prev
                    lam = (lam_high + lam_low) / 2
                else:
                    lam_low = lam_prev
                    lam = (lam_high + lam_low) / 2

        else:
            # pdb.set_trace()
            lam = (lam_high_prev + lam_low_prev) / 2

    end = time.time()
    return algo_obj, iters, np.array(subproblem_times), times, objs, lams

#### EXPERIMENTS ####


def get_networks(experiment, tails, exponents, types, test=False):

    networks = []
    netgen_base = 'netgen'
    goto_base = 'goto'

    extension = 'pickle'

    for tail in tails:
        for atype in types:
            for exponent in exponents:
                net_name = netgen_base + atype + exponent + tail
                cur_run = experiment + '/' + netgen_base + \
                    '/' + net_name[:net_name.find('.')]
                cur_run = os.path.join(
                    EXPERIMENT_PATH, cur_run + "." + extension)
                if not test:
                    if not os.path.isfile(cur_run):
                        networks.append(net_name)
                else:
                    networks.append(net_name)
        if experiment == 'varying_lambar' or experiment == 'base_vs_reliable':
            continue
        else:
            for atype in types:
                if atype == '_lo_8_' or atype == '_lo_sr_':
                    continue
                for exponent in exponents:
                    net_name = goto_base + atype + exponent + tail
                    cur_run = experiment + '/' + goto_base + \
                        '/' + net_name[:net_name.find('.')]
                    cur_run = os.path.join(
                        EXPERIMENT_PATH, cur_run + "." + extension)
                    if not test:
                        if not os.path.isfile(cur_run):
                            networks.append(net_name)
                    else:
                        networks.append(net_name)

    return networks


def small_test_case():
    # testcase
    if test:
        lambar = 10
        G = MCF_DiGraph(lambar)
        G.nxg = nx.DiGraph()

        G.m = 3
        G.n = 3
        G.mu = np.array([0.0, 0.0, 0.0])
        G.var = np.array([0.5, 0.5, 1.0])
        G.b = np.array([1.0, 0, -1.0])
        G.cap = np.array([1.0, 1.0, 1.0])
        G.rows = [0, 0, 1, 1, 2, 2, ]
        G.cols = [0, 2, 0, 1, 1, 2]
        G.values = np.array([1.0, 1.0, -1.0, 1.0, -1.0, -1.0])

        solver_obj, solver_elapsed, solver_soln = cvxpy_solve(G)
        print('Solver finished in {} seconds, with objective {}'.format(
            solver_elapsed, solver_obj))

        pdb.set_trace()

        bs_obj, bs_elapsed, bs_soln, bs_iters = bs_cvxpy(
            G, lp=False, low=0, high=lambar, test=test)

        pdb.set_trace()


def graph_family_experiment(networks, lambar, record=True):
    experiment = 'graph_families'

    for network in networks:
        G = MCF_DiGraph(lambar)
        G.nxg = nx.DiGraph()
        print(network)

        if network.find('goto') >= 0:
            generator = 'goto'
            filename = 'networks/goto/' + network

        else:
            generator = 'netgen'
            filename = 'networks/netgen/' + network

        G = load_data(networkFileName=filename, G=G, generator=generator)

        run_name = experiment + '/' + generator + \
            '/' + network[:network.find('.')]
        cplex_saved = experiment + '/' + generator + '/' + 'cplex' + \
            '/' + network[:network.find('.')] + "_cplex_results"

        if not os.path.isfile(os.path.join(EXPERIMENT_PATH, cplex_saved) + '.pickle'):
            # print('mosek in progress..')
            # import cvxpy as cp
            # solver_obj, solver_elapsed, soln = cvxpy_any(G, solver='MOSEK')

            solver_obj, solver_elapsed, solver_soln, solver_infeas = cvxpy_solve(G)
        else:
            filename = cplex_saved
            data_dic = load_run(filename)
            solver_elapsed = data_dic['solver_elapsed']
            solver_obj = data_dic['solver_obj']
            solver_infeas = data_dic['solver_infeas']

        print('Solver finished in {} seconds, with objective {}'.format(
            solver_elapsed, solver_obj))

        # getting bounds:
        _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(G, lam=0, lp=True)
        var_cost = np.multiply(G.var, soln).dot(soln)
        lam_low = G.lambar / (2.0 * np.sqrt(var_cost))

        lam_high = G.lambar / 2.0

        nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs = nr_cvxpy(
            G, lp=False, low=lam_low, high=lam_high, lam_init=lam_low)

        nr_elapsed += lb_elapsed

        print('NR finished in {} seconds, with objective {}'.format(
            nr_elapsed, nr_obj))

        _, ub_elapsed, soln, _, _ = cvxpy_solve_additive(
            G, lam=0, lp=False, bound_lam=True)
        var_cost = np.multiply(G.var, soln).dot(soln)
        lam_high = G.lambar / (2.0 * np.sqrt(var_cost))

        bs_obj, bs_elapsed, _, bs_iter_objs, bs_iter_elapsed, bs_avg_infeas, bs_lams, bs_fs = bs_cvxpy(
            G, lp=False, low=lam_low, high=lam_high)

        bs_elapsed += lb_elapsed + ub_elapsed
        print('BSC finished in {} seconds, with objective {}'.format(
            bs_elapsed, bs_obj))

        t = PrettyTable(['Method', 'Soln Time',
                         '# Iters', 'Obj', 'Rel_Gap', 'Infeas'])
        t.add_row(['CPLEX', round(solver_elapsed, 2), '-',
                   solver_obj, '-', round(solver_infeas, 3)])
        t.add_row(['NR', round(nr_elapsed, 2), len(nr_iter_elapsed),
                   nr_obj, abs(solver_obj - nr_obj) / nr_obj, round(nr_avg_infeas, 3)])
        t.add_row(['BSC', round(bs_elapsed, 2), len(bs_iter_elapsed), bs_obj, abs(
            solver_obj - bs_obj) / bs_obj, round(bs_avg_infeas, 3)])
        print(t)
        pdb.set_trace()
        if record:

            solver_iter_objs, solver_iter_elapsed, presolve_time, ordering_time = parse_cplex_log()

            keys = ['network_name', 'G', 'solver_elapsed', 'solver_obj', 'solver_infeas',
                    'nr_elapsed', 'nr_iter_elapsed', 'nr_iter_objs', 'nr_obj',
                    'nr_xi_times', 'nr_var_times', 'mean', 'var', 'nr_avg_infeas', 'nr_lams', 'nr_fs', 'lambar', 'm', 'n',
                    'elapsed_lower_bound', 'elapsed_upper_bound',
                    'bs_obj', 'bs_iter_elapsed', 'bs_iter_objs', 'bs_elapsed', 'bs_avg_infeas', 'bs_lams', 'bs_fs', 'solver_iter_objs', 'solver_iter_elapsed', 'presolve_time', 'ordering_time']

            values = [network[:network.find('.')], G, solver_elapsed, solver_obj, solver_infeas,
                      nr_elapsed, nr_iter_elapsed, nr_iter_objs, nr_obj, nr_xi_times,
                      nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs, G.lambar, G.m, G.n, lb_elapsed, ub_elapsed,
                      bs_obj, bs_iter_elapsed, bs_iter_objs, bs_elapsed, bs_avg_infeas, bs_lams, bs_fs, solver_iter_objs, solver_iter_elapsed, presolve_time, ordering_time]

            run_dict = dict(zip(keys, values))

            if not os.path.isfile(cplex_saved):
                os.makedirs(cplex_saved, exist_ok=True)
                save_run(cplex_saved, run_dict)

            save_run(run_name, run_dict)


def varying_lambda_experiment(networks, lambars, record=True):
    experiment = 'varying_lambar'

    for network in networks:
        for lambar in lambars:
            G = MCF_DiGraph(lambar)
            G.nxg = nx.DiGraph()
            print(network)

            if network.find('goto') >= 0:
                generator = 'goto'
                filename = 'networks/goto/' + network

            else:
                generator = 'netgen'
                filename = 'networks/netgen/' + network

            G = load_data(networkFileName=filename, G=G, generator=generator)

            lam_dir = str(lambar).replace('.','')
            run_name = experiment + '/' + generator + \
                '/' + network[:network.find('.')] + '_' + lam_dir


            solver_weird = False
            solver_obj, solver_elapsed, solver_soln, solver_infeas = cvxpy_solve(
                G)
            print('Solver finished in {} seconds, with objective {}'.format(
                solver_elapsed, solver_obj))

            if solver_infeas > 2:
                solver_weird = True


            # getting bounds:
            _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(G, lam=0, lp=True)
            var_cost = np.multiply(G.var, soln).dot(soln)
            lam_low = G.lambar / (2.0 * np.sqrt(var_cost))

            lam_high = G.lambar / 2.0

            nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs = nr_cvxpy(
                G, lp=False, low=lam_low, high=lam_high, lam_init=lam_low)

            nr_elapsed += lb_elapsed

            print('NR finished in {} seconds, with objective {}'.format(
                nr_elapsed, nr_obj))

            _, ub_elapsed, soln, _, _ = cvxpy_solve_additive(
                G, lam=0, lp=False, bound_lam=True)
            var_cost = np.multiply(G.var, soln).dot(soln)
            lam_high = G.lambar / (2.0 * np.sqrt(var_cost))

            bs_obj, bs_elapsed, _, bs_iter_objs, bs_iter_elapsed, bs_avg_infeas, bs_lams, bs_fs = bs_cvxpy(
                G, lp=False, low=lam_low, high=lam_high)

            bs_elapsed += lb_elapsed + ub_elapsed
            print('BSC finished in {} seconds, with objective {}'.format(
                bs_elapsed, bs_obj))

            t = PrettyTable(['Method', 'Soln Time',
                             '# Iters', 'Obj', 'Rel_Gap', 'Infeas'])
            t.add_row(['CPLEX', round(solver_elapsed, 2), '-',
                       solver_obj, '-', round(solver_infeas, 3)])
            t.add_row(['NR', round(nr_elapsed, 2), len(nr_iter_elapsed),
                       nr_obj, abs(solver_obj - nr_obj) / nr_obj, round(nr_avg_infeas, 3)])
            t.add_row(['BSC', round(bs_elapsed, 2), len(bs_iter_elapsed), bs_obj, abs(
                solver_obj - bs_obj) / bs_obj, round(bs_avg_infeas, 3)])
            print(t)

            if record:

                solver_iter_objs, solver_iter_elapsed, presolve_time, ordering_time = parse_cplex_log()

                keys = ['network_name', 'G', 'solver_elapsed', 'solver_obj', 'solver_infeas',
                        'nr_elapsed', 'nr_iter_elapsed', 'nr_iter_objs', 'nr_obj',
                        'nr_xi_times', 'nr_var_times', 'mean', 'var', 'nr_avg_infeas', 'nr_lams', 'nr_fs', 'lambar', 'm', 'n',
                        'elapsed_lower_bound', 'elapsed_upper_bound',
                        'bs_obj', 'bs_iter_elapsed', 'bs_iter_objs', 'bs_elapsed', 'bs_avg_infeas', 'bs_lams', 'bs_fs', 'solver_iter_objs', 'solver_iter_elapsed', 'presolve_time', 'ordering_time']

                values = [network[:network.find('.')], G, solver_elapsed, solver_obj, solver_infeas,
                          nr_elapsed, nr_iter_elapsed, nr_iter_objs, nr_obj, nr_xi_times,
                          nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs, G.lambar, G.m, G.n, lb_elapsed, ub_elapsed,
                          bs_obj, bs_iter_elapsed, bs_iter_objs, bs_elapsed, bs_avg_infeas, bs_lams, bs_fs, solver_iter_objs, solver_iter_elapsed, presolve_time, ordering_time]

                run_dict = dict(zip(keys, values))

                save_run(run_name, run_dict)


def base_vs_reliable_expr(networks, lambars):

    experiment = 'base_vs_reliable'

    for network in networks:
        for lambar in lambars:
            G = MCF_DiGraph(lambar)
            G.nxg = nx.DiGraph()
            print(network)

            if network.find('goto') >= 0:
                generator = 'goto'
                filename = 'networks/goto/' + network

            else:
                generator = 'netgen'
                filename = 'networks/netgen/' + network

            G = load_data(networkFileName=filename, G=G, generator=generator)
            
            lam_dir = str(lambar).replace('.','')

            run_name = experiment + '/' + generator + \
                '/' + network[:network.find('.')] + '_' + lam_dir
            
            solver_weird = False
            reliable_obj, reliable_elapsed, _, solver_infeas = cvxpy_solve(G)

            if solver_infeas > 2:
                solver_weird = True


            print('reliable_obj is {}'.format(reliable_obj))
            orig_var = copy.deepcopy(G.var)

            G.var = np.zeros(G.m)

            _, base_elapsed, base_x, _ = cvxpy_solve(G)
            base_x = base_x[:-1]
            var_cost = np.multiply(orig_var, base_x).dot(base_x)
            base_obj = G.mu.dot(base_x) + G.lambar * (np.sqrt(var_cost))

            print('base_obj is {}'.format(base_obj))

            t = PrettyTable(['Approach', 'Objective', 'Rel-Gap'])
            t.add_row(['Reliable', reliable_obj, '-'])
            t.add_row(['Base', base_obj, (base_obj - reliable_obj)/reliable_obj])
            print(t)

            keys = ['network_name', 'G', 'base_elapsed', 'base_obj',
                    'solver_obj', 'solver_elapsed', 'lambar', 'm', 'n', 'solver_weird']

            values = [network[:network.find('.')], G, base_elapsed, base_obj,
                      reliable_obj, reliable_elapsed, G.lambar, G.m, G.n, solver_weird]

            run_dict = dict(zip(keys, values))

            save_run(run_name, run_dict)

test = False
if test:
    small_test_case()


if test:
    tails = ['a.min']
    exponents = ['10']
    types = ['_8_']
    networks = get_networks(experiment, tails, exponents, types, test=True)
    lambar = 10
    graph_family_experiment(networks, lambar, record=False)


lambar = 10
tails = ['a.min']
exponents = ['12']
types = ['_sr_']
experiment = 'graph_families'
networks = get_networks(experiment, tails, exponents, types)
print(networks)
graph_family_experiment(networks, lambar)

#### check the output log for above^


# print('starting lambar experiments')
# lambars = [0.001, 0.1, 10, 1000]
# tails = ['a.min', 'b.min', 'c.min', 'd.min', 'e.min']
# exponents = ['11']
# types = ['_lo_8_', '_lo_sr_', '_sr_', '_8_']
# experiment = 'varying_lambar'
# networks = get_networks(experiment, tails, exponents, types)
# varying_lambda_experiment(networks, lambars)


# print('starting graph family experiments for exponent 14')
# tails = ['a.min', 'b.min', 'c.min', 'd.min', 'e.min']
# exponents = (np.arange(10, 15)).astype(str)
# types = ['_lo_8_', '_lo_sr_', '_sr_', '_8_']
# experiment = 'graph_families'
# networks = get_networks(experiment, tails, exponents, types)
# lambar = 10
# graph_family_experiment(networks, lambar)


# print('starting base vs reliable experiments')
# lambars = [0.1, 10, 1000]
# exponents = ['11']
# types = ['_sr_', '_8_']
# tails = ['a.min']
# experiment = 'base_vs_reliable'
# networks = get_networks(experiment, tails, exponents, types)
# base_vs_reliable_expr(networks, lambars)

