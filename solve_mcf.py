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

# import cvxpy as cp
import cplex
from prettytable import PrettyTable



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
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)
    prob.parameters.mip.display.set(0)

    prob.solve()

    obj = prob.solution.get_objective_value()
    soln = np.array(prob.solution.get_values())

    elapsed = time.time() - start

    return obj, elapsed, soln

def estimate_quad(G):


    start = time.time()

    prob = cplex.Cplex()


    underestimator = G.cap * G.var
    prob.variables.add(obj=underestimator, lb=np.zeros(G.m), ub=G.cap)

    prob.linear_constraints.add(rhs=G.b, senses='E' * G.n)
    prob.linear_constraints.set_coefficients(zip(G.rows, G.cols, G.values))

    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)
    prob.parameters.barrier.display.set(0)
    
    prob.solve()

    obj = prob.solution.get_objective_value()
    soln = np.array(prob.solution.get_values())

    # c.write("model.lp")

    elapsed = time.time() - start

    return obj, elapsed, soln, prob

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
                prob.parameters.barrier.convergetol.set(1e-1)
                prob.parameters.barrier.limits.iteration.set(14)
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

    return obj, elapsed, soln, prob


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
        prob.parameters.mip.display.set(0)
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

    return elapsed, soln, prob


def bs_cvxpy(G, low=0, high=1000, prob=None, lp=False, test=False):
    print('bsc')

    ###test
    test=test
    if test:
        warm_start = False
        for i in np.linspace(low, high/2, 300):
            lam = i
            obj_, elapsed, x, prob = cvxpy_solve_additive(G, lam, prob=prob, warm_start=warm_start, lp=lp, bound_lam=False)
            var_cost = np.multiply(G.var, x).dot(x)
            obj = G.mu.dot(x) + G.lambar * (np.sqrt(var_cost))
            
            f1 = lam - G.lambar / (2.0 * np.sqrt(var_cost))
            f2 = lam  - G.lambar / np.sqrt(var_cost)
            

            t = PrettyTable(['lam_obj', 'lambar_obj', 'time', 'derivative_eq', 'func_eq', 'lambda'])
            t.add_row([obj_, obj, time.time()-start, f1, f2, lam])
            print(t)

            if abs(f1) < 1e-1 or abs(f2) < 1e-1:
                print('now debug')
                print(f1, f2, obj_, obj)
                pdb.set_trace()

        pdb.set_trace()
    start = time.time()


    stop_tol = 1e-4

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

    while not found:
        iters += 1
        mid_prev = mid
        mid = (high + low) / 2.0

        if iters > 1 and warm_start == False:
            warm_start = True

        obj, elapsed, x, prob = cvxpy_solve_additive(G, mid, prob=prob, warm_start=warm_start, lp=lp)
        var_cost = np.multiply(G.var, x).dot(x)
        obj = G.mu.dot(x) + G.lambar * (np.sqrt(var_cost))
        print(obj, time.time() - start, f, mid)
        
        iter_objs.append(obj)
        iter_elapsed.append(time.time() - start)

        if lp:
            f = mid - G.lambar * np.multiply(G.var, x) / np.sqrt(var_cost)
        else:
            f = mid - G.lambar / (2.0 * np.sqrt(var_cost))

        if lp:
            if np.all(f) < stop_tol or iters > 9:
                found = True
                break
        else:
            if abs(mid_prev - mid) < stop_tol:
                found = True
                break

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

    elapsed = time.time() - start
    return obj, elapsed, x, iter_objs, iter_elapsed,


def nr_cvxpy(G, low=0, high=1000, prob=None, lp=False, lam_init=None):
    print('nr')
    start = time.time()
    stop_tol = 1e-4

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

    while not found:
        iters += 1

        if iters > 1 and warm_start == False:
            warm_start = True
            warm_start_xi = True

        obj, elapsed, x, prob = cvxpy_solve_additive(
            G, lam, prob=prob, warm_start=warm_start, lp=lp)
        iter_var_times.append(elapsed)

        var_cost = np.multiply(G.var, x).dot(x)

        cur_mean = G.mu.dot(x)
        obj = cur_mean + G.lambar * (np.sqrt(var_cost))
        iter_objs.append(obj)
        iter_elapsed.append(time.time() - start)
        
        mean.append(cur_mean)
        var.append(var_cost)
        print(obj, time.time() - start, f, lam)

        if lp:
            f = lam - G.lambar * np.multiply(G.var, x) / np.sqrt(var_cost)
        else:
            f = lam - G.lambar / (2.0 * np.sqrt(var_cost))

        if lp:
            if np.all(f) < stop_tol or iters > 5:
                found = True
                break
        else:
            if abs(lam_prev - lam) < stop_tol:
                found = True
                break


        elapsed_xi, xi, prob_xi = cvxpy_solve_xi(
            G, x, lam, prob=prob_xi, warm_start=warm_start_xi, lp=lp, iters=iters)
        iter_xi_times.append(elapsed_xi)

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


    elapsed = time.time() - start
    return obj, elapsed, x, iter_objs, iter_elapsed, iter_xi_times, iter_var_times, mean, var


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

test=False
#testcase
if test:
    lambar = 10
    G = MCF_DiGraph(lambar)
    G.nxg = nx.DiGraph()

    G.m = 3
    G.n = 3
    G.mu = np.array([0.0, 0.0, 0.0])
    G.var = np.array([0.5, 0.5, 1.0])
    G.b = np.array([1.0, 0 , -1.0])
    G.cap = np.array([1.0, 1.0, 1.0])
    G.rows = [0,0,1,1,2,2,]
    G.cols = [0,2,0,1,1,2]
    G.values = np.array([1.0,1.0,-1.0,1.0,-1.0,-1.0])

    solver_obj, solver_elapsed, solver_soln = cvxpy_solve(G)
    print('Solver finished in {} seconds, with objective {}'.format(solver_elapsed, solver_obj))

    pdb.set_trace()

    bs_obj, bs_elapsed, bs_soln, bs_iters = bs_cvxpy(
                G, lp=False, low=0, high=lambar, test=test)

    pdb.set_trace()

### testcase end


netgen_base = 'netgen'

tails = ['a.min', 'b.min', 'c.min', 'd.min']

exponents = (np.arange(10, 15)).astype(str)
# exponents = (np.arange(10, 16)).astype(str)

types = ['_lo_8_', '_lo_sr_', '_sr_', '_8_']

networks = []

experiment = 'graph_families'
extension = 'pickle'


goto_base = 'goto'

for tail in tails:
    for atype in types:
        for exponent in exponents:
            net_name = netgen_base + atype + exponent + tail
            cur_run = experiment + '/' + netgen_base + '/' + net_name[:net_name.find('.')] 
            cur_run = os.path.join(EXPERIMENT_PATH, cur_run + "." + extension)
            if not os.path.isfile(cur_run):
                networks.append(net_name)

    for atype in types[2:]:
        for exponent in exponents:
            net_name = goto_base + atype + exponent + tail
            cur_run = experiment + '/' + goto_base + '/' + net_name[:net_name.find('.')] 
            cur_run = os.path.join(EXPERIMENT_PATH, cur_run + "." + extension)
            if not os.path.isfile(cur_run):
                networks.append(net_name)


def graph_family_experiment(networks, lambar):
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

        solver_obj, solver_elapsed, _ = cvxpy_solve(G)
        print('Solver finished in {} seconds, with objective {}'.format(solver_elapsed, solver_obj))

        ## getting bounds:
        _, lb_elapsed, soln, prob = cvxpy_solve_additive(G, lam=0, lp=True)
        var_cost = np.multiply(G.var, soln).dot(soln)
        lam_low = G.lambar/(2.0*np.sqrt(var_cost))

        lam_high = G.lambar/2.0

        nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var = nr_cvxpy(
            G, lp=False, low=lam_low, high=lam_high, lam_init=lam_low)
        
        nr_elapsed += lb_elapsed

        print('NR finished in {} seconds, with objective {}'.format(nr_elapsed, nr_obj))

        _, ub_elapsed, soln, prob = cvxpy_solve_additive(G, lam=0, lp=False, bound_lam=True)
        var_cost = np.multiply(G.var, soln).dot(soln)
        lam_high = G.lambar/(2.0*np.sqrt(var_cost))

        bs_obj, bs_elapsed, _, bs_iter_objs, bs_iter_elapsed = bs_cvxpy(
            G, lp=False, low=lam_low, high=lam_high)

        bs_elapsed += lb_elapsed + ub_elapsed
        print('BSC finished in {} seconds, with objective {}'.format(bs_elapsed, bs_obj))


        t = PrettyTable(['Method', 'Solution Time',
                         '# Iterations', 'Objective', 'Rel_Gap'])
        t.add_row(['CPLEX', solver_elapsed, 1, solver_obj, 0])
        t.add_row(['NR', nr_elapsed, len(nr_iter_elapsed),
                   nr_obj, abs(solver_obj - nr_obj) / nr_obj])
        t.add_row(['BSC', bs_elapsed, len(bs_iter_elapsed), bs_obj, abs(solver_obj - bs_obj)/bs_obj])
        print(t)


        run_name = experiment + '/' + generator + '/' + network[:network.find('.')]

        keys = ['network_name', 'G', 'solver_elapsed', 'solver_obj',
            'nr_elapsed', 'nr_iter_elapsed', 'nr_iter_objs', 'nr_obj',
            'nr_xi_times', 'nr_var_times', 'mean', 'var', 'lambar', 'm', 'n',
            'elapsed_lower_bound', 'elapsed_upper_bound',
            'bs_obj', 'bs_iter_elapsed', 'bs_iter_objs', 'bs_elapsed']

        values = [network[:network.find('.')], G, solver_elapsed, solver_obj,
            nr_elapsed, nr_iter_elapsed, nr_iter_objs, nr_obj, nr_xi_times,
            nr_var_times, mean, var, G.lambar, G.m, G.n, lb_elapsed, ub_elapsed,
            bs_obj, bs_iter_elapsed, bs_iter_objs, bs_elapsed]

        run_dict = dict(zip(keys, values)) 

        save_run(run_name , run_dict)


lambar = 10

graph_family_experiment(networks, lambar)
pdb.set_trace()










# experiment = 'varying_lambar'

# lambars = [0.001, 0.01, 0.1, 1, 10, 100]

# networks = ['netgen_8_16a.min', 'netgen_lo_8_16a.min', 'netgen_sr_16a.min']

# for network in networks:
#     for lambar in lambars:
#         G = MCF_DiGraph(lambar)
#         G.nxg = nx.DiGraph()
#         print(network)

#         solver_obj, solver_elapsed, _ = cvxpy_solve(G)
#         print('Solver finished in {} seconds, with objective {}'.format(solver_elapsed, solver_obj))

#         nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var = nr_cvxpy(
#             G, lp=False, high=G.lambar / 2.0)
#         print('NR finished in {} seconds, with objective {}'.format(nr_elapsed, nr_obj))

#         from prettytable import PrettyTable
#         t = PrettyTable(['Method', 'Solution Time',
#                          '# Iterations', 'Objective', 'Rel_Gap'])
#         t.add_row(['CPLEX', solver_elapsed, 1, solver_obj, 0])
#         t.add_row(['NewtonRaphson', nr_elapsed, len(nr_iter_elapsed),
#                    nr_obj, abs(solver_obj - nr_obj) / solver_obj])
#         print(t)

#         run_name = experiment + '/' + generator + '/' + network[:network.find('.')]
#         run_dict = {}
#         run_dict['network_name'] = network[:network.find('.')]
#         run_dict['G'] = G
#         run_dict['solver_elapsed'] = solver_elapsed
#         run_dict['solver_obj'] = solver_obj
#         run_dict['nr_elapsed'] = nr_elapsed
#         run_dict['nr_iter_elapsed'] = nr_iter_elapsed
#         run_dict['nr_iter_objs'] = nr_iter_objs
#         run_dict['nr_obj'] = nr_obj
#         run_dict['nr_obj'] = nr_xi_times
#         run_dict['nr_obj'] = nr_var_times
#         run_dict['mean'] = mean
#         run_dict['var'] = var


#         save_run(run_name , run_dict)

# pdb.set_trace()


# experiment = 'base_vs_reliable'

# lambar = 0.7

# for network in networks:

#     G = MCF_DiGraph(lambar)
#     G.nxg = nx.DiGraph()
#     print(network)

#     reliable_obj, _, _ = cvxpy_solve(G)
#     print('reliable_obj is {}'.format(reliable_obj))
#     orig_var = copy.deepcopy(G.var)

#     G.var = np.zeros(G.m)
    
#     _, _, base_x = cvxpy_solve(G)
#     var_cost = np.multiply(orig_var, base_x).dot(base_x)
#     base_obj = G.mu.dot(base_x) + G.lambar * (np.sqrt(var_cost))

#     print('reliable_obj is {}'.format(base_obj))

#     from prettytable import PrettyTable
#     t = PrettyTable(['Approach', 'Objective'])
#     t.add_row(['Reliable', reliable_obj])
#     t.add_row(['Base', base_obj])
#     print(t)

#     run_name = experiment + '/' + generator + '/' + network[:network.find('.')]
#     run_dict = {}
#     run_dict['network_name'] = network[:network.find('.')]
#     run_dict['G'] = G
#     run_dict['base_obj'] = base_obj
#     run_dict['reliable_obj'] = reliable_obj

# pdb.set_trace()



