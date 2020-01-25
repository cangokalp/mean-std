from utils import *
from MCF_DiGraph import *

import matplotlib.pyplot as plt
import pandas
import random
import itertools
import copy
from tabulate import tabulate
import cProfile
import pstats
# sys.path.append('/usr/local/Cellar/graph-tool/2.27_1/lib/python3.7/site-packages/')
# import graph_tool as gt

import cvxpy as cp
# from cvxopt import matrix, solvers, spmatrix
import cplex
# import docplex
# from gurobipy import *


def cvxpy_solve(G):
    start = time.time()

    lambar = G.lambar
    x = cp.Variable(G.m)

    constraints = [0 <= x, x <= G.cap, G.A@x == G.b]

    objective = cp.Minimize(G.mu * x + lambar * cp.norm(cp.multiply(G.sigma, x), 2))
    prob = cp.Problem(objective, constraints)
    
    # mosek_params = {
    #                 "MSK_IPAR_INTPNT_MAX_ITERATIONS": 1000,
    #                 "MSK_DPAR_INTPNT_CO_TOL_REL_GAP":1e-12
    #                 }
    mosek_params = {}
    result = prob.solve(solver='MOSEK', verbose=False)

    # prob.unpack_results(cvx.ECOS, solver_output)
    # result = prob.solve(solver='MOSEK', verbose=True, mosek_params={,'MSK_DPAR_INTPNT_CO_TOL_INFEAS':1e-30})
    
    soln = x.value
    obj = objective.value

    elapsed = time.time() - start
   
    return obj, elapsed, soln


def cvxpy_solve_lin(G, lam, prob=None, warm_start=False):

    start = time.time()

    if not warm_start:
        weight = cp.Parameter(nonneg=True)
        x = cp.Variable(G.m)

        constraints = [0 <= x, x <= G.cap, G.A@x == G.b]
        # objective = cp.Minimize(G.mu*x + weight*np.ones(G.m)*x)

        P = np.multiply(G.var,np.eye(G.m))
        objective = cp.Minimize(G.mu*x + weight*cp.quad_form(x,P))

        prob = cp.Problem(objective, constraints)

    # weight.value = lam
    prob.parameters()[0].value = lam
    prob.solve(solver='MOSEK', verbose=False, warm_start=warm_start)

    # soln = x.value
    soln = [variable.value for variable in prob.variables()][0]
    obj = prob.value
    elapsed = time.time() - start
        
    return obj, elapsed, soln, prob


def cvxpy_solve_xi(G, soln, lam):

    start = time.time()

    x_ = soln
    xi = cp.Variable(G.m)

    constraints = [G.A@xi == 0]


    diff = x_ - G.cap

    x_zero = np.argwhere(abs(x_) < 1e-6).ravel()
    x_nonzero = np.argwhere(x_ > 1e-6).ravel()
    x_u = np.argwhere(abs(diff) < 1e-6).ravel()

    x_btw = list(set(x_nonzero).difference(set(x_u)))

    constraints.append(xi[x_zero] == 0)
    # constraints.append(xi[x_zero] <= G.cap[x_zero])

    constraints.append(xi[x_u] <= 0)
    constraints.append(xi[x_u] >= -x_[x_u])

    constraints.append(xi[x_btw] <= G.cap[x_btw] - x_[x_btw])
    constraints.append(xi[x_btw] >= - x_[x_btw])


    # objective = cp.Minimize(np.ones(G.m)*xi)
    P = np.multiply(G.var/100,np.eye(G.m))
    objective = cp.Minimize(lam*cp.quad_form(xi,P) + 2*np.multiply(G.var/100,x_)*xi)

    prob = cp.Problem(objective, constraints)

    result = prob.solve(solver='GUROBI', verbose=True)

    soln = xi.value
    obj = objective.value

    elapsed = time.time() - start
    return elapsed, soln 



def bs_cvxpy(G, low=0, high=500.0, prob=None):

    start = time.time()

    kwargs = {}
    kwargs['stop_tol'] = 1e-5

    f = 100
    found = False
    iters = 0

    # high = np.ones(G.m) * high
    # low = np.ones(G.m) * low

    if prob is not None:
        warm_start = True
    else:
        warm_start = False
    
    while not found:
        iters += 1
        mid = (high + low) / 2.0

        if iters > 1 and warm_start == False:
            warm_start = True

        obj, elapsed_1, x, prob = cvxpy_solve_lin(G, mid, prob=prob, warm_start=warm_start)
        real_obj = G.mu.dot(x) + G.lambar*(np.sqrt(x.dot(np.multiply(G.var,np.eye(G.m))).dot(x)))
        
        # elapsed_2, xi = cvxpy_solve_xi(G, x)
        var_cost = x.dot(np.multiply(G.var,np.eye(G.m))).dot(x)

        f = mid - G.lambar / (2.0 * np.sqrt(var_cost))        
        # f = mid - G.lambar*(np.multiply(G.var,x).dot(xi))/(np.sqrt(var_cost)*sum(xi))

        print(obj, real_obj, mid, f)

        if abs(f) < kwargs['stop_tol']:  
            found = True
            break

        if np.sign(f) == np.sign(1):
            high = mid
        else:
            low = mid


    elapsed = time.time() - start

    return real_obj, elapsed, x




def nr_cvxpy(G):
    start = time.time()
    kwargs = {}
    kwargs['stop_tol'] = 1e-9

    lam = 0.5
    f = 100
    found = False
    iters = 0
    warm_start = False
    prob = None

    while not found:
        iters += 1
        obj, elapsed, soln, prob = cvxpy_solve_lin(G, lam, prob=prob, warm_start=warm_start)
        
        var_cost = x.dot(np.multiply(G.var,np.eye(G.m))).dot(x)
        real_obj = G.mu.dot(x) + G.lambar*(np.sqrt(var_cost))

        f = lam - G.lambar / (2 * np.sqrt(var_cost))

        if abs(f) < kwargs['stop_tol'] or iters > 10:  
            found = True
            break

        elapsed, xi = cvxpy_solve_xi(G, soln, lam)

        xi_cost = np.multiply(G.var, x).dot(xi)

        f_lam_der = 1 + (G.lambar*xi_cost)/(2*var_cost**(-3.0/2.0))
        
        lam = lam - f/f_lam_der
        
        print(obj, real_obj, lam, f)
        pdb.set_trace()
    elapsed = time.time() - start
    print('nr_iter_num: ', iters)
    return obj, elapsed


def get_upper_bound(G, R, lambar, tol, muscalar, disc_tot=0, nmcc_tot=0):
    lam_high = lambar
    G.set_lambda(lam_high)
    G.set_muscalar(muscalar)
    start = time.time()
    varcost = 1e6
    discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
        G, R, varcost, fully=False, tol=tol)
    # solve_var(G, lam_high)
    # discount, nmcc_time, scc_time = solve_mcf_sa(G, R, fully=True)
    # discount, nmcc_time, scc_time = solve_mcf_sa(G, R)

    disc_tot += discount
    nmcc_tot += nmcc_time
    sigma_lam = np.sqrt(varcost)
    f_high = lam_high - float(lambar) / float((2 * sigma_lam))
    while f_high <= 0:
        lam_high = 2 * lam_high
        G.set_lambda(lam_high)
        discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
            G, R, varcost, fully=False, tol=tol)
        # solve_var(G, lam_high)

        # discount, nmcc_time, scc_time = solve_mcf_sa(G, R, fully=True)

        # discount, nmcc_time, scc_time = solve_mcf_sa(G, R)
        disc_tot += discount
        nmcc_tot += nmcc_time
        sigma_lam = np.sqrt(varcost)
        f_high = lam_high - float(lambar) / float((2 * sigma_lam))

    end = time.time()
    elapsed = end - start - disc_tot + nmcc_tot
    return lam_high, f_high, elapsed, varcost


def get_sigma_cost(G, lam):
    G.set_lambda(lam=lam)
    solve_var(G, lam)
    return np.sqrt(G.var_cost())


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


def NR_solve(**kwargs):
    print('-----NR------')

    if kwargs['G'] != None:
        G = copy.deepcopy(kwargs['G'])
    else:
        G = create_random_graph(**kwargs)

    start = time.time()
    if kwargs['R'] != None:
        R = copy.deepcopy(kwargs['R'])
    else:
        R = G.build_res_network()

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

    # lam_high, f_lam, lam_bound_time,varcost = get_upper_bound(G, R, lambar, subproblem_tol, muscalar)
    # uppertime = time.time()-start_uptime
    # lam = lam_high
    # uppertime = 0
    # lam_bound_time = 0

    iters = 0
    found = False
    subproblem_times = []
    # gap=[]
    # gap_perc=[]
    times = []
    objs = []
    sub_start = time.time()
    # discount = time.time()
    G.set_muscalar(kwargs['mu_scalar'])

    # sigmacost = np.sqrt(varcost)
    # algo_obj = (lambar * sigmacost + G.mu_cost() *
    #             kwargs['mu_scalar']) * kwargs['obj_scalar']
    # # algo_obj= (lambar * sigmacost + G.mu_cost()) * kwargs['obj_scalar']

    # gap_perc.append(abs(1 - algo_obj / kwargs['cvx_obj']))
    # gap.append(abs(algo_obj - kwargs['cvx_obj']))
    # objs.append(algo_obj)
    # disc_tot += time.time() - discount
    subproblem_times.append(nmcc_time_start)
    times.append(time.time() - start - disc_tot + nmcc_time_start)

    elapsed = 0
    cvx_elapsed = 0
    xicost = 1e6
    # f = lam - lambar / (2 * sigmacost)
    f = 100
    nfultol = 1e-2
    nfultol2 = 1e-1

    while not found:
        iters += 1
        G.set_lambda(lam=lam)
        sub_start = time.time()
        # if abs(f_lam) <= precision_start_tol:
        discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
            G, R, varcost, fully=True, nfullytol=nfultol, var_cost_ful=kwargs['varcostful'], vartop=G.vartop, difftol=1e-2)
        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time
        if f <= kwargs['precision_start_tol']:
            nfultol = 1e-4
            nfultol2 = 1e-1

        # else:
        #     if iters != 1:
        #         discount, nmcc_time, scc_time, varcost=solve_mcf_sa(
        #             G, R, varcost, fully=False, tol=subproblem_tol)
        #         disc_tot += discount
        #         nmcc_time_tot += nmcc_time
        #         scc_time_tot += scc_time
        #     else:
        #         pass

        subproblem_times.append(
            time.time() - sub_start - discount + nmcc_time - scc_time)
        disc_me = time.time()
        sigmacost = np.sqrt(varcost)

        # algo_obj = (lam * varcost + G.mu_cost() *
        #             kwargs['mu_scalar']) * kwargs['obj_scalar']
        algo_obj = (G.lambar * sigmacost + G.mu_cost())

        # algo_obj= (lambar * sigmacost + G.mu_cost()) * kwargs['obj_scalar']
        # print(algo_obj)

        # gap.append(abs(algo_obj - kwargs['cvx_obj']))
        # gap_perc.append(abs(1 - algo_obj / kwargs['cvx_obj']))
        objs.append(algo_obj)
        times.append(time.time() - start - elapsed -
                     disc_tot + nmcc_time_tot - scc_time_tot + nmcc_time_start)
        # print('algo_obj and f and lam', algo_obj, f, lam)
        disc_tot += time.time() - disc_me

        sigmacost = np.sqrt(varcost)
        f = lam - lambar / (2 * sigmacost)
        # print(f_lam)
        # kwargs['stop_tol']
        lams.append(lam)

        if abs(f) < kwargs['stop_tol']:  # or algo_obj<=kwargs['cvx_obj']:
            found = True
            break

        if iters == 1:
            G.find_feasible_flow_xi()
            R_xi = G.build_res_network_xi()

        # if abs(f_lam) < precision_start_tol:
        #     discount, nmcc_time, scc_time, xicost=solve_xi_mmc(
        #         G, R_xi, xicost, fully=True)
        # else:
        #     discount, nmcc_time, scc_time, xicost=solve_xi_mmc(
        #         G, R_xi, xicost, fully=False, tol=subproblem_tol)

        new_flow, new_cost = solve_xi_new(G, lam)

        old_flow, old_cost = solve_xi(G, lam)

        discount, nmcc_time, scc_time, xicost = solve_xi_mmc(
            G, R_xi, xicost, fully=True, nfullytol=1e-3, tol=kwargs['sensitivity_tol'], difftol=1e0, vartop=G.vartop)
        pdb.set_trace()

        discount, nmcc_time, scc_time, xicost = solve_xi_mmc(
            G, R_xi, xicost, fully=True, nfullytol=nfultol2, tol=kwargs['sensitivity_tol'], difftol=1e0, vartop=G.vartop)

        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time

        if kwargs['way'] == 'newxicost':
            discount, nmcc_time, scc_time, xicost = solve_xi_mmc(
                G, R_xi, xicost, fully=True, nfullytol=nfultol2, tol=kwargs['sensitivity_tol'], difftol=1e0, vartop=G.vartop)

            disc_tot += discount
            nmcc_time_tot += nmcc_time
            scc_time_tot += scc_time

        elif kwargs['way'] == 'oldformul':
            old_flow, old_cost = solve_xi(G, lam)
            xicost = old_cost

        elif kwargs['way'] == 'newformul':
            new_flow, new_cost = solve_xi_new(G, lam)
            xicost = new_cost

        f_lam_der = 1 + (lambar / 2) * (varcost**(-3 / 2)) * xicost
        # f_lam_der = 1
        # print(f, f_lam_der, f/f_lam_der)

        # lam = lam - f / f_lam_der
        lam = lam - f
        print('lam: ', lam)

    end = time.time()
    # print('lam: ', lam)
    print(end - start)
    print(lam)
    pdb.set_trace()
    # cur_r_dict = nx.to_dict_of_dicts(G.nxg)
    # cc = [v for v in cur_r_dict.keys()]
    # list_ofed = []
    # for c in cc:
    #     if 10000 in cur_r_dict[c].keys():
    #         list_ofed.append(c)

    # flows = []
    # for e in list_ofed:
    #     flows.append(G.nxg[e][10000]['flow'])

    # flows = sum(np.array(flows))
    # print(flows)
    return algo_obj, iters, np.array(subproblem_times), times, objs, lams


def bisect_solve(**kwargs):
    print('-----BISECTION------')
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
    subproblem_tol = kwargs['subproblem_tol_bs']
    precision_start_tol = kwargs['precision_start_tol_bs']
    muscalar = kwargs['mu_scalar']
    found = False
    iters = 0

    disc_tot = 0
    nmcc_time_tot = 0
    scc_time_tot = 0
    lams = []

    varcost = kwargs['start_varcost_bs']
    nmcc_time_start = kwargs['add_time_bs']
    lam_high = kwargs['lam_high_bs']
    lam = kwargs['startlam_bs']
    low = kwargs['lam_low']

    # sigmacost = np.sqrt(varcost)
    high = lam_high
    f = 100

    objs = []
    subproblem_times = []
    times = []
    G.set_muscalar(kwargs['mu_scalar'])

    # discount = time.time()
    # algo_obj = (lambar * sigmacost + G.mu_cost() *
    #             kwargs['mu_scalar']) * kwargs['obj_scalar']

    # objs.append(algo_obj)
    # disc_tot += time.time() - discount
    subproblem_times.append(nmcc_time_start)
    times.append(time.time() - start - disc_tot + nmcc_time_start)
    nfultol = 1e-2

    while not found:
        iters += 1
        mid = (high + low) / 2.0
        sub_start = time.time()
        G.set_lambda(lam=mid)

        if abs(f) <= kwargs['precision_start_tol_bs']:
            nfultol = 1e-4

        discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
            G, R, varcost, fully=True, nfullytol=nfultol, var_cost_ful=kwargs['varcostful'], vartop=kwargs['var_top'], mutop=kwargs['mu_top'], difftol=1e-2)
        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time

        subproblem_times.append(
            time.time() - sub_start - discount + nmcc_time - scc_time)

        sigmacost = np.sqrt(varcost)

        disc_me = time.time()
        algo_obj = (lambar * sigmacost + G.mu_cost())

        # algo_obj = (mid * varcost + G.mu_cost() *
        #             kwargs['mu_scalar']) * kwargs['obj_scalar']
        times.append(time.time() - start - disc_tot +
                     nmcc_time_tot - scc_time_tot + nmcc_time_start)

        objs.append(algo_obj)
        lams.append(mid)

        print('algo_obj: {}, f: {}, lam: {}'.format(algo_obj, f, mid))

        disc_tot += time.time() - disc_me
        f = mid - float(lambar) / float((2 * sigmacost))

        if abs(f) < kwargs['stop_tol']:
            found = True
            break

        if np.sign(f) == np.sign(1):
            high = mid
        else:
            low = mid

        mid_prev = mid

    end = time.time()

    # cur_r_dict = nx.to_dict_of_dicts(G.nxg)
    # cc = [v for v in cur_r_dict.keys()]
    # list_ofed = []
    # endnode = len(G.nxg.nodes())
    # for c in cc:
    #     if endnode in cur_r_dict[c].keys():
    #         list_ofed.append(c)

    # flows = []
    # for e in list_ofed:
    #     flows.append(G.nxg[e][endnode]['flow'])

    # sflows = sum(np.array(flows))
    # print('totflow: ', sflows)
    # x = nx.get_edge_attributes(G.nxg, 'flow')

    return algo_obj, iters, np.array(subproblem_times), times, objs, lams, 5


def run_experiment(num_nodes, lams, mus, variances, d_tops, experiment_name, multi_arcs=False, seeds=9 * np.arange(10), repeat=10, test4dense=False, fixmsk=False, fix_bs=False, start_range=0, test=False):
    num_nodes = num_nodes
    lams = lams
    mus = mus
    variances = variances
    parameterszip = list(itertools.product(
        num_nodes, lams, mus, variances, d_tops))
    experiment_name = experiment_name + '/'
    lamBar = lams[0]
    for parameter_set in parameterszip:
        n = parameter_set[0]
        lamBar = parameter_set[1]
        mu_top = parameter_set[2]
        var_top = parameter_set[3]
        d_top = parameter_set[4]
        cur_lam = parameter_set[1]
        if multi_arcs:
            num_arcs = [n * 2, n * 8]
        else:
            if test4dense:
                num_arcs = [n * 8]
            else:
                num_arcs = [n * 4]

        for narcs in num_arcs:

            params = {}
            # if lamBar / 5 == 1 or lamBar / 1 == 1:
            #     scale = 1e1
            # elif lamBar / 10 == 1:
            #     scale = 1e2
            # elif lamBar / 100 == 1:
            #     scale = 1e3
            # else:
            #     scale = 1.0
            scale = 1.0
            params['mu_top'] = mu_top
            params['var_top'] = var_top
            params['d_top'] = d_top
            params['scale'] = 1
            params['lambar'] = lamBar / scale
            params['num_nodes'] = n
            params['sensitivity_tol'] = 1e-3
            params['stop_tol'] = 1e-5
            params['precision_start_tol'] = 1e-2
            params['precision_start_tol_bs'] = 1e-2
            params['subproblem_tol'] = 1e-1
            params['subproblem_tol_bs'] = 1e-1
            params['cvx_obj'] = None
            params['mu_scalar'] = 1.0 / scale
            params['var_scalar'] = 1.0  # scale*scale
            params['obj_scalar'] = scale
            params['num_arcs'] = narcs
            params['varcostnoful'] = 1e0
            params['varcostful'] = 1e0
            params['xicostnoful'] = 1e0
            seeds = seeds

            times_cvx_list = []
            times_bs_list = []
            times_nr_list = []
            gap_cvx_list = []
            gap_bs_list = []
            gap_nr_list = []

            orig_var_top = var_top
            print(repeat)
            for i in range(start_range, repeat):
                print(i)
                first = True

                gap_cvx = []
                gap_perc_cvx = []
                times_cvx = []
                seed = seeds[i]
                params['seed'] = seed
                lamstr = str(lamBar / scale)
                lamstr = lamstr.replace(".", "-")

                save_extension = lamstr + '_' + \
                    str(mu_top) + str(orig_var_top) + '_' + str(narcs) + \
                    '_' + str(seed) + str(params['d_top'])

                if experiment_name == 'varying_lams/':
                    save_extension = str(mu_top) + str(var_top) + '_' + str(narcs) + \
                        '_' + str(seed) + str(params['d_top'])

                if experiment_name == 'varying_lams/':

                    #     UG=read_and_create_graph(**params)
                    #     save(UG, 'graph' + save_extension, n)
                    # else:
                    #     UG = load('graph' + save_extension, n)
                    # UG=read_and_create_graph(**params)

                    if first:
                        print('first')
                        UG = create_random_graph(**params)
                        save(UG, 'graph' + save_extension, n)
                    else:
                        try:
                            UG = load('graph' + save_extension, n)
                            print('graph loaded')
                        except:
                            print('graph does not exists, creating now...')
                            UG = create_random_graph(**params)
                            save(UG, 'graph' + save_extension, n)
                    first = False

                    print('graph is ready')
                    params['vartop'] = var_top
                    params['G'] = UG

                    res_time = time.time()
                    UR = UG.build_res_network()
                    res_time = time.time() - res_time
                    res_time = 0
                    params['R'] = None

                    if lamBar == 0:
                        varcost = None
                        lamstr = str(round(lamBar, 6))
                        print(params['lambar'])

                        lamstr = lamstr.replace('.', '')
                        # discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa_lin(
                        #     UG, UR, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-3, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                        # mucost = UG.mu_cost()
                        # varcost = np.sqrt(UG.var_cost())
                        # print('mu and var: ', mucost, varcost)
                        # save(mucost, 'mu_cost_' + lamstr + '_' + save_extension , n, experiment_name)
                        # save(varcost, 'var_cost_' + lamstr + '_' + save_extension, n, experiment_name)
                        print('lam is 0')
                        cvx_mosek_obj, cvx_elapsed, mosek_overhead, x_msk, mucost, varcost = mosek_solve_lin(
                            **params)
                        print('mosek done in: ', cvx_elapsed)
                        print('mosek obj: ', cvx_mosek_obj)
                        print('mu and var: ', mucost, varcost)
                        save(mucost, 'mu_cost_' + lamstr + '_' +
                             save_extension, n, experiment_name)
                        save(varcost, 'var_cost_' + lamstr + '_' +
                             save_extension, n, experiment_name)

                    else:
                        print('lam not 0')
                        print(params['lambar'])
                        lamstr = str(round(lamBar, 6))

                        lamstr = lamstr.replace('.', '')
                        pdb.set_trace()
                        cvx_mosek_obj, cvx_elapsed, mosek_overhead, x_msk, mucost, varcost = cplex_solve(
                            **params)
                        print('mosek done in: ', cvx_elapsed)
                        print('mosek obj: ', cvx_mosek_obj)
                        print('mu and var: ', mucost, varcost)
                        save(mucost, 'mu_cost_' + lamstr + '_' +
                             save_extension, n, experiment_name)
                        save(varcost, 'var_cost_' + lamstr + '_' +
                             save_extension, n, experiment_name)
                    continue

                if experiment_name == 'varying_lams_2/':
                    print('cur_lam: {}, lam: {}'.format(cur_lam, lamBar))
                    lamstr = str(cur_lam / scale)
                    lamstr = lamstr.replace(".", "-")

                    #     UG=read_and_create_graph(**params)
                    #     save(UG, 'graph' + save_extension, n)
                    # else:
                    #     UG = load('graph' + save_extension, n)
                    # UG=read_and_create_graph(**params)

                    # if first:
                    #     print('first')
                    #     UG = create_random_graph(**params)
                    #     save(UG, 'graph' + save_extension, n)
                    # else:
                    try:
                        UG = load('graph' + save_extension, n)
                        print('graph loaded')
                    except:
                        print('graph does not exists, creating now...')
                        UG = create_random_graph(**params)
                        save(UG, 'graph' + save_extension, n)

                    first = False

                    multip = cur_lam / lamBar

                    vardict = nx.get_edge_attributes(UG.nxg, 'var')
                    for key in vardict:
                        vardict[key] = vardict[key] * multip**2

                    var_top = orig_var_top * multip**2

                    save_extension = lamstr + '_' + \
                        str(mu_top) + str(orig_var_top) + '_' + str(narcs) + \
                        '_' + str(seed) + str(params['d_top'])

                    nx.set_edge_attributes(UG.nxg, vardict, 'var')

                    params['vartop'] = var_top
                    params['G'] = UG

                    res_time = time.time()
                    UR = UG.build_res_network()
                    res_time = time.time() - res_time
                    res_time = 0
                    params['R'] = None

                    print('lam not 0')
                    print(params['lambar'])
                    params['cur_lam'] = cur_lam
                    lamstr = str(round(lamBar, 6))

                    lamstr = lamstr.replace('.', '')
                    cvx_mosek_obj, cvx_elapsed, mosek_overhead, x_msk, mucost, varcost = mosek_solve(
                        **params)
                    print('mosek done in: ', cvx_elapsed)
                    print('mosek obj: ', cvx_mosek_obj)

                    varcost = None
                    decoyG2 = copy.deepcopy(UG)
                    decoyR2 = copy.deepcopy(UR)

                    discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa_lin(
                        decoyG2, decoyR2, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-4, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                    del decoyG2
                    del decoyR2
                    del UR
                    lam_low = (lamBar / scale) / (2.0 * np.sqrt(varhigh))

                    lam_high = (lamBar / scale) / 2.0
                    lam = (lam_high + lam_low) / 2

                    print('lamhigh: ', lam_high)
                    print('lamlow: ', lam_low)
                    print('addtime: ', res_time + nmcc_time2)
                    params['start_varcost'] = None
                    params['add_time'] = res_time + nmcc_time2
                    params['lam_high'] = lam_high
                    params['lam_low'] = lam_low
                    params['startlam'] = lam

                    params['start_varcost_bs'] = None
                    params['add_time_bs'] = res_time + nmcc_time2
                    params['lam_high_bs'] = lam_high
                    params['startlam_bs'] = lam
                    params['max_iters'] = None

                    params['way'] = 'none'
                    nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                        **params)
                    print('nr done in: ', times_nr[-1])
                    print('nr obj: ', nr_obj)
                    continue
                    # params['way'] = 'newformul'
                    # nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                    #     **params)
                    # print('nr done in: ', times_nr[-1])
                    # print('nr obj: ', nr_obj)

                    # params['way'] = 'newxicost'
                    # nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                    #     **params)
                    # print('nr done in: ', times_nr[-1])
                    # print('nr obj: ', nr_obj)

                    # nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                    #     **params)
                    # print('nr done in: ', times_nr[-1])
                    # print('nr obj: ', nr_obj)

                    # alg3_obj, alg3_iters, alg3_subproblem_times,  times_alg3, alg3_objs, alg3_lams = alg3(
                    #     **params)
                    # print('alg3 done in: ', times_alg3[-1])
                    # print('alg3 obj: ', alg3_obj)

                    # bs_obj, bs_iters, bs_subproblem_times, times_bs, bs_objs, bs_lams, x_bs = bisect_solve(
                    #     **params)
                    # print('bs done in: ', times_bs[-1])
                    # print('bs obj: ', bs_obj)

                    save(cvx_mosek_obj, 'cvx_mosek_obj' +
                         save_extension, n, experiment_name)
                    save(mosek_overhead, 'mosek_overhead_time' +
                         save_extension, n, experiment_name)

                    save(times_bs, 'times_bs' +
                         save_extension, n, experiment_name)
                    save(bs_subproblem_times, 'bs_subproblem_times' +
                         save_extension, n, experiment_name)
                    save(bs_objs, 'bs_objs' +
                         save_extension, n, experiment_name)
                    save(bs_lams, 'bs_lams' +
                         save_extension, n, experiment_name)

                    save(times_nr, 'times_nr' +
                         save_extension, n, experiment_name)
                    save(nr_subproblem_times, 'nr_subproblem_times' +
                         save_extension, n, experiment_name)
                    save(nr_objs, 'nr_objs' +
                         save_extension, n, experiment_name)
                    save(nr_lams, 'nr_lams' +
                         save_extension, n, experiment_name)

                    save(times_alg3, 'times_alg3' +
                         save_extension, n, experiment_name)
                    save(alg3_subproblem_times, 'alg3_subproblem_times' +
                         save_extension, n, experiment_name)
                    save(alg3_objs, 'alg3_objs' +
                         save_extension, n, experiment_name)
                    save(alg3_lams, 'alg3_lams' +
                         save_extension, n, experiment_name)

                    print('node, lambda_bar, mu, var, arcno, seed: ',
                          n, cur_lam, mu_top, var_top, narcs, seeds[i])

                    print(tabulate([['NR', nr_obj, times_nr[-1], nr_iters], ['ALG3', alg3_obj, times_alg3[-1], alg3_iters], ['BS', bs_obj, times_bs[-1], bs_iters], [
                        'CVX', cvx_mosek_obj, cvx_elapsed, 1]], headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations']))

                    continue

                if experiment_name == 'varying_variance':
                    for vv in [1, 10, 100]:

                        try:
                            UG = load('graph' + save_extension, n)
                            print('graph loaded')
                        except:
                            print('graph does not exists, creating now...')
                            UG = create_random_graph(**params)
                            save(UG, 'graph' + save_extension, n)
                        print('graph is ready')

                        vardict = nx.get_edge_attributes(UG.nxg, 'var')
                        for key in vardict:
                            vardict[key] = vardict[key] * vv

                        var_top = orig_var_top * vv

                        save_extension = lamstr + '_' + \
                            str(mu_top) + str(var_top) + '_' + str(narcs) + \
                            '_' + str(seed) + str(params['d_top'])

                        nx.set_edge_attributes(UG.nxg, vardict, 'var')

                        print('node, lambda_bar, mu, var, arcno, seed: ',
                              n, lamBar, mu_top, var_top, narcs, seeds[i])
                        params['vartop'] = var_top
                        params['G'] = UG

                        # try:
                        #     UG = load('graph' + save_extension, n)
                        #     print('graph loaded')
                        # except:
                        print('graph does not exists, creating now...')
                        # UG = create_random_graph(**params)
                        save(UG, 'graph' + save_extension, n)
                        print('graph is ready')

                        res_time = time.time()
                        UR = UG.build_res_network()
                        res_time = time.time() - res_time
                        res_time = 0
                        params['R'] = None

                if fixmsk:
                    cvx_mosek_obj, cvx_elapsed, mosek_overhead = mosek_solve(
                        **params)
                    print('mosek done in: ', cvx_elapsed)
                    print('mosek obj: ', cvx_mosek_obj)
                    save(cvx_mosek_obj, 'cvx_mosek_obj' +
                         save_extension, n, experiment_name)
                    save(mosek_overhead, 'mosek_overhead_time' +
                         save_extension, n, experiment_name)
                else:

                    # varcost = None
                    # decoyG = copy.deepcopy(UG)
                    # decoyG.set_lambda(lam=0.1)
                    # decoyR = copy.deepcopy(UR)
                    # nx.set_edge_attributes(decoyG.nxg, 0, 'mu')
                    # discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
                    #     decoyG, decoyR, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-3, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                    # del decoyG
                    # del decoyR
                    # lam_high = (lamBar / scale) / (2 * np.sqrt(varcost))

                    # varhigh,d = mosek_solve_lin(**params)
                    # print('elapsed linear: ', d, varhigh)
                    # varcost = None
                    # decoyG2 = copy.deepcopy(UG)
                    # decoyG2.set_lambda(lam=0.7)
                    # decoyR2 = copy.deepcopy(UR)

                    # discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa(
                    # decoyG2, decoyR2, varcost, fully=True,
                    # tol=params['subproblem_tol'], lamsearch=False,
                    # nfullytol=1e-7, vartop=params['var_top'],
                    # mutop=params['mu_top'], difftol=1e-2)

                    # algo_obj = (0.7 * decoyG2.var_cost() + decoyG2.mu_cost())
                    # print(algo_obj)

                    # c, flow = capacity_scaling(UG.nxg, weight='mu', lam=0.7)
                    # feasible_flow = {}
                    # for i in UG.nxg.nodes():
                    #     for j in flow[i].keys():
                    #         if j != 't':
                    #             # if abs(round(flow[i][j],4) - flow[i][j]) <= 1e-10:
                    #             #     flow[i][j] = round(flow[i][j],4)
                    #             feasible_flow[i, j] = flow[i][j]

                    # UG.set_flow_attributes(feasible_flow)
                    # algo_obj = (0.7 * UG.var_cost() + UG.mu_cost())
                    # print('cap scaling: ', algo_obj)

                    # print('time and cost: ', nmcc_time2, varhigh)
                    #############

                    # decoyG2 = copy.deepcopy(UG)
                    # decoyG2.set_lambda(lam=0.1)
                    # decoyR2 = copy.deepcopy(UR)
                    # c,d = capacity_scaling(decoyG2.nxg, weight='mu')
                    # varcost = None
                    # m = UG.nxg.number_of_edges()
                    # orig = nx.get_edge_attributes(decoyG2.nxg,'mu')
                    # nx.set_edge_attributes(decoyG2.nxg, 0, 'var')
                    # discount, nmcc_time, scc_time, varcost = solve_mcf_sa(
                    #     decoyG2, decoyR2, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-2, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                    # nx.set_edge_attributes(decoyG2.nxg, orig, 'mu')

                    # varcost_h = m * (params['var_top']
                    #              * params['d_top'] * params['d_top'])
                    # lam_low = (lamBar / scale) / (2 * np.sqrt(varcost_h))

                    UG = create_random_graph(**params)
                    save(UG, 'graph' + save_extension, n)
                    params['vartop'] = var_top
                    params['G'] = UG

                    res_time = time.time()
                    UR = UG.build_res_network()
                    res_time = time.time() - res_time
                    res_time = 0
                    params['R'] = None

                    try:
                        cvx_mosek_obj, cvx_elapsed, mosek_overhead, x_msk, mucost, varcost = cplex_solve(
                            **params)
                        print('mosek done in: ', cvx_elapsed)
                        print('mosek obj: ', cvx_mosek_obj)
                    except:
                        continue

                    # continue
                    # cvx_elapsed=5
                    # cvx_mosek_obj=5
                    # mosek_overhead=5
                    # cvx_mosek_obj = 5
                    # cvx_elapsed = 5
                        # if lamBar ==0:
                        #     varcost = None

                        #     discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa_lin(
                        #         UG, UR, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-6, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                        #     mucost = UG.mu_cost()
                        #     varcost = np.sqrt(UG.var_cost())
                        #     print('mu and var: ', mucost, varcost)
                        #     save(mucost, 'mu_cost_' + lamstr + '_' + save_extension , n, experiment_name)
                        #     save(varcost, 'var_cost_' + lamstr + '_' + save_extension, n, experiment_name)

                        #     break

                        # else:
                        #     varcost = None
                        #     decoyG2 = copy.deepcopy(UG)
                        #     decoyR2 = copy.deepcopy(UR)

                        #     discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa_lin(
                        #         decoyG2, decoyR2, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-4, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                        #     del decoyG2
                        #     del decoyR2
                        #     del UR
                        #     lam_low = (lamBar / scale) / (2.0 * np.sqrt(varhigh))
                        #     lam_high = (lamBar/scale)/2.0
                        #     lam = (lam_high + lam_low) / 2

                        #     print('lamhigh: ', lam_high)
                        #     print('lamlow: ', lam_low)
                        #     print('addtime: ', res_time + nmcc_time2)
                        #     params['start_varcost'] = varcost
                        #     params['add_time'] = res_time + nmcc_time2
                        #     params['lam_high'] = lam_high
                        #     params['lam_low'] = lam_low
                        #     params['startlam'] = lam

                        #     params['start_varcost_bs'] = varcost
                        #     params['add_time_bs'] = res_time + nmcc_time2
                        #     params['lam_high_bs'] = lam_high
                        #     params['startlam_bs'] = lam
                        #     params['max_iters'] = None

                        #     bs_obj, bs_iters, bs_subproblem_times, times_bs, bs_objs, bs_lams, x_bs,mucost, varcost = bisect_solve_exp(
                        #         **params)

                        #     print('mu and var: ', mucost, varcost)
                        #     save(mucost, 'mu_cost_' + lamstr + '_' + save_extension , n, experiment_name)
                        #     save(varcost, 'var_cost_' + lamstr + '_' + save_extension, n, experiment_name)
                        #     break

                    varcost = None
                    decoyG2 = copy.deepcopy(UG)
                    decoyR2 = copy.deepcopy(UR)

                    discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa_lin(
                        decoyG2, decoyR2, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-4, vartop=params['var_top'], mutop=params['mu_top'], difftol=1e-2)
                    del decoyG2
                    del decoyR2
                    del UR
                    lam_low = (lamBar / scale) / (2.0 * np.sqrt(varhigh))
                    # nmcc_time2 = 0
                    # lam_low=0
                    lam_high = (lamBar / scale) / 2.0
                    lam = (lam_high + lam_low) / 2

                    print('lamhigh: ', lam_high)
                    print('lamlow: ', lam_low)
                    print('addtime: ', res_time + nmcc_time2)
                    params['start_varcost'] = None
                    params['add_time'] = res_time + nmcc_time2
                    params['lam_high'] = lam_high
                    params['lam_low'] = lam_low
                    params['startlam'] = lam

                    params['start_varcost_bs'] = None
                    params['add_time_bs'] = res_time + nmcc_time2
                    params['lam_high_bs'] = lam_high
                    params['startlam_bs'] = lam
                    params['max_iters'] = None

                    if fix_bs:
                        bs_obj, bs_iters, bs_subproblem_times, times_bs, bs_objs, bs_lams, x_bs = bisect_solve(
                            **params)
                        save(times_bs, 'times_bs' +
                             save_extension, n, experiment_name)
                        save(bs_subproblem_times, 'bs_subproblem_times' +
                             save_extension, n, experiment_name)
                        save(bs_objs, 'bs_objs' +
                             save_extension, n, experiment_name)
                        save(bs_lams, 'bs_lams' +
                             save_extension, n, experiment_name)

                        print('bs done in: ', times_bs[-1])
                        print('bs obj: ', bs_obj)
                        # x_msk = np.array(list(x_msk.values()))
                        # x_bs = np.array(list(x_bs.values()))
                        pdb.set_trace()
                        print('max indiv difference: ', max(abs(x_msk - x_bs)))
                        max_ind = np.argmax(abs(x_msk - x_bs))
                        print('max diff corresponding value for msk {}, for bs {}'.format(
                            x_msk[max_ind], x_bs[max_ind]))
                    else:

                        params['way'] = 'oldformul'
                        nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                            **params)
                        print('nr done in: ', times_nr[-1])
                        print('nr obj: ', nr_obj)

                        params['way'] = 'newformul'
                        nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                            **params)
                        print('nr done in: ', times_nr[-1])
                        print('nr obj: ', nr_obj)

                        params['way'] = 'newxicost'
                        nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
                            **params)
                        print('nr done in: ', times_nr[-1])
                        print('nr obj: ', nr_obj)

                        bs_obj, bs_iters, bs_subproblem_times, times_bs, bs_objs, bs_lams, x_bs = bisect_solve(
                            **params)
                        print('bs done in: ', times_bs[-1])
                        print('bs obj: ', bs_obj)

                        alg3_obj, alg3_iters, alg3_subproblem_times,  times_alg3, alg3_objs, alg3_lams = alg3(
                            **params)
                        print('alg3 done in: ', times_alg3[-1])
                        print('alg3 obj: ', alg3_obj)

                        print('================||||||||=================')
                        save(cvx_mosek_obj, 'cvx_mosek_obj' +
                             save_extension, n, experiment_name)
                        save(mosek_overhead, 'mosek_overhead_time' +
                             save_extension, n, experiment_name)

                        save(times_bs, 'times_bs' +
                             save_extension, n, experiment_name)
                        save(bs_subproblem_times, 'bs_subproblem_times' +
                             save_extension, n, experiment_name)
                        save(bs_objs, 'bs_objs' +
                             save_extension, n, experiment_name)
                        save(bs_lams, 'bs_lams' +
                             save_extension, n, experiment_name)

                        save(times_nr, 'times_nr' +
                             save_extension, n, experiment_name)
                        save(nr_subproblem_times, 'nr_subproblem_times' +
                             save_extension, n, experiment_name)
                        save(nr_objs, 'nr_objs' +
                             save_extension, n, experiment_name)
                        save(nr_lams, 'nr_lams' +
                             save_extension, n, experiment_name)

                        save(times_alg3, 'times_alg3' +
                             save_extension, n, experiment_name)
                        save(alg3_subproblem_times, 'alg3_subproblem_times' +
                             save_extension, n, experiment_name)
                        save(alg3_objs, 'alg3_objs' +
                             save_extension, n, experiment_name)
                        save(alg3_lams, 'alg3_lams' +
                             save_extension, n, experiment_name)
                        print('node, lambda_bar, mu, var, arcno, seed: ',
                              n, lamBar, mu_top, var_top, narcs, seeds[i])

                        print(tabulate([['NR', nr_obj, times_nr[-1], nr_iters], ['ALG3', alg3_obj, times_alg3[-1], alg3_iters], ['BS', bs_obj, times_bs[-1], bs_iters], [
                            'CVX', cvx_mosek_obj, cvx_elapsed, 1]], headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations']))

params = {}
lambar = 0.9
G = MCF_DiGraph(lambar)
G.nxg = nx.DiGraph()

filename = 'networks/netgen_8_12a.min'
if filename.find('goto') >=0:
    generator = 'goto'
else:
    generator = 'netgen' 

G = load_data(networkFileName=filename, G=G, generator=generator)

solver_obj, solver_elapsed, soln1 = cvxpy_solve(G)
print('solver')
print(solver_obj, solver_elapsed)

obj, bound_elapsed, x, prob = cvxpy_solve_lin(G, 0, prob=None, warm_start=False)
var_cost = x.dot(np.multiply(G.var,np.eye(G.m))).dot(x)
low = G.lambar/(2.0*np.sqrt(var_cost))

# elapsed, xi = cvxpy_solve_xi(G, soln)
# low = G.lambar*(np.multiply(G.var,soln).dot(xi))/(np.sqrt(var_cost)*sum(xi))

obj, bound_elapsed_2, x, prob = cvxpy_solve_lin(G, 1000, prob=None, warm_start=False)
var_cost = x.dot(np.multiply(G.var,np.eye(G.m))).dot(x)
high = G.lambar/(2.0*np.sqrt(var_cost))

# elapsed2, xi = cvxpy_solve_xi(G, soln)
# high = G.lambar*(np.multiply(G.var,soln).dot(xi))/(np.sqrt(var_cost)*sum(xi))

bound_elapsed = bound_elapsed+bound_elapsed_2
print(high, low, bound_elapsed)

bs_obj, bs_elapsed, soln2 = bs_cvxpy(G, high=high, low=low, prob=prob)
print('bs')
print(bs_obj, bs_elapsed)
print('total elapsed: ', bs_elapsed + bound_elapsed)
pdb.set_trace()

nr_obj, nr_elapsed = nr_cvxpy(G)
print('nr')
print(nr_obj, nr_elapsed)





pdb.set_trace()




scale = 1
params['scale'] = 1
params['lambar'] = lambar / scale
params['sensitivity_tol'] = 1e-3
params['stop_tol'] = 1e-5
params['precision_start_tol'] = 1e-2
params['precision_start_tol_bs'] = 1e-2
params['subproblem_tol'] = 1e-1
params['subproblem_tol_bs'] = 1e-1
params['cvx_obj'] = None
params['mu_scalar'] = 1.0 / scale
params['var_scalar'] = 1.0  # scale*scale
params['obj_scalar'] = scale
params['varcostnoful'] = 1e0
params['varcostful'] = 1e0
params['xicostnoful'] = 1e0


res_time = time.time()


R = G.build_res_network()
res_time = time.time() - res_time
params['G'] = G
params['R'] = None

varcost = None
decoyG2 = copy.deepcopy(G)
decoyR2 = copy.deepcopy(R)

discount, nmcc_time2, scc_time, varhigh = solve_mcf_sa_lin(
    decoyG2, decoyR2, varcost, fully=True, tol=params['subproblem_tol'], lamsearch=False, nfullytol=1e-4, vartop=G.vartop, mutop=G.mutop, difftol=1e-2)
del decoyG2
del decoyR2
del UR
lam_low = (lamBar / scale) / (2.0 * np.sqrt(varhigh))
# nmcc_time2 = 0
# lam_low=0
lam_high = (lamBar / scale) / 2.0
lam = (lam_high + lam_low) / 2

print('lamhigh: ', lam_high)
print('lamlow: ', lam_low)
print('addtime: ', res_time + nmcc_time2)
params['start_varcost'] = None
params['add_time'] = res_time + nmcc_time2
params['lam_high'] = lam_high
params['lam_low'] = lam_low
params['startlam'] = lam

params['start_varcost_bs'] = None
params['add_time_bs'] = res_time + nmcc_time2
params['lam_high_bs'] = lam_high
params['startlam_bs'] = lam
params['max_iters'] = None
params['way'] = 'newformul'
nr_obj, nr_iters, nr_subproblem_times, times_nr, nr_objs, nr_lams = NR_solve(
    **params)
print('nr done in: ', times_nr[-1])
print('nr obj: ', nr_obj)


