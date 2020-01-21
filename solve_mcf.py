from MCF_DiGraph import *
from utils import *
import matplotlib.pyplot as plt
import pandas
import time
import random
import itertools
import copy
from tabulate import tabulate
import cProfile
import pstats
import scipy
# sys.path.append('/usr/local/Cellar/graph-tool/2.27_1/lib/python3.7/site-packages/')
# import graph_tool as gt


# Solvers
import mosek.fusion as mf
from mosek.fusion import Expr, Set, Domain, Var, ObjectiveSense, Matrix
import cvxpy as cp
# from cvxopt import matrix, solvers, spmatrix
# import cplex
# import docplex
# from gurobipy import *


graph_loc = '/Users/cgokalp/repos/dev/msmcf/residual_graph.lgf'
nmcc_loc = '/Users/cgokalp/repos/dev/msmcf/msmcf/nmcc.txt'


def solve_xi_mosek(G, lam):

    cvx_time_st = time.time()
    G.set_lambda(lam=lam)

    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    x_ = np.zeros(m)
    mu = np.zeros(m)
    var = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m, n))
    A = np.zeros((n, m))

    rows = []
    cols = []
    values = []
    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu', 0)
        # sigma[i] = np.sqrt(e.get('var', 0))
        var[i] = e.get('var', 0)
        x_[i] = e.get('flow', 0)
        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)
        # A[u-1, i] = 1
        # A[v-1, i] = -1
        # A[u, i] = 1
        # A[v, i] = -1
        # rows.append(u - 1)
        rows.append(u)

        cols.append(i)
        values.append(1)
        # rows.append(v - 1)
        rows.append(v)

        cols.append(i)
        values.append(-1)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

        arc_data = {'Start': [], 'End': []}

    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)

    x = model.variable("x", m)
    model.constraint(x, Domain.lessThan(cap))

    # x_n = model.variable("x_n", m)
    # decoy = model.variable("decoy", 1)
    # decoy_2 = Var.vstack(decoy, x_n)

    model.objective("myobj",
                    ObjectiveSense.Minimize, Expr.dot(mu, x))
    # model.constraint(Expr.sub(x_n, Expr.mulElm(sigma, x)),
    #                  Domain.equalsTo(0.0))
    # model.constraint(decoy_2, Domain.inQCone())
    A = Matrix.sparse(n, m, rows, cols, values)
    model.constraint(Expr.mul(A, x), Domain.equalsTo(0.0))

    constraints = []
    for i in range(m):
        if abs(x_[i]) < 1e-9:
            constraints.append(xi[i] == 0)
            # constraints.append(xi[i]<=cap[i])

        elif x_[i] == cap[i]:
            constraints.append(-x_[i] <= xi[i])
            constraints.append(xi[i] <= 0)
        else:
            constraints.append(-x_[i] <= xi[i])
            constraints.append(xi[i] <= cap[i] - x_[i])

    for i in range(n):
        preds = []
        succs = []
        for key, vals in arc_dict.items():
            if vals[0] == i:
                preds.append(key)
            if vals[1] == i:
                succs.append(key)

        constraint = sum(xi[p] for p in preds) - sum(xi[s] for s in succs) == 0
        constraints.append(constraint)

    # Construct the problem.
    objective = cp.Minimize(G.lam * var.T * xi**2 +
                            2 * np.multiply(var, x_).T * xi)
    # print(sum(np.multiply(var,x_)))
    # print(sum(x_))
    prob = cp.Problem(objective, constraints)
    cvx_time_st = time.time()

    result = prob.solve(solver='MOSEK', verbose=False)
    print("status:", prob.status)

    cvx_soln = xi.value
    cvx_obj = objective.value
    cvx_elapsed = time.time() - cvx_time_st
    print(cvx_elapsed)
    # cost = 0
    # for u,v,e in self.nxg.edges(data=True):
    #     cost += e.get('xi', 0) * e.get('var', 0) * e.get('flow', 0)
    # return cost

    return cvx_soln, np.multiply(var, x_).dot(cvx_soln)


def mosek_solve_lin(**kwargs):

    if kwargs['G'] != None:
        G = copy.deepcopy(kwargs['G'])
    else:
        G = create_random_graph(**kwargs)
    lambar = kwargs['lambar']

    G.set_lambda(lam=lambar)
    cvx_time_st = time.time()

    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    mu = np.zeros(m)
    var = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m, n))
    A = np.zeros((n, m))

    rows = []
    cols = []
    values = []
    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu', 0)
        # sigma[i] = np.sqrt(e.get('var', 0))
        var[i] = e.get('var', 0)

        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)
        # A[u-1, i] = 1
        # A[v-1, i] = -1
        # A[u, i] = 1
        # A[v, i] = -1
        # rows.append(u - 1)
        rows.append(u)

        cols.append(i)
        values.append(1)
        # rows.append(v - 1)
        rows.append(v)

        cols.append(i)
        values.append(-1)
        i += 1
    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

    with mf.Model() as model:

        const_time = time.time()

        x = model.variable("x", m, Domain.greaterThan(0.0))
        # x_n = model.variable("x_n", m)
        # decoy = model.variable("decoy", 1)
        # decoy_2 = Var.vstack(decoy, x_n)

        model.objective("myobj",
                        ObjectiveSense.Minimize, Expr.dot(mu, x))
        # model.constraint(Expr.sub(x_n, Expr.mulElm(sigma, x)),
        #                  Domain.equalsTo(0.0))
        model.constraint(x, Domain.lessThan(cap))
        # model.constraint(decoy_2, Domain.inQCone())
        A = Matrix.sparse(n, m, rows, cols, values)
        model.constraint(Expr.sub(Expr.mul(A, x), d), Domain.equalsTo(0.0))

        # if kwargs['max_iters'] == None:
        # model.setSolverParam("intpntCoTolRelGap", 1.0e-10)
        # model.setSolverParam("intpntCoTolPfeas", 1.0e-10)
        # model.setSolverParam("intpntCoTolDfeas", 1.0e-10)
        # model.setSolverParam("intpntCoTolMuRed", 1.0e-10)
        # model.setSolverParam("intpntMaxIterations", 100000)
        # else:
        #     model.setSolverParam("intpntMaxIterations", kwargs['max_iters'])

        # import logging

        # logging = Logger(**kwargs)
        # model.setLogHandler(logging)

        # model.setSolverParam("logIntpnt", 1000)
        # model.setSolverParam("log", 1000)
        # model.setSolverParam("logFile", 1000)

        # tm = model.getSolverDoubleInfo("optimizerTime")
        # om = model.getSolverDoubleInfo("intpntPrimalObj")
        # it = model.getSolverIntInfo("intpntIter")
        # print('Time: {0}\nIterations: {1}\nObjective {2}'.format(tm, it, om))
        const_time = time.time() - const_time
        solve_time = time.time()

        model.solve()

        solve_time = time.time() - solve_time
        mosek_overhead = const_time
        # tm = model.getSolverDoubleInfo("optimizerTime")
        # om = model.getSolverDoubleInfo("intpntPrimalObj")
        # it = model.getSolverIntInfo("intpntIter")
        # print('Time: {0}\nIterations: {1}\nObjective {2}'.format(tm, it, om))
        # pdb.set_trace()
        # print(model.primalObjValue())
        # print('const_time: ', const_time)
        # print('total_time_elapsed: ', solve_time + const_time)
        # print(model.getProblemStatus())

        cvx_obj = model.primalObjValue()
        # x = x.level()
        # try:
        #     i = 0
        #     for u, v, e in G.nxg.edges(data=True):
        #         cur_x = x[i]
        #         if abs(round(cur_x, 4) - cur_x) <= 1e-8:
        #             cur_x = round(cur_x, 4)
        #         G.nxg[u][v]['flow'] = cur_x
        #         x[i] = cur_x
        #         i += 1

        #     cur_r_dict = nx.to_dict_of_dicts(G.nxg)
        #     cc = [v for v in cur_r_dict.keys()]
        #     list_ofed = []
        #     endnode = len(G.nxg.nodes())
        #     for c in cc:
        #         if endnode in cur_r_dict[c].keys():
        #             list_ofed.append(c)

        #     flows = []
        #     for e in list_ofed:
        #         flows.append(G.nxg[e][endnode]['flow'])
        #     flows = sum(np.array(flows))
        # except:
        #     pass

        # x = nx.get_edge_attributes(G.nxg, 'flow')
        # cvx_obj = G.tot_cost()
        cvx_elapsed = solve_time + const_time

        # cvx_obj = model.primalObjValue()
        # x = x.level()
        # try:
        #     i = 0
        #     for u, v, e in G.nxg.edges(data=True):
        #         cur_x = x[i]
        #         if abs(round(cur_x, 4) - cur_x) <= 1e-8:
        #             cur_x = round(cur_x, 4)
        #         G.nxg[u][v]['flow'] = cur_x
        #         x[i] = cur_x
        #         i += 1

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
        # flows = sum(np.array(flows))
        # print('going to the last node: ', flows)
        # except:
        #     pass

        # var = G.var_cost()
        # cvx_elapsed = solve_time + const_time

        return cvx_obj, cvx_elapsed, mosek_overhead, x.level(), mu.dot(x.level()), np.sqrt(var.dot(x.level()**2))


def gurobipy_solve(**kwargs):
    if kwargs['G'] != None:
        G = copy.deepcopy(kwargs['G'])
    else:
        G = create_random_graph(**kwargs)
    lambar = kwargs['lambar']
    G.set_lambda(lam=lambar)

    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    mu = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m, n))

    i = 0
    arc_dict = {}
    A = np.zeros((n, m))

    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu', 0)
        sigma[i] = np.sqrt(e.get('var', 0))
        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)
        A[u - 1, i] = 1
        A[v - 1, i] = -1
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1
    model = Model('gurobi')

    # Create variables

    A = scipy.sparse.csc_matrix((values, (rows, cols)))

    x = model.addVars(m, name="x")
    x_n = model.addVars(m, name="x_n")
    decoy = model.addVar()
    expr2 = lambar * decoy

    # obj

    mobj = sum(x[k] * mu[k] for k in range(len(x)))
    model.setObjective((kwargs['mu_scalar'] * mobj + expr2)
                       * kwargs['obj_scalar'], GRB.MINIMIZE)

    # Arc capacity constraints
    model.addConstrs((x[i] <= cap[i] for i in range(m)), "upper")
    model.addConstrs((0 <= x[i] for i in range(m)), "lower")
    model.addConstrs((x_n[i] == sigma[i] * x[i]
                      for i in range(m)), "sigma_thing")
    model.addConstr(decoy * decoy >=
                    (sum(x_n[k] * x_n[k] for k in range(len(x_n)))))

    constraint_time = time.time()

    model.addMConstrs(A * x == d)

    const_time = time.time() - constraint_time

    # Compute optimal solution
    if kwargs['max_iters'] == None:
        model.Params.BarQCPConvTol = 1e-15
        model.Params.BarIterLimit = 150
        model.Params.BarConvTol = 1e-15
        model.Params.FeasibilityTol = 1e-9
    # else:
    #     model.Params.BarQCPConvTol = 1e-9
    #     model.Params.BarConvTol = 1e-9
    #     model.Params.FeasibilityTol = 1e-9
    #     model.Params.BarIterLimit = kwargs['max_iters']

    model.Params.OutputFlag = False
    cvx_time_st = time.time()
    model.update()

    model.optimize()

    # Print solution
    # print(model.status)
    cvx_obj = round(model.ObjVal, kwargs['cutoff'])
    cvx_elapsed = time.time() - cvx_time_st
    return cvx_obj, cvx_elapsed, const_time


def cvxpy_solve(G):


    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    mu = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m, n))
    x = cp.Variable(m)
    theta = cp.Variable(1)
    # x_n=cp.Variable(m)
    rows = []
    values = []
    cols = []
    var = np.zeros(m)
    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu', 0)
        sigma[i] = np.sqrt(e.get('var', 0))
        var[i] = e.get('var', 0)

        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)
        rows.append(u - 1)
        # rows.append(u)
        cols.append(i)
        values.append(1)
        rows.append(v - 1)
        # rows.append(v)
        cols.append(i)
        values.append(-1)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

    kwargs = {}
    kwargs['obj_scalar'] = 1
    kwargs['mu_scalar'] = 1
    kwargs['num_nodes'] = n 
    kwargs['lambar'] = 0.5
    kwargs['num_arcs'] = m
    kwargs['mu_top'] = 1
    kwargs['var_top'] = 1
    kwargs['seed'] = 1

    lambar = G.lambar

    A = scipy.sparse.csc_matrix((values, (rows, cols)))

    constraints = [0 <= x, x <= cap, A@x == d]

    # objective=cp.Minimize(
    # (kwargs['mu_scalar'] * mu.T * x + lambar * cp.norm(x_n, 2)) * kwargs['obj_scalar'])
    objective = cp.Minimize((kwargs['mu_scalar'] * mu.T * x + lambar *
                             cp.norm(cp.multiply(sigma, x), 2)) * kwargs['obj_scalar'])
    prob = cp.Problem(objective, constraints)

    # print(objective.value)
    cvx_time_st = time.time()
    # 'SCS','ECOS','CVXOPT' - 'MOSEK', 'GUROBI', 'CPLEX'
    result = prob.solve(solver='MOSEK', verbose=True)  # gurobi mosek compare
    print(objective.value)
    # prob.unpack_results(cvx.ECOS, solver_output)
    # result = prob.solve(solver='MOSEK', verbose=True, mosek_params={'MSK_DPAR_INTPNT_CO_TOL_REL_GAP':1e-20,'MSK_DPAR_INTPNT_CO_TOL_INFEAS':1e-30,'MSK_IPAR_INTPNT_MAX_ITERATIONS':1000})
    cvx_soln = x.value
    cvx_obj = objective.value
    cvx_elapsed = time.time() - cvx_time_st
    print(prob.solver_stats.solve_time, prob.solver_stats.setup_time)
    cvx_elapsed = prob.solver_stats.solve_time
    pdb.set_trace()
   
    return cvx_obj, cvx_elapsed, cvx_soln


def mosek_solve(G):

    start_time = time.time()
    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    mu = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m, n))
    A = np.zeros((n, m))

    rows = []
    cols = []
    values = []
    i = 0
    arc_dict = {}


    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu', 0)
        sigma[i] = np.sqrt(e.get('var', 0))
        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)
        rows.append(u - 1)
        # rows.append(u)
        cols.append(i)
        values.append(1)
        rows.append(v - 1)
        # rows.append(v)
        cols.append(i)
        values.append(-1)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

    kwargs = {}
    kwargs['obj_scalar'] = 1
    kwargs['mu_scalar'] = 1
    kwargs['num_nodes'] = n 
    kwargs['lambar'] = 0.5
    kwargs['num_arcs'] = m
    kwargs['mu_top'] = 1
    kwargs['var_top'] = 1
    kwargs['seed'] = 1

    lambar = G.lambar
    with mf.Model() as model:

        const_time = time.time()

        x = model.variable("x", m, Domain.greaterThan(0.0))
        x_n = model.variable("x_n", m)
        decoy = model.variable("decoy", 1)
        decoy_2 = Var.vstack(decoy, x_n)

        model.objective("myobj",
                        ObjectiveSense.Minimize,
                        Expr.mul(kwargs['obj_scalar'],
                                 Expr.add(Expr.mul(kwargs['mu_scalar'], Expr.dot(mu, x)), Expr.mul(lambar, decoy))))
        model.constraint(Expr.sub(x_n, Expr.mulElm(sigma, x)),
                         Domain.equalsTo(0.0))
        model.constraint(x, Domain.lessThan(cap))
        model.constraint(decoy_2, Domain.inQCone())
        A = Matrix.sparse(n, m, rows, cols, values)
        model.constraint(Expr.sub(Expr.mul(A, x), d), Domain.equalsTo(0.0))

        # model.setSolverParam("intpntCoTolRelGap", 1.0e-9)
        # model.setSolverParam("intpntCoTolPfeas", 1.0e-9)
        # model.setSolverParam("intpntCoTolDfeas", 1.0e-9)
        # model.setSolverParam("intpntCoTolMuRed", 1.0e-9)
        # model.setSolverParam("intpntMaxIterations", 100000)

        import logging

        logging = Logger(**kwargs)
        model.setLogHandler(logging)

        # model.setSolverParam("logIntpnt", 1000)
        # model.setSolverParam("log", 1000)
        # model.setSolverParam("logFile", 1000)

        # tm = model.getSolverDoubleInfo("optimizerTime")
        # om = model.getSolverDoubleInfo("intpntPrimalObj")
        # it = model.getSolverIntInfo("intpntIter")
        # print('Time: {0}\nIterations: {1}\nObjective {2}'.format(tm, it, om))
        const_time = time.time() - const_time

        solve_time = time.time()

        model.solve()

        solve_time = time.time() - solve_time

        mosek_overhead = const_time
        # tm = model.getSolverDoubleInfo("optimizerTime")
        # om = model.getSolverDoubleInfo("intpntPrimalObj")
        # it = model.getSolverIntInfo("intpntIter")
        # print('Time: {0}\nIterations: {1}\nObjective {2}'.format(tm, it, om))
        # pdb.set_trace()
        # print(model.primalObjValue())
        # print('const_time: ', const_time)
        # print('total_time_elapsed: ', solve_time + const_time)
        # print(model.Model.getProblemStatus())

        obj = model.primalObjValue()
        x = x.level()
        # pdb.set_trace()
        # try:
        #     i = 0
        #     for u, v, e in G.nxg.edges(data=True):
        #         cur_x = x[i]
        #         if abs(round(cur_x, 4) - cur_x) <= 1e-8:
        #             cur_x = round(cur_x, 4)
        #         G.nxg[u][v]['flow'] = cur_x
        #         x[i] = cur_x
        #         i += 1

        #     cur_r_dict = nx.to_dict_of_dicts(G.nxg)
        #     cc = [v for v in cur_r_dict.keys()]
        #     list_ofed = []
        #     endnode = len(G.nxg.nodes())
        #     for c in cc:
        #         if endnode in cur_r_dict[c].keys():
        #             list_ofed.append(c)

        #     flows = []
        #     for e in list_ofed:
        #         flows.append(G.nxg[e][endnode]['flow'])
        #     flows = sum(np.array(flows))
        # except:
        #     pass

        # x = nx.get_edge_attributes(G.nxg, 'flow')
        # cvx_obj = G.tot_cost()
        elapsed = solve_time + const_time

        return obj, elapsed, mosek_overhead, x, mu.dot(x), decoy.level()[0]


def solve_xi_mmc(G, R, xicost, discount=0, nmcc_tot_ms=0, scc_time=0, tol=1e-3, fully=False, nfullytol=1e-3, difftol=1e-5, xicostnoful=1e0, vartop=1e2):
    prev_augment_amount = -1.0
    nmcc_exists = True
    iters = 0
    consistent = False
    consistent2 = False
    consistent3 = False
    here2 = -10
    prev_xi_cost = 10
    here = -10
    while nmcc_exists:
        iters += 1
        st_write = time.time()

        strconnected = True
        # scc_time_start = time.time()
        # try:
        #     scc = max(nx.strongly_connected_components(G), key=len)
        #     strconnected = False
        # except:
        #     pass
        # scc_time += time.time() - scc_time_start

        if not strconnected:
            # scc_time_start = time.time()
            # scc = [c for c in sorted(nx.strongly_connected_components(R.nxg),key=len, reverse=True)]
            # scc_time += time.time() - scc_time_start
            # scc_count = len(scc)
            # maxim = -np.inf
            # for i in range(scc_count):
            #     lenlist = len(scc[i])
            #     if lenlist >  maxim:
            #         maxim = max(maxim, lenlist)
            #         max_i = i
            # nodelist = scc[max_i]
            nodelist = scc
            # pdb.set_trace()

            cur_g = R.nxg.copy()
            allnodes = np.arange(len(R.nxg.nodes)) + 1
            removals = list(set(allnodes) - set(nodelist))
            cur_g.remove_nodes_from(removals)

                # manual
                # R.DG_2_dict(G)
                # if R.nmcc_exists:
                #     R.set_nmcc(G)
                #     nmcc = R.get_nmcc()
                # else:
                #     nmcc = []

            output_graph(filename=graph_loc, g=cur_g, G=G, ptype='xi')
            solve_nmcc()
            nmcc = read_nmcc(filename=nmcc_loc)  
            
            resid_nodelist = nodelist
            nmcc_decoy = []
            for n in nmcc:
                nmcc_decoy.append(resid_nodelist[n])
            nmcc = nmcc_decoy
        else:
            output_graph(filename=graph_loc, g=R, G=G, ptype='xi')
            solve_nmcc()
            nmcc = read_nmcc(filename=nmcc_loc)  

        if len(nmcc) < 3:
            break
        # if float(nmcc_cost)/vartop <= -1e-7 :
        if float(nmcc_cost) >= 0:
            # print('nmcc cost is bigger than 0')
            break

        R.set_nmcc_lemon(nmcc)
        G.set_nmcc_lemon(nmcc)
        nmcc_tot_ms += nmcc_time_ms
        discount += time.time() - st_write

        delta = R.find_delta()
        xsi_star = R.find_xsi_star_xi(G)
        augment_amount = min(xsi_star, delta)
        print(augment_amount)

        # if fully:
        # if abs(augment_amount) <= nfullytol and xsi_star <= nfullytol:
        if abs(augment_amount) <= nfullytol and xsi_star <= nfullytol:

            if not consistent:
                here = iters
                consistent = True
            elif consistent and not consistent2:
                if here == iters - 1:
                    consistent2 = True
                else:
                    consistent = False
            elif consistent and consistent2:
                if here == iters - 2:
                    print('xsi_star very small')
                    break
                else:
                    consistent = False
                    consistent2 = False
        # # if xicost ==0:
        # #     break
        # if iters >= 20:
        #     if iters%20==0:
        #         xicost=G.xi_cost()
        #         print('xicost: ', xicost)
        #         if xicost == 0:
        #             break
        #         if abs(100*(prev_xi_cost/xicost -1)) < difftol:
        #             if consistent3 == False:
        #                 here2 = iters
        #                 consistent3 = True
        #             elif consistent3 == True and here2 == iters-20:
        #                 break
        #             else:
        #                 consistent3 = False
        #         prev_xi_cost = xicost

        #     if iters%10==0:
        #         xicost=G.xi_cost()
        #         if abs(np.sqrt(prev_xicost) - np.sqrt(xicost)) < xicostnoful:
        #             if consistent == False:
        #                 here = iters
        #                 consistent = True
        #             elif consistent == True and here == iters-10:
        #                 break
        #             else:
        #                 consistent = False
        #         prev_xicost = xicost
        #
        #
        # else:
        #     if abs(augment_amount) < tol:
        #         break
        #     if iters%5==0:
        #         xicost=G.xi_cost()
        #         if abs(100*(np.array(prev_xicost)/(xicost) -1)) < 1e-1:
        #
        #         # if abs(np.sqrt(prev_xicost) - np.sqrt(xicost)) < xicostnoful:
        #             if consistent == False:
        #                 here = iters
        #                 consistent = True
        #             elif consistent == True and here == iters-5:
        #                 break
        #             else:
        #                 consistent = False
        #         prev_xicost = xicost

        R.augment_flow(augment_amount)
        G.adjust_flow_xi(nmcc, augment_amount)
    # print('xi iters: ', iters)
    xicost = G.xi_cost()
    nmcc_time = nmcc_tot_ms * 0.001
    return discount, nmcc_time, scc_time, xicost


def solve_mcf_sa(G, R, varcost=None, discount=0, nmcc_tot_ms=0, scc_time=0, difftol=1e-6, tol=1e-2, fully=False, nfullytol=5e-4, lamsearch=False, var_cost_noful_tol=5e0, var_cost_ful=5e0, nr_run=False, vartop=1, mutop=1e2):
    prev_augment_amount = -1.0
    nmcc_exists = True
    iters = 0
    consistent = False
    consistent2 = False
    consistent3 = False
    here = -10
    here2 = -10
    prev_var_cost = 10

    while nmcc_exists:
        iters += 1
        st_write = time.time()

        strconnected = True

        # scc_time_start = time.time()

        # can=nx.is_strongly_connected(R.nxg)
        # if not can:
        # strconnected = False
        # print('not connected')
        # try:
        #     scc = max(nx.strongly_connected_components(G), key=len)
        #
        # except:
        #     pass
        # scc_time += time.time() - scc_time_start

        if not strconnected:
            scc_time_start = time.time()
            scc = [c for c in sorted(
                nx.strongly_connected_components(R.nxg), key=len, reverse=True)]
            scc_time += time.time() - scc_time_start
            scc_count = len(scc)
            maxim = -np.inf
            for i in range(scc_count):
                lenlist = len(scc[i])
                if lenlist > maxim:
                    maxim = max(maxim, lenlist)
                    max_i = i
            nodelist = scc[max_i]
            cur_g = copy.deepcopy(R.nxg)
            allnodes = np.arange(len(R.nxg.nodes)) + 1
            removals = list(set(allnodes) - set(nodelist))
            cur_g.remove_nodes_from(removals)

            output_graph(filename=graph_loc, g=cur_g, G=G, ptype='sa')
            solve_nmcc()
            nmcc = read_nmcc(filename=nmcc_loc)  
        else:
            output_graph(filename=graph_loc, g=R, G=G, ptype='sa')
            solve_nmcc()
            nmcc = read_nmcc(filename=nmcc_loc)   

        if len(nmcc) < 3:
            print('len nmcc is 2')
            break

        divider = 10**(len(str(vartop)) - 1)

        if float(nmcc_cost) >= 0:  # *1e4/divider >= -1e-4:
            # print('nmcc cost bigger than 0')
            break

        R.set_nmcc_lemon(nmcc)
        G.set_nmcc_lemon(nmcc)
        nmcc_tot_ms += nmcc_time_ms
        discount += time.time() - st_write

        delta = R.find_delta()
        xsi_star = R.find_xsi_star(G)

        augment_amount = min(xsi_star, delta)
        # print(nmcc_cost, delta, xsi_star)
        # if not lamsearch:
        #     if fully:
        # print('nmcc_cost is: {}, delta is: {}, and xsi_star is: {}'.format(
        #                 float(nmcc_cost), delta, xsi_star))

        if abs(augment_amount) <= nfullytol and xsi_star <= nfullytol:
            # print('nmcc_cost is: {}, delta is: {}, and xsi_star is: {}'.format(
            #             float(nmcc_cost), delta, xsi_star))

            # if iters >= 10:
            #     if iters %10 == 0:
            if not consistent:
                here = iters
                consistent = True
            elif consistent and not consistent2:
                if here == iters - 1:
                    consistent2 = True
                else:
                    consistent = False
            elif consistent and consistent2:
                if here == iters - 2:
                    # print('xsi star is very small')
                    # print('nmcc_cost is: {}, delta is: {}, and xsi_star is: {}'.format(
                    #     float(nmcc_cost), delta, xsi_star))
                    break
                else:
                    consistent = False
                    consistent2 = False
        # if iters >= 20:
        #     if iters%20==0:
        #         varcost=G.var_cost()
        #         if abs(100*(np.array(np.sqrt(prev_var_cost))/np.sqrt(varcost) -1)) < difftol:
        #             if consistent3 == False:
        #                 here2 = iters
        #                 consistent3 = True
        #             elif consistent3 == True and here2 == iters-20:
        #                 break
        #             else:
        #                 consistent3 = False
        #         prev_var_cost = varcost
        #         #
        #
        #     else:
        #         if abs(augment_amount) <= tol:
        #             # if consistent == False:
        #                 # here = iters
        #                 # consistent = True
        #             # elif consistent == True and here == iters-1:
        #             break
        #             # else:
        #                 # consistent = False
        #
        #
        # if lamsearch:
        #
        #
        #     # if abs(augment_amount) <= tol:
        #     #     if consistent == False:
        #     #         here = iters
        #     #         consistent = True
        #     #     elif consistent == True and here == iters-1:
        #     #         break
        #     #     else:
        #     #         consistent = False
        #             #
        #     # if float(nmcc_cost) > -10:
        #     #     break
        #     if nr_run:
        #         if abs(augment_amount) <= tol:
        #             break
        #
        #     if iters%5==0:
        #         varcost=G.var_cost()
        #         if abs(100*(np.array(np.sqrt(prev_var_cost))/np.sqrt(varcost) -1)) < 1e-1:
        #
        #         # if abs(np.sqrt(prev_var_cost) - np.sqrt(varcost)) < var_cost_noful_tol:
        #             if consistent == False:
        #                 here = iters
        #                 consistent = True
        #             elif consistent == True and here == iters-5:
        #                 break
        #             else:
        #                 consistent = False
        #         prev_var_cost = varcost
        #
        # # else:
        # #     break

        R.augment_flow(augment_amount)
        G.adjust_flow(nmcc, augment_amount)

    # print('iters: ', iters)
    varcost = G.var_cost()
    nmcc_time = nmcc_tot_ms * 0.001
    return discount, nmcc_time, scc_time, varcost


def solve_mcf_sa_lin(G, R, varcost=None, discount=0, nmcc_tot_ms=0, scc_time=0, difftol=1e-1, tol=1e-2, fully=False, nfullytol=5e-4, lamsearch=False, var_cost_noful_tol=5e0, var_cost_ful=5e0, nr_run=False, vartop=1e1, mutop=1e2):
    prev_augment_amount = -1.0
    nmcc_exists = True
    iters = 0
    consistent = False
    consistent2 = False
    consistent3 = False
    here = -10
    here2 = -10
    prev_var_cost = 10

    while nmcc_exists:
        iters += 1
        st_write = time.time()
        
        strconnected = True

        # scc_time_start = time.time()

        # can=nx.is_strongly_connected(R.nxg)
        # if not can:
        # strconnected = False
        # print('not connected')
        # try:
        #     scc = max(nx.strongly_connected_components(G), key=len)
        #
        # except:
        #     pass
        # scc_time += time.time() - scc_time_start

        if not strconnected:
            scc_time_start = time.time()
            scc = [c for c in sorted(
                nx.strongly_connected_components(R.nxg), key=len, reverse=True)]
            scc_time += time.time() - scc_time_start
            scc_count = len(scc)
            maxim = -np.inf
            for i in range(scc_count):
                lenlist = len(scc[i])
                if lenlist > maxim:
                    maxim = max(maxim, lenlist)
                    max_i = i
            nodelist = scc[max_i]
            cur_g = copy.deepcopy(R.nxg)
            allnodes = np.arange(len(R.nxg.nodes)) + 1
            removals = list(set(allnodes) - set(nodelist))
            cur_g.remove_nodes_from(removals)

            output_graph(filename=graph_loc, g=cur_g, G=G)
            solve_nmcc()
            nmcc = read_nmcc(filename=nmcc_loc)     

        else:
            output_graph(filename=graph_loc, g=R, G=G)
            solve_nmcc()
            nmcc = read_nmcc(filename=nmcc_loc)            


        if len(nmcc) < 3:
            print('len nmcc is 2')
            break

        if float(nmcc_cost) >= 0:
            print('nmcc cost bigger than 0')
            break

        R.set_nmcc_lemon(nmcc)
        G.set_nmcc_lemon(nmcc)
        nmcc_tot_ms += nmcc_time_ms
        discount += time.time() - st_write

        delta = R.find_delta()
        # xsi_star = R.find_xsi_star(G)

        augment_amount = delta
        # if not lamsearch:
        #     if fully:
        if abs(augment_amount) <= nfullytol:

            if not consistent:
                here = iters
                consistent = True
            elif consistent and not consistent2:

                if here == iters - 1:
                    consistent2 = True
                else:
                    consistent = False
            elif consistent and consistent2:

                if here == iters - 2:
                    print('xsi star is very small')
                    print('nmcc_cost is: {}, delta is: {}, and xsi_star is: {}'.format(
                        float(nmcc_cost), delta, xsi_star))

                    break
                else:
                    consistent = False
                    consistent2 = False

                # if iters >= 20:
                #     if iters%20==0:
                #         varcost=G.var_cost()
                #         print(varcost)
                #         print(xsi_star)
                #         print('diff: ', abs(100*(np.array(np.sqrt(prev_var_cost))/np.sqrt(varcost) -1)))
                #         if abs(100*(np.array(np.sqrt(prev_var_cost))/np.sqrt(varcost) -1)) < difftol:
                #             if consistent3 == False:
                #                 here2 = iters
                #                 consistent3 = True
                #             elif consistent3 == True and here2 == iters-20:
                #                 break
                #             else:
                #                 consistent3 = False
                #         prev_var_cost = varcost
        #         #
        #
        #     else:
        #         if abs(augment_amount) <= tol:
        #             # if consistent == False:
        #                 # here = iters
        #                 # consistent = True
        #             # elif consistent == True and here == iters-1:
        #             break
        #             # else:
        #                 # consistent = False
        #
        #
        # if lamsearch:
        #
        #
        #     # if abs(augment_amount) <= tol:
        #     #     if consistent == False:
        #     #         here = iters
        #     #         consistent = True
        #     #     elif consistent == True and here == iters-1:
        #     #         break
        #     #     else:
        #     #         consistent = False
        #             #
        #     # if float(nmcc_cost) > -10:
        #     #     break
        #     if nr_run:
        #         if abs(augment_amount) <= tol:
        #             break
        #
        #     if iters%5==0:
        #         varcost=G.var_cost()
        #         if abs(100*(np.array(np.sqrt(prev_var_cost))/np.sqrt(varcost) -1)) < 1e-1:
        #
        #         # if abs(np.sqrt(prev_var_cost) - np.sqrt(varcost)) < var_cost_noful_tol:
        #             if consistent == False:
        #                 here = iters
        #                 consistent = True
        #             elif consistent == True and here == iters-5:
        #                 break
        #             else:
        #                 consistent = False
        #         prev_var_cost = varcost
        #
        # # else:
        # #     break

        R.augment_flow(augment_amount)
        G.adjust_flow(nmcc, augment_amount)

    print('iters: ', iters)
    varcost = G.var_cost()
    nmcc_time = nmcc_tot_ms * 0.001
    return discount, nmcc_time, scc_time, varcost


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
lambar = 0.7
G = MCF_DiGraph(lambar)
G.nxg = nx.DiGraph()
G = load_data(networkFileName='networks/goto_8_16a.min', G=G, generator='goto')
# solver_obj, solver_elapsed, solver_overhead, solver_x, solver_mu, solver_var = mosek_solve(G)
#solver_obj = 212132961996.7943
#solver_time = 66.10268497467041
st = time.time()
G.find_feasible_flow()
print(time.time() - st)
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









# num_nodes = [250]
# for n in num_nodes:
#     experiment_name = 'varying_lams'
#     num_nodes = [n]
#     lams = np.linspace(0, 5, num=50)
#     # lams = [0.1, 1, 10, 100]
#     mus = [50]
#     variances = [50]
#     d_tops = [1]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
# experiment_name, multi_arcs=False, test4dense=False, fix_bs=False,
# start_range=1, repeat=2)

# pdb.set_trace()
# num_nodes = [250]
# for n in num_nodes:
#     experiment_name = 'varying_lams'
#     num_nodes = [n]
#     lams = np.linspace(0, 5, num=50)
#     # lams = [0.1, 1, 10, 100]
#     mus = [50]
#     variances = [50]
#     d_tops = [1]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
# experiment_name, multi_arcs=False, test4dense=True, fix_bs=False,
# start_range=1)

# pdb.set_trace()


# num_nodes = [12]
# for n in num_nodes:
#     experiment_name = 'performance'
#     num_nodes = [n]
#     # lams = np.linspace(0, 5, num=50)
#     lams = [0.1, 5]
#     mus = [50]
#     variances = [50]
#     d_tops = [10]
#     pdb.set_trace()
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
#                    experiment_name, multi_arcs=False, test4dense=False, fix_bs=False, repeat=1)

# pdb.set_trace()


# num_nodes = [10000]
# for n in num_nodes:
#     experiment_name = 'performance'
#     num_nodes = [n]
#     lams = [0.7]
#     mus = [50]
#     d_tops = [10]
#     variances = [20]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
#                    experiment_name, multi_arcs=False, fix_bs=False, fixmsk=False, repeat=10)


# num_nodes = [25000]
# for n in num_nodes:
#     experiment_name = 'network_density'
#     num_nodes = [n]
#     lams = [0.7]
#     mus = [50]
#     d_tops = [1]
#     variances = [20]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
#                    experiment_name, multi_arcs=True, fix_bs=False, repeat=10)

# num_nodes = [5000, 50000]
# for n in num_nodes:
#     experiment_name = 'performance'
#     num_nodes = [n]
#     lams = [0.7]
#     mus = [50]
#     d_tops = [1]
#     variances = [20]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
# experiment_name, multi_arcs=False, fix_bs=False, fixmsk=False,
# repeat=10)


# num_nodes = [250]
# for n in num_nodes:
#     experiment_name = 'varying_lams'
#     num_nodes = [n]
#     lams = np.linspace(0, 5, num=50)
#     # lams = [0.1, 1, 10, 100]
#     print(lams)
#     mus = [50]
#     variances = [50]
#     d_tops = [1]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
# experiment_name, multi_arcs=False, test4dense=False, fix_bs=False,
# repeat=1)


# num_nodes = [250]
# for n in num_nodes:
#     experiment_name = 'varying_lams'
#     num_nodes = [n]
#     lams = np.linspace(0, 5, num=50)
#     # lams = [0.1, 1, 10, 100]
#     print(lams)
#     mus = [50]
#     variances = [50]
#     d_tops = [1]
#     run_experiment(num_nodes, lams, mus, variances, d_tops,
# experiment_name, multi_arcs=False, test4dense=True, fix_bs=False,
# repeat=1)
