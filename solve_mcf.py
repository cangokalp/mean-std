from utils import *
import pandas
import random
import itertools
import copy
import cProfile
import pstats
from prettytable import PrettyTable
import caffeine
from cplex.callbacks import BarrierCallback


class TimeLimit_LoggingCallback(BarrierCallback):

    def __call__(self):
        self.duals.append(self.get_dual_objective_value())
        self.primals.append(self.get_objective_value())
        self.is_feasible.append(self.is_primal_feasible())
        self.primal_infeas.append(self.get_primal_infeasibility())
        timeused = self.get_time() - self.starttime
        self.iter_times.append(timeused)

        dual = self.get_dual_objective_value()
        primal = self.get_objective_value()


class TimeLimitCallback(BarrierCallback):

    def __call__(self):

        dual = self.get_dual_objective_value()
        primal = self.get_objective_value()
        rel_gap = abs(dual - primal) / abs(min(primal, dual))
        if rel_gap < self.acceptablegap:
            self.aborted = 1
            self.abort()


def cvxpy_solve(G, cvxpy=False, fusion=False, run_name='.', acceptablegap=1e-8):
    start = time.time()

    if fusion:
        G.A_msk = mf.Matrix.sparse(G.n, G.m, G.rows, G.cols, G.values)

        with mf.Model() as model:

            x = model.variable("x", G.m, mf.Domain.greaterThan(0.0))
            x_n = model.variable("x_n", G.m)
            decoy = model.variable("decoy", 1)
            decoy_2 = mf.Var.vstack(decoy, x_n)
            model.objective("myobj",
                            mf.ObjectiveSense.Minimize,
                            mf.Expr.add(mf.Expr.dot(G.mu, x), mf.Expr.mul(lambar, decoy)))
            model.constraint(mf.Expr.sub(x_n, mf.Expr.mulElm(np.sqrt(G.var), x)),
                             mf.Domain.equalsTo(0.0))
            model.constraint(x, mf.Domain.lessThan(G.cap))
            model.constraint(decoy_2, mf.Domain.inQCone())
            model.constraint(mf.Expr.sub(mf.Expr.mul(
                G.A_msk, x), G.b), mf.Domain.equalsTo(0.0))

            model.setSolverParam("intpntCoTolRelGap", 1.0e-8)
            model.setSolverParam("intpntCoTolPfeas", 1.0e-8)
            model.setSolverParam("intpntCoTolDfeas", 1.0e-8)
            model.setSolverParam("intpntCoTolMuRed", 1.0e-8)
            model.setSolverParam("intpntMaxIterations", 100000)

            import logging
            log_path = os.path.join(EXPERIMENT_PATH, 'mosek')
            log_path = os.path.join(log_path, run_name)

            logging = Logger(run_dir=log_path)
            model.setLogHandler(logging)
            model.setSolverParam("logIntpnt", 1000)
            model.setSolverParam("log", 1000)
            model.setSolverParam("logFile", 1000)
            model.solve()

            obj = model.primalObjValue()
            x = x.level()
            # mean_cost = G.mu.dot(x)
            # var_cost = G.var.dot(x)
            elapsed = time.time() - start
            print(obj, elapsed)
            pdb.set_trace()
            return obj, elapsed, x, None

    if cvxpy:
        x = cp.Variable(G.m)
        constraints = [0 <= x, x <= G.cap, G.A@x == G.b]
        objective = cp.Minimize(
            G.mu.T * x + G.lambar * cp.norm(cp.multiply(np.sqrt(G.var), x), 2))
        prob = cp.Problem(objective, constraints)

        # result = prob.solve(solver=solver, verbose=True)  # gurobi mosek
        # compare

        result = prob.solve(solver='MOSEK', verbose=False, mosek_params={'MSK_IPAR_LOG_INTPNT': 2, 'MSK_DPAR_INTPNT_QO_TOL_REL_GAP': 1e-14, 'MSK_DPAR_INTPNT_CO_TOL_MU_RED': 1e-12,
                                                                         'MSK_DPAR_INTPNT_CO_TOL_PFEAS': 1e-12, 'MSK_DPAR_INTPNT_CO_TOL_REL_GAP': 1e-12, 'MSK_DPAR_INTPNT_CO_TOL_INFEAS': 1e-12, 'MSK_IPAR_INTPNT_MAX_ITERATIONS': 10000})

        cvx_soln = x.value
        cvx_obj = objective.value

        # print(prob.solver_stats.solve_time, prob.solver_stats.setup_time)
        cvx_elapsed = prob.solver_stats.solve_time
        elapsed = time.time() - start

        return cvx_obj, elapsed, cvx_soln, None

    else:

        prob = cplex.Cplex()
        cplex_cb = prob.register_callback(TimeLimit_LoggingCallback)
        cplex_cb.acceptablegap = acceptablegap
        cplex_cb.aborted = 0
        cplex_cb.starttime = prob.get_time()

        cplex_cb.duals = []
        cplex_cb.primals = []
        cplex_cb.is_feasible = []
        cplex_cb.primal_infeas = []
        cplex_cb.iter_times = []
        cplex_cb.feasible_dual = []
        cplex_cb.feasible_primal = []

        x_names = ['x' + str(i) for i in range(G.m)]
        lin_obj = G.mu
        prob.variables.add(obj=lin_obj, lb=np.zeros(G.m),
                           ub=G.cap, names=x_names)

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
        prob.parameters.barrier.qcpconvergetol.set(1e-12)
        # prob.parameters.barrier.algorithm.set(2)
        # prob.parameters.barrier.limits.corrections.set(100)

        prob.solve()
        sol = prob.solution

        obj = prob.solution.get_objective_value()
        soln = np.array(prob.solution.get_values())

        m = prob.solution.quality_metric
        max_x, max_infeas = prob.solution.get_float_quality(
            [m.max_x, m.max_primal_infeasibility])

        elapsed = time.time() - start
        status = sol.status[sol.get_status()]
        if status != 'optimal' and status != 'abort_user' and status != 'num_best':
            pdb.set_trace()

        cplex_cb.primals.append(obj)
        cplex_cb.iter_times.append(elapsed)
        feasible = prob.solution.is_primal_feasible()
        return status, feasible, obj, elapsed, soln, max_infeas, cplex_cb.duals, cplex_cb.primals, cplex_cb.is_feasible, cplex_cb.primal_infeas, cplex_cb.iter_times, cplex_cb.feasible_dual, cplex_cb.feasible_primal


def cvxpy_solve_additive(G, lam, prob=None, warm_start=False, lp=False, solver='MOSEK', bound_lam=False, cvxpy=False, acceptablegap=1e-6, callback=True):

    start = time.time()

    if cvxpy:
        start = time.time()
        x = cp.Variable(G.m)
        constraints = [0 <= x, x <= G.cap, G.A@x == G.b]
        if lp:
            objective = cp.Minimize(G.mu * x)

        elif bound_lam:
            objective = cp.Minimize(G.var * x**2)

        else:
            objective = cp.Minimize(G.mu.T * x + lam * G.var * x**2)

        prob = cp.Problem(objective, constraints)

        result = prob.solve(solver='MOSEK', verbose=False)
        cvx_soln = x.value
        cvx_obj = objective.value

        cvx_elapsed = prob.solver_stats.solve_time
        elapsed = time.time() - start

        return cvx_obj, elapsed, cvx_soln, None, None

    else:

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
            prob.linear_constraints.set_coefficients(
                zip(G.rows, G.cols, G.values))

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

        if callback:
            timelim_cb = prob.register_callback(TimeLimitCallback)
            timelim_cb.acceptablegap = acceptablegap
            timelim_cb.aborted = 0

        prob.solve()

        obj = prob.solution.get_objective_value()
        soln = np.array(prob.solution.get_values())

        elapsed = time.time() - start

        m = prob.solution.quality_metric
        max_x, max_infeas = prob.solution.get_float_quality(
            [m.max_x, m.max_primal_infeasibility])

        return obj, elapsed, soln, prob, max_infeas


def cvxpy_solve_xi(G, soln, lam, prob=None, warm_start=False, lp=False, solver='MOSEK', iters=0, cvxpy=False, acceptablegap=1e-3, callback=True):

    start = time.time()
    trial = False
    cvxpy = False
    if cvxpy:

        # xi = cp.Variable(G.m)
        # x_ = soln
        # constraints = [G.A@xi == 0, -x_ <= xi, xi <= G.cap]

        # diff = x_ - G.cap
        # x_zero = np.argwhere(abs(x_) < 1e-6).ravel().astype(int)
        # x_nonzero = np.argwhere(x_ > 1e-6).ravel().astype(int)
        # x_u = np.argwhere(abs(diff) < 1e-6).ravel().astype(int)
        # x_btw = np.array(list(set(x_nonzero).difference(set(x_u)))).astype(int)
        # diff_ = G.cap[x_btw] - x_[x_btw]

        # # if len(x_zero) > 0:
        # # constraints.append(xi[x_zero] == np.zeros(len(x_zero)))

        # if len(x_u) > 0:
        #     constraints.append(xi[x_u] <= np.zeros(len(x_u)))
        #     constraints.append(xi[x_u] >= -x_[x_u])

        # if len(x_btw) > 0:
        #     constraints.append(xi[x_btw] <= diff_)
        #     constraints.append(xi[x_btw] >= -x_[x_btw])

        # objective = cp.Minimize(
        #     G.lam * G.var.T * cp.square(xi) + 2 * np.multiply(G.var, x_).T * xi)
        # prob = cp.Problem(objective, constraints)

        # result = prob.solve(solver='MOSEK', verbose=False)
        # cvx_soln = xi.value
        # cvx_obj = objective.value

        # cvx_elapsed = prob.solver_stats.solve_time
        # elapsed = time.time() - start
        # return elapsed, cvx_soln, None, None
        import cvxpy as cp
        x_zero = np.argwhere(abs(soln) < 1e-4).ravel().astype(int)
        x_nonzero = np.argwhere(soln > 1e-4).ravel().astype(int)
        diff = abs(np.array(soln) - np.array(G.cap))/np.minimum(np.array(G.cap), np.array(soln)) * 100
        x_u = np.argwhere(abs(diff) < 1e-2).ravel().astype(int)
        x_btw = np.array(list(set(x_nonzero).difference(set(x_u)))).astype(int)


        xi = cp.Variable(len(x_btw))
        x_ = soln

        Gcols = []
        Grows = []
        Gvals = []
        construct_dict = {}
        construct_count = 0
        nodeset = []

        for ell in x_btw:
            (u, v) = G.arc_dict[ell]
            
            Grows.append(u - 1)
            Gcols.append(construct_count)
            Gvals.append(1.0)
            Grows.append(v - 1)
            Gcols.append(construct_count)
            Gvals.append(-1.0)
            construct_count += 1
            nodeset.append(u-1)
            nodeset.append(v-1)


        nodeset = list(set(nodeset))

        nGrows = [nodeset.index(i) for i in Grows]
        nGA = scipy.sparse.csc_matrix((Gvals, (nGrows, Gcols)))       

        constraints = [nGA@xi == 0, -np.inf <= xi, xi <= np.inf]

        objective = cp.Minimize(
            G.lam * G.var[x_btw].T * cp.square(xi) + 2 * np.multiply(G.var[x_btw], x_[x_btw]).T * xi)
        prob = cp.Problem(objective, constraints)

        result = prob.solve(solver='CPLEX', verbose=False)
        cvx_soln = xi.value
        cvx_obj = objective.value

        cvx_elapsed = prob.solver_stats.solve_time
        elapsed = time.time() - start
        
        construct_soln = x_
        construct_soln[x_btw] = cvx_soln
        construct_soln[~x_btw] = 0

        return elapsed, construct_soln, None, None
    else:
        if trial:
            warm_start = False
            diff = soln - G.cap
            x_ = soln

            x_zero = np.argwhere(abs(soln) < 1e-4).ravel().astype(int)
            x_nonzero = np.argwhere(soln > 1e-4).ravel().astype(int)
            diff = abs(np.array(soln) - np.array(G.cap))/np.minimum(np.array(G.cap), np.array(soln)) * 100
            x_u = np.argwhere(abs(diff) < 1e-2).ravel().astype(int)
            x_btw = np.array(list(set(x_nonzero).difference(set(x_u)))).astype(int)

            diff_ = G.cap[x_btw] - x_[x_btw]

            cols = np.r_[x_btw, x_btw]

            rows = np.arange(len(cols))
            if not warm_start:

                prob = cplex.Cplex()
                if lp:
                    obj = np.ones(len(x_btw))
                    prob.variables.add(obj=obj, lb=-x_[x_btw], ub=diff_)
                else:
                    # prob.variables.add(lb=-x_[x_btw], ub=diff_)
                    prob.variables.add(lb=-np.inf * np.ones(len(x_btw)), ub=np.inf * np.ones(len(x_btw)))

                    lin_coeffs = 2 * G.var[x_btw] * x_[x_btw]
                    quad_coeffs = lam * G.var[x_btw] * 2
                    prob.objective.set_linear(
                        [(int(i), j) for i, j in zip(np.arange(len(x_btw)), lin_coeffs)])
                    prob.objective.set_quadratic(quad_coeffs)

                Gcols = []
                Grows = []
                Gvals = []
                construct_dict = {}
                construct_count = 0
                nodeset = []

                for ell in x_btw:
                    (u, v) = G.arc_dict[ell]
                    
                    Grows.append(u - 1)
                    Gcols.append(construct_count)
                    Gvals.append(1.0)
                    Grows.append(v - 1)
                    Gcols.append(construct_count)
                    Gvals.append(-1.0)
                    construct_count += 1
                    nodeset.append(u-1)
                    nodeset.append(v-1)


                nodeset = list(set(nodeset))

                nGrows = [nodeset.index(i) for i in Grows]

                prob.linear_constraints.add(rhs=np.zeros(
                    len(nodeset)), senses='E' * len(nodeset))
                prob.linear_constraints.set_coefficients(
                    zip(nGrows, Gcols, Gvals))

                prob.objective.set_sense(prob.objective.sense.minimize)
                prob.set_log_stream(None)
                prob.set_error_stream(None)
                prob.set_warning_stream(None)
                prob.set_results_stream(None)
                prob.parameters.barrier.display.set(0)

        else:
            diff = soln - G.cap
            x_ = soln
            x_zero = np.argwhere(abs(soln) < 1e-5).ravel().astype(int)
            x_nonzero = np.argwhere(soln > 1e-5).ravel().astype(int)
            x_u = np.argwhere(abs(diff) < 1e-5).ravel().astype(int)
            x_btw = np.array(
                list(set(x_nonzero).difference(set(x_u)))).astype(int)

            diff_ = G.cap[x_btw] - x_[x_btw]

            b = np.r_[np.zeros(len(x_zero)), np.zeros(
                len(x_u))]

            # cols = np.r_[x_zero, x_u]
            # senses = np.r_[['E'] * len(x_zero), ['E'] * len(x_u)]

            # rows = np.arange(len(cols))
            # values = np.ones(len(rows))

            if not warm_start:

                prob = cplex.Cplex()
                if lp:
                    obj = np.ones(G.m)
                    prob.variables.add(obj=obj, lb=-x_ *
                                       np.ones(G.m), ub=G.cap)
                else:
                    # prob.variables.add(lb=-np.inf * np.ones(G.m), ub=G.cap)
                    prob.variables.add(lb=-np.inf * np.ones(G.m), ub=np.inf * np.ones(G.m))

                    lin_coeffs = 2 * G.var * x_
                    quad_coeffs = lam * G.var * 2
                    prob.objective.set_linear(
                        [(int(i), j) for i, j in zip(np.arange(G.m), lin_coeffs)])
                    if lam!=0:
                        prob.objective.set_quadratic(quad_coeffs)

                prob.linear_constraints.add(
                    rhs=np.zeros(G.n), senses='E' * G.n)
                prob.linear_constraints.set_coefficients(
                    zip(G.rows, G.cols, G.values))

                prob.variables.set_lower_bounds(zip([int(i) for i in x_btw], np.zeros(len(x_btw))))
                prob.variables.set_upper_bounds(zip([int(i) for i in x_btw], np.zeros(len(x_btw))))

                # prob.linear_constraints.add(rhs=b, senses=senses)
                # prob.linear_constraints.set_coefficients(
                #     zip(rows.astype(int), cols.astype(int), values))
                # prob.linear_constraints.set_coefficients(
                    # zip([int(i) for i in rows], [int(i) for i in cols], values))

                prob.objective.set_sense(prob.objective.sense.minimize)
                prob.set_log_stream(None)
                prob.set_error_stream(None)
                prob.set_warning_stream(None)
                prob.set_results_stream(None)
                prob.parameters.barrier.display.set(0)
            else:

                prob.variables.set_lower_bounds(zip([int(i) for i in range(G.m)], np.ones(G.m)*(-np.inf)))
                prob.variables.set_upper_bounds(zip([int(i) for i in range(G.m)], np.ones(G.m)*(np.inf)))

                prob.variables.set_lower_bounds(zip([int(i) for i in x_btw], np.zeros(len(x_btw))))
                prob.variables.set_upper_bounds(zip([int(i) for i in x_btw], np.zeros(len(x_btw))))


                # prob.linear_constraints.delete(
                #     G.n, prob.linear_constraints.get_num() - 1)
                # prob.linear_constraints.add(rhs=b, names=names, senses=senses)
                # prob.linear_constraints.set_coefficients(
                #     zip(rows.astype(int), cols.astype(int), values))

                if not lp:
                    lin_coeffs = 2 * G.var * x_
                    quad_coeffs = lam * G.var * 2
                    prob.objective.set_linear(
                        [(int(i), j) for i, j in zip(np.arange(G.m), lin_coeffs)])
                    prob.objective.set_quadratic(quad_coeffs)

        if callback:
            timelim_cb = prob.register_callback(TimeLimitCallback)
            timelim_cb.acceptablegap = acceptablegap
            timelim_cb.aborted = 0
        # pdb.set_trace()
        prob.solve()

        obj = prob.solution.get_objective_value()
        soln = np.array(prob.solution.get_values())
        construct_soln = soln
        # construct_soln = x_
        # construct_soln[x_btw] = soln
        # construct_soln[~x_btw] = 0

        elapsed = time.time() - start

        m = prob.solution.quality_metric
        max_x, max_infeas = prob.solution.get_float_quality(
            [m.max_x, m.max_primal_infeasibility])

        return elapsed, construct_soln, prob, max_infeas


def bs_cvxpy(G, low=0, high=1000, prob=None, lp=False, solver_weird=False, cvxpy=False, solver_obj=None):
    print('bsc')

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
    acceptablegap = 1e-8

    while not found:
        iters += 1
        mid_prev = mid
        mid = (high + low) / 2.0

        if iters > 1 and warm_start == False:
            warm_start = True

        lams.append(mid)
        print(mid)
        obj, elapsed, x, prob, max_infeas = cvxpy_solve_additive(
            G, mid, prob=prob, warm_start=warm_start, lp=lp, cvxpy=cvxpy, acceptablegap=acceptablegap, callback=False)

        infeas.append(max_infeas)
        var_cost = np.multiply(G.var, x).dot(x)
        obj = G.mu.dot(x) + G.lambar * (np.sqrt(var_cost))

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
        if iters > 30 or 100*abs(obj-solver_obj)/min(obj,solver_obj) < 1e-3:
            break

    elapsed = time.time() - start
    if cvxpy:
        return obj, elapsed, x, iter_objs, iter_elapsed, None, lams, fs

    else:
        infeas = np.array(infeas)
        return obj, elapsed, x, iter_objs, iter_elapsed, infeas.mean(), lams, fs


def nr_cvxpy(G, low=0, high=1000, prob=None, lp=False, lam_init=None, solver_weird=False, cvxpy=False, solver_obj=None):
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
    acceptablegap = 1e-8
    callback = False
    while not found:
        iters += 1

        if iters > 1:
            warm_start = True
            warm_start_xi = True
        print(lam)
        lams.append(lam)

        obj, elapsed, x, prob, max_infeas = cvxpy_solve_additive(
            G, lam, prob=prob, warm_start=warm_start, lp=lp, cvxpy=cvxpy, acceptablegap=acceptablegap, callback=callback)
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
            G, x, lam, prob=prob_xi, warm_start=warm_start_xi, lp=lp, iters=iters, cvxpy=cvxpy, acceptablegap=1e-8, callback=callback)
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

        if iters > 9 or 100*abs(obj-solver_obj)/min(obj,solver_obj) < 1e-3:
            break

    elapsed = time.time() - start

    if cvxpy:
        return obj, elapsed, x, iter_objs, iter_elapsed, iter_xi_times, iter_var_times, mean, var,  None, lams, fs

    else:
        infeas = np.array(infeas)
        return obj, elapsed, x, iter_objs, iter_elapsed, iter_xi_times, iter_var_times, mean, var,  infeas.mean(), lams, fs


def hybrid(G, low=0, high=1000, prob=None, lp=False, solver_weird=False, cvxpy=False, solver_obj=None):

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


def get_networks(experiment, tails, exponents, types, test=False, fam='netgen', lams=[]):

    networks = []
    fam_base = fam

    extension = 'pickle'
    for tail in tails:
        for exponent in exponents:
            if experiment == 'varying_lambar':
                for atype in types:
                    for lambar in lams:
                        lam_dir = str(lambar).replace('.', '_')
                        net_name = fam_base + atype + exponent + tail
                        cur_run = experiment + '/' + fam_base + \
                            '/' + net_name[:net_name.find('.')] + '_' + lam_dir
                        cur_run = os.path.join(
                            EXPERIMENT_PATH, cur_run + "." + extension)
                        if not test:
                            if not os.path.isfile(cur_run):
                                if not net_name in networks:
                                    networks.append(net_name)
                        else:
                            if not net_name in networks:
                                networks.append(net_name)

            else:

                for atype in types:
                    net_name = fam_base + atype + exponent + tail
                    cur_run = experiment + '/' + fam_base + \
                        '/' + net_name[:net_name.find('.')]
                    cur_run = os.path.join(
                        EXPERIMENT_PATH, cur_run + "." + extension)
                    # if not test:
                    #     if not os.path.isfile(cur_run):
                    #         networks.append(net_name)
                    # else:
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


def graph_family_experiment(networks, lambar, record=True, cvxpy=False, fusion=False, which='cp'):
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

        cvxpy_alg = cvxpy


        if network.find('15') >= 0:

            cp_saved = experiment + '/' + generator + '/' + \
                network[:network.find('.')] + "_cp_results"

            if not os.path.isfile(os.path.join(EXPERIMENT_PATH, cp_saved) + '.pickle'):

                status, feasible, solver_obj, solver_elapsed, solver_soln, solver_infeas, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal = cvxpy_solve(
                    G, cvxpy=cvxpy, fusion=fusion, run_name=run_name)

                keys = ['solver_status', 'solver_feasible', 'solver_elapsed', 'solver_obj', 'solver_duals', 'solver_primals',
                    'solver_is_feasible', 'solver_primal_infeas', 'solver_iter_times', 'solver_feasible_dual', 'solver_feasible_primal', 'solver_infeas']
                values = [status, feasible, solver_elapsed, solver_obj, solver_duals, solver_primals, solver_is_feasible,
                    solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal, solver_infeas]

                run_dict = dict(zip(keys, values))

                os.makedirs(cp_saved, exist_ok=True)

                save_run(cp_saved, run_dict)
 
            if which == 'nr':

                nr_saved = experiment + '/' + generator + '/' + \
                    network[:network.find('.')] + "_nr_results"

                if not os.path.isfile(os.path.join(EXPERIMENT_PATH, nr_saved) + '.pickle'):
                    filename = cp_saved
                    data_dic = load_run(filename)
                    solver_obj = data_dic['solver_obj']

                    _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(
                        G, lam=0, lp=True, cvxpy=cvxpy_alg)
                    var_cost = np.multiply(G.var, soln).dot(soln)
                    lam_low = G.lambar / (2.0 * np.sqrt(var_cost))

                    nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs = nr_cvxpy(
                        G, lp=False, low=lam_low, high=G.lambar / 2.0, lam_init=lam_low, cvxpy=cvxpy_alg, solver_obj=solver_obj)

                    keys = ['nr_elapsed', 'nr_iter_elapsed', 'nr_iter_objs', 'nr_obj',
                                'nr_xi_times', 'nr_var_times', 'mean', 'var', 'nr_avg_infeas', 'nr_lams', 'nr_fs', 'elapsed_lower_bound']
                    values = [nr_elapsed, nr_iter_elapsed, nr_iter_objs, nr_obj, nr_xi_times,
                        nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs, lb_elapsed]

                    run_dict = dict(zip(keys, values))

                    os.makedirs(nr_saved, exist_ok=True)

                    save_run(nr_saved, run_dict)
                continue

            if which == 'bs':
                bs_saved = experiment + '/' + generator + '/' + \
                    network[:network.find('.')] + "_bs_results"

                if not os.path.isfile(os.path.join(EXPERIMENT_PATH, bs_saved) + '.pickle'):
                    filename = cp_saved
                    data_dic = load_run(filename)
                    solver_obj = data_dic['solver_obj']


                    _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(
                        G, lam=0, lp=True, cvxpy=cvxpy_alg)
                    var_cost = np.multiply(G.var, soln).dot(soln)
                    lam_low = G.lambar / (2.0 * np.sqrt(var_cost))


                    _, ub_elapsed, soln, _, _ = cvxpy_solve_additive(
                        G, lam=0, lp=False, bound_lam=True, cvxpy=cvxpy_alg, acceptablegap=1e-2)
                    var_cost = np.multiply(G.var, soln).dot(soln)
                    lam_high = G.lambar / (2.0 * np.sqrt(var_cost))

                    bs_obj, bs_elapsed, _, bs_iter_objs, bs_iter_elapsed, bs_avg_infeas, bs_lams, bs_fs = bs_cvxpy(
                    G, lp=False, low=lam_low, high=lam_high, cvxpy=cvxpy_alg,solver_obj=solver_obj)

                    keys = ['bs_obj', 'bs_iter_elapsed', 'bs_iter_objs', 'bs_elapsed', 'bs_avg_infeas',
                        'bs_lams', 'bs_fs', 'elapsed_lower_bound', 'elapsed_upper_bound', ]
                    values = [bs_obj, bs_iter_elapsed, bs_iter_objs, bs_elapsed,
                        bs_avg_infeas, bs_lams, bs_fs, lb_elapsed, ub_elapsed]

                    run_dict = dict(zip(keys, values))

                    os.makedirs(bs_saved, exist_ok=True)

                    save_run(bs_saved, run_dict)
                continue


        # continue
        if not os.path.isfile(os.path.join(EXPERIMENT_PATH, cplex_saved) + '.pickle'):
            status, feasible, solver_obj, solver_elapsed, solver_soln, solver_infeas, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal = cvxpy_solve(
                G, cvxpy=cvxpy, fusion=fusion, run_name=run_name)

        else:
            filename = cplex_saved
            data_dic = load_run(filename)

            status = data_dic['solver_status']
            feasible = data_dic['solver_feasible']
            solver_elapsed = data_dic['solver_elapsed']
            solver_obj = data_dic['solver_obj']
            solver_duals = data_dic['solver_duals']
            solver_primals = data_dic['solver_primals']
            solver_is_feasible =data_dic['solver_is_feasible']
            solver_primal_infeas = data_dic['solver_primal_infeas']
            solver_iter_times = data_dic['solver_iter_times']
            solver_feasible_dual = data_dic['solver_feasible_dual']
            solver_feasible_primal = data_dic['solver_feasible_primal']
            solver_infeas = data_dic['solver_infeas']

        print('Solver finished in {} seconds, with objective {}'.format(
            solver_elapsed, solver_obj))


        # _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(
        #     G, lam=0, lp=True, cvxpy=cvxpy_alg)
        # var_cost = np.multiply(G.var, soln).dot(soln)
        # lam_low = G.lambar / (2.0 * np.sqrt(var_cost))
        # print('lb_elapsed: ', lb_elapsed)
        lb_elapsed = 0
        lam_low = 1e-6
        print(lam_low)
# 

        if not os.path.isfile(os.path.join(EXPERIMENT_PATH, run_name) + '.pickle'):


            nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs = nr_cvxpy(
                    G, lp=False, low=lam_low, high=G.lambar / 2.0, lam_init=lam_low, cvxpy=cvxpy_alg, solver_obj=solver_obj)
            nr_elapsed += lb_elapsed

            can_nr_elapsed = np.array(nr_iter_elapsed) + lb_elapsed
            can_nr_gaps = abs(np.array(nr_iter_objs) - solver_obj)/np.minimum(np.array(nr_iter_objs), solver_obj)

            print(can_nr_elapsed)
            print(can_nr_gaps)

            print('NR finished in {} seconds, with objective {}'.format(
                nr_elapsed, nr_obj))
        else:
            filename = run_name
            data_dic = load_run(filename)

            nr_elapsed = data_dic['nr_elapsed']
            nr_iter_elapsed = data_dic['nr_iter_elapsed']
            nr_iter_objs = data_dic['nr_iter_objs']
            nr_obj = data_dic['nr_obj']
            nr_xi_times = data_dic['nr_xi_times']
            nr_var_times = data_dic['nr_var_times']
            mean = data_dic['mean']
            var = data_dic['var']
            nr_avg_infeas = data_dic['nr_avg_infeas']
            nr_lams = data_dic['nr_lams']
            nr_fs = data_dic['nr_fs']





        if not os.path.isfile(os.path.join(EXPERIMENT_PATH, run_name) + '.pickle'):

            # _, ub_elapsed, soln, _, _ = cvxpy_solve_additive(
            #     G, lam=0, lp=False, bound_lam=True, cvxpy=cvxpy_alg, acceptablegap=1e-2)
            # var_cost = np.multiply(G.var, soln).dot(soln)
            # lam_high = G.lambar / (2.0 * np.sqrt(var_cost))
            ub_elapsed = 0
            lam_high = G.lambar / 2.0
            print(lam_high)

            print('ub_elapsed: ', ub_elapsed)
            bs_obj, bs_elapsed, _, bs_iter_objs, bs_iter_elapsed, bs_avg_infeas, bs_lams, bs_fs = bs_cvxpy(
            G, lp=False, low=lam_low, high=lam_high, cvxpy=cvxpy_alg, solver_obj=solver_obj)


            can_bs_elapsed = np.array(bs_iter_elapsed) + lb_elapsed + ub_elapsed
            can_bs_gaps = abs(np.array(bs_iter_objs) - solver_obj)/np.minimum(np.array(bs_iter_objs), solver_obj)

            print(can_bs_elapsed)
            print(can_bs_gaps)



            bs_elapsed += lb_elapsed + ub_elapsed
            print('BSC finished in {} seconds, with objective {}'.format(
                bs_elapsed, bs_obj))

        else:

            filename = run_name
            data_dic = load_run(filename)

            lb_elapsed = data_dic['elapsed_lower_bound']
            ub_elapsed = data_dic['elapsed_upper_bound']
            bs_obj = data_dic['bs_obj']
            bs_iter_elapsed = data_dic['bs_iter_elapsed']
            bs_iter_objs = data_dic['bs_iter_objs']
            bs_elapsed = data_dic['bs_elapsed']
            bs_avg_infeas =data_dic['bs_avg_infeas']
            bs_lams = data_dic['bs_lams']
            bs_fs = data_dic['bs_fs']



        # _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(
        #     G, lam=0, lp=True, cvxpy=cvxpy_alg)
        # var_cost = np.multiply(G.var, soln).dot(soln)
        # lam_low = G.lambar / (2.0 * np.sqrt(var_cost))
        # print('lam_low: ', lam_low)





           




        if cvxpy:
            t = PrettyTable(['Method', 'Soln Time', 
                             '# Iters', 'Obj', 'Rel_Gap'])
            t.add_row(['CPLEX', round(solver_elapsed, 2), '-',
                       solver_obj, '-'])
            t.add_row(['NR', round(nr_elapsed, 2), len(nr_iter_elapsed),
                       nr_obj, abs(solver_obj - nr_obj) / nr_obj])
            t.add_row(['BSC', round(bs_elapsed, 2), len(bs_iter_elapsed), bs_obj, abs(
                solver_obj - bs_obj) / bs_obj])
        else:
            t = PrettyTable(['Method', 'Soln Time', 'Status', 'Feasible',
                             '# Iters', 'Obj', 'Rel_Gap', 'Infeas'])
            t.add_row(['CPLEX', round(solver_elapsed, 2), status, feasible, '-',
                       solver_obj, '-', round(solver_infeas, 3)])
            t.add_row(['NR', round(nr_elapsed, 2), '-', '-', len(nr_iter_elapsed),
                       nr_obj, abs(solver_obj - nr_obj) / nr_obj, round(nr_avg_infeas, 3)])
            t.add_row(['BSC', round(bs_elapsed, 2), '-', '-', len(bs_iter_elapsed), bs_obj, abs(
                solver_obj - bs_obj) / bs_obj, round(bs_avg_infeas, 3)])
        
        print(t)


        if record:

            # solver_iter_objs, solver_iter_elapsed, presolve_time, ordering_time = parse_cplex_log()

            keys = ['network_name', 'solver_status', 'solver_feasible', 'solver_elapsed', 'solver_obj', 'solver_duals', 'solver_primals', 'solver_is_feasible', 'solver_primal_infeas', 'solver_iter_times', 'solver_feasible_dual', 'solver_feasible_primal', 'solver_infeas',
                    'nr_elapsed', 'nr_iter_elapsed', 'nr_iter_objs', 'nr_obj',
                    'nr_xi_times', 'nr_var_times', 'mean', 'var', 'nr_avg_infeas', 'nr_lams', 'nr_fs', 'lambar', 'm', 'n',
                    'elapsed_lower_bound', 'elapsed_upper_bound',
                    'bs_obj', 'bs_iter_elapsed', 'bs_iter_objs', 'bs_elapsed', 'bs_avg_infeas', 'bs_lams', 'bs_fs']

            values = [network[:network.find('.')], status, feasible, solver_elapsed, solver_obj, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal, solver_infeas,
                      nr_elapsed, nr_iter_elapsed, nr_iter_objs, nr_obj, nr_xi_times,
                      nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs, G.lambar, G.m, G.n, lb_elapsed, ub_elapsed,
                      bs_obj, bs_iter_elapsed, bs_iter_objs, bs_elapsed, bs_avg_infeas, bs_lams, bs_fs]

            run_dict = dict(zip(keys, values))

            if not os.path.isfile(cplex_saved):
                os.makedirs(cplex_saved, exist_ok=True)
                if cvxpy:
                    save_run(cplex_saved, run_dict, prefix='mosek')
                else:
                    save_run(cplex_saved, run_dict)
            if cvxpy:
                save_run(run_name, run_dict, prefix='mosek')
            else:
                save_run(run_name, run_dict)

def varying_lambda_experiment(networks, lambars, record=True):
    experiment = 'varying_lambar'

    for network in networks:
        for lambar in lambars:
            G = MCF_DiGraph(lambar)
            G.nxg = nx.DiGraph()
            print(network, lambar)
            if network.find('goto') >= 0:
                generator = 'goto'
                filename = 'networks/goto/' + network

            else:
                generator = 'netgen'
                filename = 'networks/netgen/' + network

            G = load_data(networkFileName=filename, G=G, generator=generator)

            lam_dir = str(lambar).replace('.', '_')
            run_name = experiment + '/' + generator + \
                '/' + network[:network.find('.')] + '_' + lam_dir

            solver_weird = False

            status, feasible, solver_obj, solver_elapsed, solver_soln, solver_infeas, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal = cvxpy_solve(
                G, cvxpy=cvxpy, fusion=fusion, run_name=run_name)

            print('Solver finished in {} seconds, with objective {}'.format(
                solver_elapsed, solver_obj))

            if solver_infeas > 2:
                solver_weird = True

            # getting bounds:
            _, lb_elapsed, soln, _, _ = cvxpy_solve_additive(G, lam=0, lp=True)
            var_cost = np.multiply(G.var, soln).dot(soln)
            lam_low = G.lambar / (2.0 * np.sqrt(var_cost))
            lam_high = G.lambar / (2.0 * np.sqrt(var_cost))


            nr_obj, nr_elapsed, _, nr_iter_objs, nr_iter_elapsed, nr_xi_times, nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs = nr_cvxpy(
                G, lp=False, low=lam_low, high=lam_high, lam_init=lam_low, solver_obj=solver_obj)

            nr_elapsed += lb_elapsed

            print('NR finished in {} seconds, with objective {}'.format(
                nr_elapsed, nr_obj))

            _, ub_elapsed, soln, _, _ = cvxpy_solve_additive(
                G, lam=0, lp=False, bound_lam=True, acceptablegap=1e-3, callback=False)
            var_cost = np.multiply(G.var, soln).dot(soln)
            lam_high = G.lambar / (2.0 * np.sqrt(var_cost))
            bs_obj, bs_elapsed, _, bs_iter_objs, bs_iter_elapsed, bs_avg_infeas, bs_lams, bs_fs = bs_cvxpy(
                G, lp=False, low=lam_low, high=lam_high, solver_obj=solver_obj)

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

                # solver_iter_objs, solver_iter_elapsed, presolve_time, ordering_time = parse_cplex_log()

                keys = ['network_name', 'solver_status', 'solver_feasible', 'solver_elapsed', 'solver_obj', 'solver_duals', 'solver_primals', 'solver_is_feasible', 'solver_primal_infeas', 'solver_iter_times', 'solver_feasible_dual', 'solver_feasible_primal', 'solver_infeas',
                        'nr_elapsed', 'nr_iter_elapsed', 'nr_iter_objs', 'nr_obj',
                        'nr_xi_times', 'nr_var_times', 'mean', 'var', 'nr_avg_infeas', 'nr_lams', 'nr_fs', 'lambar', 'm', 'n',
                        'elapsed_lower_bound', 'elapsed_upper_bound',
                        'bs_obj', 'bs_iter_elapsed', 'bs_iter_objs', 'bs_elapsed', 'bs_avg_infeas', 'bs_lams', 'bs_fs']

                values = [network[:network.find('.')], status, feasible, solver_elapsed, solver_obj, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal, solver_infeas,
                          nr_elapsed, nr_iter_elapsed, nr_iter_objs, nr_obj, nr_xi_times,
                          nr_var_times, mean, var, nr_avg_infeas, nr_lams, nr_fs, G.lambar, G.m, G.n, lb_elapsed, ub_elapsed,
                          bs_obj, bs_iter_elapsed, bs_iter_objs, bs_elapsed, bs_avg_infeas, bs_lams, bs_fs]



                run_dict = dict(zip(keys, values))

                save_run(run_name, run_dict)
                del run_dict
                del values

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

            lam_dir = str(lambar).replace('.', '_')

            run_name = experiment + '/' + generator + \
                '/' + network[:network.find('.')] + '_' + lam_dir

            solver_weird = False
            status, feasible, reliable_obj, reliable_elapsed, solver_soln, solver_infeas, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal = cvxpy_solve(
                G, cvxpy=cvxpy, fusion=fusion, run_name=run_name)

            solver_soln = solver_soln[:-1]
            rel_mean = G.mu.dot(solver_soln)
            var_cost = np.multiply(G.var, solver_soln).dot(solver_soln)
            rel_std = np.sqrt(var_cost)

            if solver_infeas > 2:
                solver_weird = True

            print('reliable_obj is {}'.format(reliable_obj))
            orig_var = copy.deepcopy(G.var)

            G.var = np.zeros(G.m)

            status, feasible, obj, base_elapsed, base_x, solver_infeas, solver_duals, solver_primals, solver_is_feasible, solver_primal_infeas, solver_iter_times, solver_feasible_dual, solver_feasible_primal = cvxpy_solve(
                G, cvxpy=cvxpy, fusion=fusion, run_name=run_name)
            base_x = base_x[:-1]
            var_cost = np.multiply(orig_var, base_x).dot(base_x)
            base_obj = G.mu.dot(base_x) + G.lambar * (np.sqrt(var_cost))
            base_mean = G.mu.dot(base_x)
            base_std = np.sqrt(var_cost)

            print('base_obj is {}'.format(base_obj))

            t = PrettyTable(['Approach', 'Objective', 'Rel-Gap', 'Mean', 'Std'])
            t.add_row(['Reliable', reliable_obj, '-', rel_mean, rel_std])
            t.add_row(
                ['Base', base_obj, (base_obj - reliable_obj) / reliable_obj, G.mu.dot(base_x), np.sqrt(var_cost)])
            print(t)

            keys = ['network_name', 'G', 'base_elapsed', 'base_obj', 'base_mean', 'base_std', 'rel_mean', 'rel_std',
                    'solver_obj', 'solver_elapsed', 'lambar', 'm', 'n', 'solver_weird']

            values = [network[:network.find('.')], G, base_elapsed, base_obj, base_mean, base_std, rel_mean, rel_std,
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

cvxpy = False
fusion = False

# if cvxpy:
#     import cvxpy as cp
if fusion:
    import mosek.fusion as mf
else:
    import cplex



# print('starting lambar experiments')
# lambars = [0.01, 10, 1000]
# tails = ['f.min', 'g.min', 'h.min', 'j.min', 'k.min']
# exponents = ['12']
# types = ['_sr_', '_lo_8_', '_lo_sr_', '_8_']
# experiment = 'varying_lambar'
# networks = get_networks(experiment, tails, exponents, types, lams=lambars)
# varying_lambda_experiment(networks, lambars)


# print("lambar: ",10)
# print("loose interval")

# lambar = 10
# tails = ['e.min']
# exponents = (np.arange(12, 13)).astype(str)
# types = ['_sr_']
# experiment = 'graph_families'
# networks = get_networks(experiment, tails, exponents, types)
# graph_family_experiment(networks, lambar, cvxpy=cvxpy, fusion=fusion, which='cp')


lambar = 10
tails = ['e.min']
exponents = (np.arange(12, 13)).astype(str)
types = ['_sr_']
experiment = 'graph_families'
networks = get_networks(experiment, tails, exponents, types)
graph_family_experiment(networks, lambar, cvxpy=cvxpy, fusion=fusion, which='cp')


# lambar = 10
# tails = ['e.min', 'd.min', 'c.min', 'b.min', 'a.min']
# exponents = (np.arange(12, 15)).astype(str)
# types = ['_sr_', '_lo_8_', '_lo_sr_', '_8_']
# experiment = 'graph_families'
# networks = get_networks(experiment, tails, exponents, types)
# graph_family_experiment(networks, lambar, cvxpy=cvxpy, fusion=fusion)




# print('starting base vs reliable experiments')
# lambars = np.logspace(-1, 3, 200)
# exponents = ['10']
# types = ['_8_']
# tails = ['e.min']
# experiment = 'base_vs_reliable'
# networks = get_networks(experiment, tails, exponents, types)
# base_vs_reliable_expr(networks, lambars)


