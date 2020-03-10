
import cplex
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

        if self.is_primal_feasible():
            print('feasible')
            self.feasible_primal.append(primal)
            self.feasible_dual.append(dual)
            rel_gap = abs(dual - primal) / primal
            if rel_gap < self.acceptablegap:
                print('dual: ', dual)
                print('primal: ', primal)
                self.aborted = 1
                self.abort()

def cvxpy_solve(G, cvxpy=False, fusion=False, run_name='.', acceptablegap=1e-7):
    start = time.time()

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
    print(sol.status[sol.get_status()])
    obj = prob.solution.get_objective_value()
    soln = np.array(prob.solution.get_values())

    elapsed = time.time() - start

    return obj, elapsed, soln, cplex_cb.duals, cplex_cb.primals, cplex_cb.is_feasible, cplex_cb.primal_infeas, cplex_cb.iter_times, cplex_cb.feasible_dual, cplex_cb.feasible_primal
