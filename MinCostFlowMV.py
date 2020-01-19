import pandas
import pyomo
import pyomo.opt
import pyomo.environ as pe
from pyomo.environ import Suffix
from pyomo.core import Var, Constraint

class MinCostFlowMV:
    """This class implements a standard min-cost-flow model.  
    It takes as input two csv files, providing data for the nodes and the arcs of the network.  
    The nodes file should have columns: Node, Imbalance
    that specify the node name and the flow imbalance at the node.  The arcs file should have columns:
    Start, End, Cost, UpperBound, LowerBound
    that specify an arc start node, an arc end node, a cost for the arc, and upper and lower bounds for the flow."""
    
    def __init__(self, node_data, arc_data, lam=1):
        """Read in the csv data."""
        self.node_data = node_data
        self.node_data.set_index(['Node'], inplace=True)
        self.node_data.sort_index(inplace=True)
        
        self.arc_data = arc_data
        self.arc_data.set_index(['Start','End'], inplace=True)
        self.arc_data.sort_index(inplace=True)
        self.arc_set = self.arc_data.index.unique()

        self.node_set = self.node_data.index.unique()
        
        self.lam = lam
        self.createModel()
    
    def createModel(self):
        """Create the pyomo model"""
        self.m = pe.ConcreteModel()

        # Create sets
        self.m.node_set = pe.Set(initialize=self.node_set)
        self.m.arc_set = pe.Set(initialize=self.arc_set , dimen=2)

        # Create variables
        self.m.Y = pe.Var(self.m.arc_set, domain=pe.Reals)
 
        # Lower bounds rule
        def lower_bounds_rule(m, n1, n2):
            e = (n1,n2)
            return m.Y[e] >= self.arc_data.loc[e, 'LowerBound']

        # Upper bounds rule
        def upper_bounds_rule(m, n1, n2):
            e = (n1,n2)
            return m.Y[e] <= self.arc_data.loc[e, 'UpperBound']
        
        def flow_bal_rule(m, n):
            arcs = self.arc_data.reset_index()
            preds = arcs[ arcs.End == n ]['Start']
            succs = arcs[ arcs.Start == n ]['End']
            return sum(m.Y[(p,n)] for p in preds) - sum(m.Y[(n,s)] for s in succs) ==  - self.node_data.loc[n,'Imbalance']

        # Objective rule for the mean variance problem
        def obj_rule(m):
            return sum(m.Y[e]*self.arc_data.loc[e, 'Mean_Cost']  + 
                       self.lam*(m.Y[e]**2) * self.arc_data.loc[e, 'Var_Cost']  for e in self.arc_set)
        
        self.m.OBJ = pe.Objective(rule=obj_rule, sense=pe.minimize)

        self.m.FlowBal = pe.Constraint(self.m.node_set, rule=flow_bal_rule)
        self.m.UpperBound = pe.Constraint(self.m.arc_set, rule=upper_bounds_rule)
        self.m.LowerBound = pe.Constraint(self.m.arc_set, rule=lower_bounds_rule)
        
    def solve(self):
        """Solve the model."""
        solver = pyomo.opt.SolverFactory("cplex")
        self.m.slack = Suffix(direction=Suffix.IMPORT)
        results = solver.solve(self.m, tee=False, keepfiles=False, symbolic_solver_labels=True)

        self.m.solutions.store_to(results)
#         log_infeasible_constraints(self.m)
        
        if (results.solver.status != pyomo.opt.SolverStatus.ok):
            logging.warning('Check solver not ok?')
        if (results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal):  
            logging.warning('Check solver optimality?') 
        return results