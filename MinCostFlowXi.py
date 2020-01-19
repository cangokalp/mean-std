import pyomo
import pyomo.opt
import pyomo.environ as pe
from pyomo.environ import Suffix
from pyomo.core import Var, Constraint

class MinCostFlowXi:
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

        self.node_set = self.node_data.index.unique()
        self.arc_set = self.arc_data.index.unique()
        
        neg_dom_data = self.arc_data[self.arc_data['Capacity']==self.arc_data['Flow']]
#         pos_dom_data = self.arc_data[self.arc_data['Flow']==0]
        eq_zer_data = self.arc_data[self.arc_data['Flow']==0]
    
        if len(eq_zer_data > 0):
            self.eq_zer_set = eq_zer_data.index.unique()
        else:
            self.eq_zer_set = []

        
        if len(neg_dom_data > 0):
            self.neg_dom_set = neg_dom_data.index.unique()
        else:
            self.neg_dom_set = []
            
#         self.pos_dom_set = pos_dom_data.index.unique()
        # self.eq_zer_set = eq_zer_data.index.unique()
#         print(eq_zer_data)
    
        self.lam = lam
        
        self.createModel()
    
    def createModel(self):
        """Create the pyomo model"""
        self.m = pe.ConcreteModel()

        # Create sets
        self.m.neg_dom_set = pe.Set(initialize=self.neg_dom_set, dimen=2)
#         self.m.pos_dom_set = pe.Set(initialize=self.pos_dom_set, dimen=2)
        self.m.eq_zer_set = pe.Set(initialize=self.eq_zer_set, dimen=2)
        self.m.node_set = pe.Set(initialize=self.node_set)
        self.m.arc_set = pe.Set(initialize=self.arc_set , dimen=2)
        
        # Create variables
        self.m.Y = pe.Var(self.m.arc_set, domain=pe.Reals)
        
        def flow_eq_cap(m, n1, n2):
            e = (n1,n2)
#             if self.arc_data.loc[e, 'Flow'] == self.arc_data.loc[e, 'Capacity']:
            return m.Y[e] <= 0
        
#         def flow_eq_zero(m, n1, n2):
#             e = (n1,n2)
# #             if self.arc_data.loc[e, 'Flow'] == 0:
#             return m.Y[e] >= 0

        def flow_noflow(m, n1, n2):
            e = (n1,n2)
            return m.Y[e] == 0
            
        # Lower bounds rule
        def lower_bounds_rule(m, n1, n2):
            e = (n1,n2)
            return m.Y[e] >= self.arc_data.loc[e, 'LowerBound']

        # Upper bounds rule
        def upper_bounds_rule(m, n1, n2):
            e = (n1,n2)
            return m.Y[e] <= self.arc_data.loc[e, 'UpperBound']

        # Flow Balance rule
        def flow_bal_rule(m, n):
            arcs = self.arc_data.reset_index()
            preds = arcs[ arcs.End == n ]['Start']
            succs = arcs[ arcs.Start == n ]['End']
            return sum(m.Y[(p,n)] for p in preds) - sum(m.Y[(n,s)] for s in succs) == self.node_data.loc[n,'Imbalance']

        # Objective rule for the mean variance problem
        def obj_rule(m):
            return sum(self.lam * m.Y[e]**2 * self.arc_data.loc[e, 'Var_Cost'] +
                        2*m.Y[e] * self.arc_data.loc[e, 'Flow']* self.arc_data.loc[e, 'Var_Cost'] for e in self.arc_set)
        
        self.m.OBJ = pe.Objective(rule=obj_rule, sense=pe.minimize)

        self.m.FlowBal = pe.Constraint(self.m.node_set, rule=flow_bal_rule)
        self.m.FlowCap = pe.Constraint(self.m.neg_dom_set, rule=flow_eq_cap)
#         self.m.FlowZer = pe.Constraint(self.m.pos_dom_set, rule=flow_eq_zero)
        self.m.NoFlow = pe.Constraint(self.m.eq_zer_set, rule=flow_noflow)
        self.m.UpperBound = pe.Constraint(self.m.arc_set, rule=upper_bounds_rule)
        self.m.LowerBound = pe.Constraint(self.m.arc_set, rule=lower_bounds_rule)

    def solve(self):
        """Solve the model."""
        solver = pyomo.opt.SolverFactory("cplex")
        self.m.slack = Suffix(direction=Suffix.IMPORT)
        results = solver.solve(self.m, tee=False, keepfiles=False)

        self.m.solutions.store_to(results)

        
        if (results.solver.status != pyomo.opt.SolverStatus.ok):
            logging.warning('Check solver not ok?')
        if (results.solver.termination_condition != pyomo.opt.TerminationCondition.optimal):  
            logging.warning('Check solver optimality?') 
        return results