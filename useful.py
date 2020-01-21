def create_random_graph_uncap(**kwargs):
    """
    Generates random graph with the given network parameters
    """
    seed = kwargs['seed']
    num_nodes = kwargs['num_nodes']
    mu_top = kwargs['mu_top']
    var_top = kwargs['var_top'] * kwargs['var_scalar']
    d_top = kwargs['d_top']
    arcs = kwargs['num_arcs']
    arc_pnode = int(arcs / num_nodes)

    # fix random seed for reproducibility
    np.random.seed(seed)
    random.seed(seed)
    G = MCF_DiGraph(kwargs['lambar'])

    nodes = np.arange(num_nodes) + 1
    nnodes = len(nodes)

    # that a feasible flow exists
    for node in nodes:
        if node != len(nodes):
            G.nxg.add_edge(node, node + 1, capacity=np.inf,
                           mu=mu_top, var=var_top)
    max_node = nodes[-1]
    # add additional arcs, and assign capacity, mean, variance - by sampling
    # from a U[0, max_param_value]
    while np.array(list(dict(G.nxg.out_degree(nodes)).values())).mean() <= arc_pnode:
        d = round(np.random.uniform(0, d_top), 2)
        mu = round(np.random.uniform(0, mu_top), 4)
        var = round(np.random.uniform(0, var_top), 4)
        src = np.random.randint(1, max_node - 1)
        dest_list = np.arange(max_node) + 1
        dest_list = np.delete(dest_list, src - 1)
        dest = random.choice(dest_list)
        if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)):
            G.nxg.add_edge(src, dest, capacity=np.inf, mu=mu, var=var)

    # set supply demand at the nodes
    G.nxg.node[1]['demand'] = d_top
    G.nxg.node[len(nodes)]['demand'] = -d_top
    for i in range(2, len(nodes)):
        G.nxg.node[i]['demand'] = 0

    G.find_feasible_flow()
    G.set_weight_uc(kwargs['lambar'])
    return G


def read_and_create_graph(**kwargs):
    mu_top = kwargs['mu_top']
    var_top = kwargs['var_top'] * kwargs['var_scalar']
    d_top = kwargs['d_top']
    G = MCF_DiGraph(kwargs['lambar'])
    G.nxg = nx.DiGraph()

    df = pandas.read_csv('austin_net.csv')

    for index, row in df.iterrows():

        init = row['Init_Node']
        term = row['Term_Node']
        e = (init, term)
        cov_coef = np.random.uniform(0.15, 0.3)
        mu = row['Free_Flow_Time'] * 60
        std = mu * cov_coef
        var = std**2
        var = round(var, 4)
        mu = round(mu, 4)
        # cap = np.random.uniform(d_top/2, d_top)
        cap = d_top

        # if (not G.nxg.has_edge(init, term)) and (not G.nxg.has_edge(term, init)) and (init != term):
        # if init ==1 or term==6661:
        # cap = d_top

        G.nxg.add_edge(*e, capacity=cap, mu=mu, var=var)

    nodes = list(G.nxg.nodes())

    nx.set_node_attributes(G.nxg, 0, 'demand')
    # pdb.set_trace()
    G.nxg.node[1]['demand'] = d_top
    G.nxg.node[6661]['demand'] = -d_top
    # hopp = list(nx.all_simple_paths(G.nxg, 1, 6661))
    # pdb.set_trace()

    G.find_feasible_flow(d_top, mu_top, var_top)
    return G


def create_random_graph(**kwargs):

    # changes
    """
    Generates random graph with the given network parameters
    """
    seed = kwargs['seed']
    num_nodes = kwargs['num_nodes']
    mu_top = kwargs['mu_top']
    var_top = kwargs['var_top'] * kwargs['var_scalar']
    d_top = kwargs['d_top']
    arcs = kwargs['num_arcs']
    arc_pnode = int(arcs / num_nodes)

    # fix random seed for reproducibility
    np.random.seed(seed)
    random.seed(seed)
    G = MCF_DiGraph(kwargs['lambar'])

    nodes = np.arange(num_nodes) + 1
    nnodes = len(nodes)

    all_nodes = np.arange(nnodes)[1:][:-1]
    np.random.shuffle(all_nodes)
    spath = np.insert(all_nodes, 0, 0)
    spath = np.insert(spath, nnodes - 1, nnodes - 1)
    G.nxg = nx.DiGraph()
    G.nxg.add_path(spath)

    for u, v, e in G.nxg.edges(data=True):
        e['capacity'] = d_top
        e['mu'] = mu_top
        e['var'] = (mu_top * 0.3)**2
    # pdb.set_trace()
    # for j in range(arcs):
    extra = arcs - len(G.nxg.edges())

    for j in range(extra):
        # for src in G.nxg.nodes():
        #     for j in range(arc_pnode - 1):
        succ = False
        mu = np.random.uniform(mu_top)
        cov_coef = np.random.uniform(0.15, 0.3)
        std = mu * cov_coef
        var = std**2

        src = np.random.randint(0, nnodes)

        dest_list = np.arange(nnodes)
        dest_list = np.delete(dest_list, src)

        d = np.random.uniform(d_top * 0.4, d_top)
        # d=d_top

        if src == 0:
            dest_list = np.delete(dest_list, nnodes - 2)
        if src == nnodes - 1:
            dest_list = np.delete(dest_list, 0)

        while not succ:
            dest = random.choice(dest_list)
            if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)) and (src != dest):
                G.nxg.add_edge(src, dest, capacity=d, mu=mu, var=var)
                succ = True

    # for src in G.nxg.nodes():
    #     while G.nxg.degree(src) -2 < arc_pnode:
    #     # for j in range(arc_pnode):
    #         succ = False
    #         mu = np.random.uniform(mu_top)
    #         cov_coef = np.random.uniform(0.15, 0.3)
    #         std = mu*cov_coef
    #         var = std**2
    #         var = round(var, 4)

    #         dest_list = np.arange(nnodes)
    #         dest_list = np.delete(dest_list, src)

    #         if src == 0:
    #             dest_list = np.delete(dest_list, nnodes-2)
    #         if src == nnodes-1:
    #             dest_list = np.delete(dest_list, 0)

    #         while not succ:
    #             dest = random.choice(dest_list)
    #             if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)) and (src != dest) and (G.nxg.degree(dest)-2 < arc_pnode):
    #                 G.nxg.add_edge(src, dest, capacity=d_top, mu=mu, var=var)
    #                 succ = True
    print('num edges: ', len(G.nxg.edges()))

    # while np.array(list(dict(G.nxg.degree(nodes)).values())).mean() -1 <=
    # arc_pnode :

    #     d = d_top
    #     mu = np.random.uniform(mu_top)
    #     cov_coef = np.random.uniform(0.15, 0.3)
    #     std = mu*cov_coef
    #     var = std**2
    #     var = round(var, 4)

    #     src = np.random.randint(0, nnodes)
    #     dest = np.random.randint(0, nnodes)

    #     dest_list = np.arange(nnodes)
    #     dest_list = np.delete(dest_list, src)

    #     if src == 0:
    #         dest_list = np.delete(dest_list, nnodes-2)
    #     if src == nnodes-1:
    #         dest_list = np.delete(dest_list, 0)

    #     dest = random.choice(dest_list)

    #     if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)) and (src != dest):
    #         G.nxg.add_edge(src, dest, capacity=d, mu=mu, var=var)

    # for src in G.nxg.nodes():
    #     for j in range(arc_pnode):
    #         succ = False
    #         d = np.random.uniform(d_top)
    #         mu = np.random.uniform(mu_top)
    #         cov_coef = np.random.uniform(0.15, 0.3)
    #         std = mu*cov_coef
    #         var = std**2
    #         var = round(var, 4)
    #         # var = np.random.uniform(var_top)

    #         dest_list = np.arange(nnodes)
    #         dest_list = np.delete(dest_list, src)

    #         if src == 0:
    #             dest_list = np.delete(dest_list, nnodes-2)
    #         if src == nnodes-1:
    #             dest_list = np.delete(dest_list, 0)

    #         while not succ:
    #             dest = random.choice(dest_list)
    #             if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)) and (src != dest):
    #                 G.nxg.add_edge(src, dest, capacity=d_top, mu=mu, var=var)
    #                 succ = True

    G.nxg.node[0]['demand'] = d_top
    G.nxg.node[nnodes - 1]['demand'] = -d_top
    for i in range(1, nnodes - 1):
        G.nxg.node[i]['demand'] = 0
    G.init_node = 0
    G.end_node = nnodes - 1
    G.find_feasible_flow(d_top, mu_top, var_top)

    # gg = nx.connected_watts_strogatz_graph(100, 4*2, 1, tries=100, seed=None)
    # can = list(nx.all_simple_paths(gg, 0, 99))
    # pdb.set_trace()

    # gg = nx.connected_watts_strogatz_graph(num_nodes, arc_pnode*2, 0.2, tries=100, seed=None)
    # print(list(gg.selfloop_edges()))
    # print('edge number: ', len(gg.edges()))
    # G.nxg = gg

    # for u,v,e in G.nxg.edges(data=True):
    #     # if v-u != 1:
    #     # for edge in G.nxg.edges():
    #         # G.nxg.add_edge(src, dest, capacity=d, mu=mu, var=var)
    #     # can = np.random.randint(1, 3)
    #     # if can == 2 or can == 3:
    #     #     # var = round(np.random.uniform(4,var_top), 4)
    #     #     var = min(max(np.random.normal(var_top*0.9, var_top*0.2, 1)[0],2), var_top)
    #     #     mu = min(max(np.random.normal(mu_top*0.9, mu_top*0.2, 1)[0],2), mu_top)
    #     # else:
    #     #     # var = round(np.random.uniform(var_top), 4)
    #     #     mu = min(max(np.random.normal(mu_top*0.9, mu_top*0.2, 1)[0],2), mu_top)
    #     #     var = min(max(np.random.normal(var_top*0.2, var_top*0.1, 1)[0],2), var_top)

    #     # d = min(max(np.random.normal(d_top, d_top*0.5,1)[0],0), d_top)
    #     # d = d_top
    #     d = np.random.uniform(d_top)
    #     mu = np.random.uniform(mu_top)
    #     var = np.random.uniform(var_top)

    #     try:

    #         if e['capacity'] == d_top:
    #             print('hiii')
    #             continue
    #         else:
    #             e['capacity']= d
    #             e['mu'] =  mu
    #             e['var'] = var
    #     except:
    #         e['capacity']= d
    #         e['mu'] =  mu
    #         e['var'] = var

    # pdb.set_trace()

    # for node in nodes:
    #     if node != len(nodes):
    #         mu = min(max(np.random.normal(mu_top/2, mu_top*0.33, 1)[0],2), mu_top)
    #         var = min(max(np.random.normal(var_top/7, var_top*0.1, 1)[0],2), var_top)
    #         d = min(max(np.random.normal(d_top/1.2, d_top*0.33,1)[0],0), d_top)
    #         G.nxg.add_edge(node, node + 1, capacity=d,
    #                        mu=mu, var=var)
    # # edgeno = len(G.nxg.edges())
    # max_node = nodes[-1]
    # # print(edgeno)
    # # add additional arcs, and assign capacity, mean, variance - by sampling
    # # from a U[0, max_param_value]
    # while np.array(list(dict(G.nxg.degree(nodes)).values())).mean() <= arc_pnode :
    #     # d = round(np.random.uniform(d_top),4)
    #     # d = round(np.random.uniform(d_top), 2)
    #     # mu = round(np.random.uniform(4,mu_top), 4)
    #     # mu = round(np.random.uniform(mu_top), 4)
    #     # var = round(np.random.uniform(var_top), 4)

    #     src = np.random.randint(1, max_node+1)
    #     dest = np.random.randint(1, max_node+1)

    #     diff = dest - src

    #     can = np.random.randint(1, 3)
    #     if can == 1:
    #         # var = round(np.random.uniform(4,var_top), 4)
    #         var = min(max(np.random.normal(var_top/2, var_top*0.33, 1)[0],2), var_top)
    #         mu = min(max(np.random.normal(mu_top/7, mu_top*0.1, 1)[0],2), mu_top)
    #     else:
    #         # var = round(np.random.uniform(var_top), 4)
    #         mu = min(max(np.random.normal(mu_top/2, mu_top*0.33, 1)[0],2), mu_top)
    #         var = min(max(np.random.normal(var_top/7, var_top*0.1, 1)[0],2), var_top)

    #     d = min(max(np.random.normal(d_top/1.2, d_top*0.2,1)[0],0), d_top)

    #     # print(d, mu, var)

    #     # dest_list = np.arange(max_node) + 1
    #     # # dest_list = np.delete(dest_list, src - 1)
    #     # dest = random.choice(dest_list)
    #     # if dest-src == 1:
    #     #     G.nxg.add_edge(src, dest, capacity=d, mu=mu, var=var)
    #     # else:
    #     if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)) and (src != dest):
    #         G.nxg.add_edge(src, dest, capacity=d, mu=mu, var=var)

    # # pdb.set_trace()
    # # while nx.is_strongly_connected(G.nxg):
    # #     arc_pnode +=1
    # #     print(arc_pnode)
    # #     while np.array(list(dict(G.nxg.out_degree(nodes)).values())).mean() <= arc_pnode :
    # #         d = round(np.random.uniform(0, d_top), 2)
    # #         mu = round(np.random.uniform(0, mu_top), 4)
    # #         var = round(np.random.uniform(0, var_top), 4)
    # #         src = np.random.randint(1, max_node - 1)
    # #         dest_list = np.arange(max_node) + 1
    # #         dest_list = np.delete(dest_list, src - 1)
    # #         dest = random.choice(dest_list)
    # #         if (not G.nxg.has_edge(src, dest)) and (not G.nxg.has_edge(dest, src)):
    # #             G.nxg.add_edge(src, dest, capacity=d, mu=mu, var=var)
    # pdb.set_trace()

    # # set supply demand at the nodes
    # # G.nxg.node[1]['demand'] = d_top
    # G.nxg.node[0]['demand'] = d_top
    # G.nxg.node[len(nodes)-1]['demand'] = -d_top
    # # for i in range(2, len(nodes)):
    # for i in range(1, len(nodes)-1):
    # for i in range(1, len(nodes)-1):

    #     # try:
    #     G.nxg.node[i]['demand'] = 0
    #     # except:
    #     #     if (not G.nxg.has_edge(i-1, i)) and (not G.nxg.has_edge(i, i-1)):
    #     #         G.nxg.add_edge(i-1, i, capacity=d_top,
    #     #                        mu=mu_top*10000, var=var_top*10000)
    #     #     if (not G.nxg.has_edge(i, i+1)) and (not G.nxg.has_edge(i+1, i)):
    #     #         G.nxg.add_edge(i, i+1, capacity=d_top,
    #     #                        mu=mu_top*10000, var=var_top*10000)
    #     #     print('added decoy link')
    #     #     G.nxg.node[i]['demand'] = 0
    # G.init_node = 0
    # G.end_node = nnodes-1

    # G.find_feasible_flow()
    # G.set_weight_uc(kwargs['lambar'])

    # m = len(G.nxg.edges())
    # n = len(G.nxg.nodes())
    # f = open('10k_20k.tntp', 'w+')
    # f.write('<NUMBER OF ZONES>' + str(n) + '\n')
    # f.write('<NUMBER OF NODES>' + str(n) + '\n')
    # f.write('<FIRST THRU NODE> 1 \n')
    # f.write('<NUMBER OF LINKS>' + str(m) + '\n')

    # f.write('<END OF METADATA>\n')
    # f.write('~   Init node   Term node   Capacity    Length  Free Flow Time
    # Alpha   Beta    Speed limit     Toll    Type    ;\n')

    # nodelist = [n for n in G.nxg.nodes()]

    # for u, v, e in G.nxg.edges(data=True):
    #     try:
    #         curvar = G.nxg[u][v]['var']
    #         w = 2 * curvar * G.nxg[u][v]['flow'] + \
    #             2 * G.lam * G.nxg[u][v]['xi'] * curvar
    #         f.write(str(u) + '       ' + str(v) + '       ' +
    #                 str(w) + '       ' + str(u) + '\n')
    #     except:
    #         curvar = G.nxg[v][u]['var']
    #         w = 2 * curvar * G.nxg[v][u]['flow'] + \
    #             2 * G.lam * G.nxg[v][u]['xi'] * curvar
    #         w = w / 100.0
    #         f.write(str(u) + '       ' + str(v) + '       ' + str((-1.0)
    #                                                               * w) + '       ' + str(u) + '\n')

    # f.write('@attributes\n')
    # f.write('source ' + str(nodelist[0]) + '\n')
    # f.write('target ' + str(nodelist[-1]) + '\n')
    # f.close()

    return G
def cplex_solve(**kwargs):

    c = cplex.Cplex()

    if kwargs['G'] != None:
        G = copy.deepcopy(kwargs['G'])
    else:
        G = create_random_graph(**kwargs)

    lam_bar = kwargs['lam_bar']

    G.set_lambda(lam=lam_bar)
    cvx_time_st = time.time()

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
    all_uv = []
    pdb.set_trace()
    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu', 0)
        sigma[i] = np.sqrt(e.get('var', 0))
        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)

        # rows.append(u - 1)
        rows.append(u)
        cols.append(i)
        values.append(1)
        # rows.append(v - 1)
        rows.append(v)
        cols.append(i)
        values.append(-1)
        A[u - 1, i]=1
        A[v - 1, i]=-1
        all_uv.append(u)
        all_uv.append(v)

        i += 1
    tot = i
    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

    cone = {}

    c.variables.add(names = [ "x1", "x2", "x3", "x4", "x5", "x6" ])
    c.variables.set_lower_bounds([("x2", -cplex.infinity),
                                  ("x3", -cplex.infinity),
                                  ("x5", -cplex.infinity)])

    # Create objective function.
    c.objective.set_linear([("x1", 1), ("x2", 1), ("x3", 1),
                            ("x4", 1), ("x5", 1), ("x6", 1)])

    # Create constraints.
    # c.linear_constraints.add(lin_expr = [[["x1", "x2", "x5"], [1, 1, 1]],
    #                                      [["x3", "x5", "x6"], [1, 1, 1]]],
    #                          senses = ['E', 'E'],
    #                          rhs = [8, 10],
    #                          names = ["c1","c2"])
    c.quadratic_constraints.add(quad_expr = [["x1", "x2", "x3"],
                                             ["x1", "x2", "x3"],
                                             [-1,   1,    1   ]],
                                sense = 'L',
                                rhs = 0,
                                name = "q1")
    c.quadratic_constraints.add(quad_expr = [["x4", "x5"],
                                             ["x4", "x5"],
                                             [-1,   1   ]],
                                sense = 'L',
                                rhs = 0,
                                name = "q2")

    # Set up the cone map.
    # cone["x1"] = "q1"
    # cone["x2"] = NOT_CONE_HEAD
    # cone["x3"] = NOT_CONE_HEAD
    # cone["x4"] = "q2"
    # cone["x5"] = NOT_CONE_HEAD
    # cone["x6"] = NOT_IN_CONE



    c.write("model_test.lp")
    pdb.set_trace()
    c.minimize()

    # c.objective.set_sense(c.objective.sense.minimize)

    pdb.set_trace()

    my_obj = []

    c.variables.add(names=["x{0}".format(i) for i in range(m)], lb=np.zeros(m), ub=cap)

    name2idx = { n : j for j, n in enumerate(c.variables.get_names())}

    c.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [cols[j] for j in [index for index,value in enumerate(rows) if value == i]], val = [A[i,j] for j in [cols[index] for index,value in enumerate(rows) if value == i]]) for i in list(set(rows))],
                                                          senses = ["E"] * n, rhs = d)
    
    c.write("model.lp")
    #check if right
    # c.linear_constraints.add(senses=["E"]*n, rhs=d)
    # c.linear_constraints.set_coefficients(zip(rows, cols, values))



    # const_time = time.time()

    # x = model.variable("x", m, Domain.greaterThan(0.0))
    # x_n = model.variable("x_n", m)
    # decoy = model.variable("decoy", 1)
    # decoy_2 = Var.vstack(decoy, x_n)

    # model.objective("myobj",
    #                 ObjectiveSense.Minimize,
    #                 Expr.mul(kwargs['obj_scalar'],
    #                          Expr.add(Expr.mul(kwargs['mu_scalar'], Expr.dot(mu, x)), Expr.mul(lam_bar, decoy))))
    # model.constraint(Expr.sub(x_n, Expr.mulElm(sigma, x)),
    #                  Domain.equalsTo(0.0))
    # model.constraint(decoy_2, Domain.inQCone())

    # # if kwargs['max_iters'] == None:
    # model.setSolverParam("intpntCoTolRelGap", 1.0e-9)
    # model.setSolverParam("intpntCoTolPfeas", 1.0e-9)
    # model.setSolverParam("intpntCoTolDfeas", 1.0e-9)
    # model.setSolverParam("intpntCoTolMuRed", 1.0e-9)
    # model.setSolverParam("intpntMaxIterations", 100000)
    # # else:
    # #     model.setSolverParam("intpntMaxIterations", kwargs['max_iters'])

    # import logging

    # logging = Logger(**kwargs)
    # model.setLogHandler(logging)

    # model.setSolverParam("logIntpnt", 1000)
    # model.setSolverParam("log", 1000)
    # model.setSolverParam("logFile", 1000)

    # # tm = model.getSolverDoubleInfo("optimizerTime")
    # # om = model.getSolverDoubleInfo("intpntPrimalObj")
    # # it = model.getSolverIntInfo("intpntIter")
    # # print('Time: {0}\nIterations: {1}\nObjective {2}'.format(tm, it, om))
    # const_time = time.time() - const_time
    # solve_time = time.time()

    # model.solve()

    # solve_time = time.time() - solve_time
    # mosek_overhead = const_time
    # # tm = model.getSolverDoubleInfo("optimizerTime")
    # # om = model.getSolverDoubleInfo("intpntPrimalObj")
    # # it = model.getSolverIntInfo("intpntIter")
    # # print('Time: {0}\nIterations: {1}\nObjective {2}'.format(tm, it, om))
    # # pdb.set_trace()
    # # print(model.primalObjValue())
    # # print('const_time: ', const_time)
    # # print('total_time_elapsed: ', solve_time + const_time)
    # # print(model.getProblemStatus())

    # cvx_obj = model.primalObjValue()
    # x = x.level()
    # # pdb.set_trace()
    # # try:
    # #     i = 0
    # #     for u, v, e in G.nxg.edges(data=True):
    # #         cur_x = x[i]
    # #         if abs(round(cur_x, 4) - cur_x) <= 1e-8:
    # #             cur_x = round(cur_x, 4)
    # #         G.nxg[u][v]['flow'] = cur_x
    # #         x[i] = cur_x
    # #         i += 1

    # #     cur_r_dict = nx.to_dict_of_dicts(G.nxg)
    # #     cc = [v for v in cur_r_dict.keys()]
    # #     list_ofed = []
    # #     endnode = len(G.nxg.nodes())
    # #     for c in cc:
    # #         if endnode in cur_r_dict[c].keys():
    # #             list_ofed.append(c)

    # #     flows = []
    # #     for e in list_ofed:
    # #         flows.append(G.nxg[e][endnode]['flow'])
    # #     flows = sum(np.array(flows))
    # # except:
    # #     pass

    # # x = nx.get_edge_attributes(G.nxg, 'flow')
    # # cvx_obj = G.tot_cost()
    # cvx_elapsed = solve_time + const_time

    # return cvx_obj, cvx_elapsed, mosek_overhead, x, mu.dot(x), decoy.level()[0]
def cvxpy_solve_xi(G,lam):

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
    F = np.zeros((m,n))
    xi = cp.Variable(m)

    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        x_[i] = e.get('flow',0)
        mu[i] = e.get('mu',0)
        var[i] = e.get('var',0)
        sigma[i] = np.sqrt(e.get('var',0))
        cap[i] = e.get('capacity',0)
        arc_dict[i] = (u, v)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        d[i] = 0
        i += 1

        arc_data = {'Start':[], 'End':[]}

    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)

    constraints = []
    for i in range(m):
        if abs(x_[i]) < 1e-9:
            constraints.append(xi[i]==0)
            # constraints.append(xi[i]<=cap[i])

        elif x_[i]==cap[i]:
            constraints.append(-x_[i]<=xi[i])
            constraints.append(xi[i]<=0)
        else:
            constraints.append(-x_[i]<=xi[i])
            constraints.append(xi[i]<=cap[i] - x_[i])

    for i in range(n):
        preds = []
        succs = []
        for key, vals in arc_dict.items():
            if vals[0] == i:
                preds.append(key)
            if vals[1] == i:
                succs.append(key)

        constraint = sum(xi[p] for p in preds) - sum(xi[s] for s in succs) ==  0
        constraints.append(constraint)


    # Construct the problem.
    objective = cp.Minimize(G.lam*var.T*xi**2 + 2*np.multiply(var,x_).T*xi)
    # print(sum(np.multiply(var,x_)))
    # print(sum(x_))
    prob = cp.Problem(objective, constraints)
    cvx_time_st = time.time()

    result = prob.solve(solver='MOSEK',verbose=False)
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
def cvxpy_solve_xi_new(G,lam):

    cvx_time_st = time.time()
    G.set_lambda(lam=lam)

    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()
    # print(n,m)

    x_ = np.zeros(m)
    mu = np.zeros(m)
    var = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m,n))
    xi = cp.Variable(m)

    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        x_[i] = e.get('flow',0)
        mu[i] = e.get('mu',0)
        var[i] = e.get('var',0)
        sigma[i] = np.sqrt(e.get('var',0))
        cap[i] = e.get('capacity',0)
        arc_dict[i] = (u, v)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        d[i] = 0
        i += 1

        arc_data = {'Start':[], 'End':[]}

    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)

    constraints = []
    for i in range(m):
        if abs(x_[i]) < 1e-6:
            constraints.append(xi[i]==0)
        elif x_[i]==cap[i]:
            constraints.append(xi[i]<=0)

    for i in range(n):
        preds = []
        succs = []
        for key, vals in arc_dict.items():
            if vals[0] == i:
                preds.append(key)
            if vals[1] == i:
                succs.append(key)

        constraint = sum(xi[p] for p in preds) - sum(xi[s] for s in succs) ==  0

        constraints.append(constraint)


    # Construct the problem.
    objective = cp.Minimize(G.lam*var.T*xi**2 + 2*np.multiply(var,x_).T*xi)
    # print(sum(np.multiply(var,x_)))
    # print(sum(x_))
    prob = cp.Problem(objective, constraints)

    result = prob.solve(solver='MOSEK',verbose=False)
    print("status:", prob.status)

    cvx_soln = xi.value
    cvx_obj = objective.value
    cvx_elapsed = time.time() - cvx_time_st

    # i = 0
    # arc_dict = {}
    # summ = 0
    # for u, v, e in G.nxg.edges(data=True):
    #     G.nxg[u][v]['xi'] = cvx_soln[i]
    #     summ += cvx_soln[i]

    #     i += 1
    return cvx_soln, np.multiply(var, x_).dot(cvx_soln)

    # for k in range(m):
    #     constraints.append(x_n[k] == sigma[k] * x[k])

    # for i in range(n):
    #     preds=[]
    #     succs=[]
    #     for key, vals in arc_dict.items():
    #         if vals[0] == i + 1:
    #             preds.append(key)
    #         if vals[1] == i + 1:
    #             succs.append(key)
    #     constraint=sum(x[p] for p in preds) - sum(x[s]
    #                                                 for s in succs) == d[i]
    #     constraints.append(constraint)

# def solve_xi(G,lam):

#     cvx_time_st = time.time()
#     G.set_lambda(lam=lam)

#     m = G.nxg.number_of_edges()
#     n = G.nxg.number_of_nodes()
#     # print(n,m)

#     x_ = np.zeros(m)
#     mu = np.zeros(m)
#     var = np.zeros(m)
#     sigma = np.zeros(m)
#     cap = np.zeros(m)
#     d = np.zeros(n)
#     F = np.zeros((m,n))
#     xi = cp.Variable(m)

#     i = 0
#     arc_dict = {}
#     for u, v, e in G.nxg.edges(data=True):
#         x_[i] = round(e.get('flow',0), 5)
#         mu[i] = e.get('mu',0)
#         var[i] = e.get('var',0)
#         sigma[i] = np.sqrt(e.get('var',0))
#         cap[i] = e.get('capacity',0)
#         arc_dict[i] = (u, v)
#         i += 1

#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         d[i] = dat['demand']
#         i += 1

#         arc_data = {'Start':[], 'End':[]}
#     for u, v, e in G.nxg.edges(data=True):
#         arc_data['Start'].append(u)
#         arc_data['End'].append(v)

#     constraints = []
#     for i in range(m):
#         if abs(x_[i]) < 1e-6:
#             constraints.append(xi[i]==0)
#         elif x_[i]==cap[i]:
#             constraints.append(-x_[i]<=xi[i])
#             constraints.append(-xi[i]<=0)
#         else:
#             constraints.append(-x_[i]<=xi[i])
#             constraints.append(xi[i]<=cap[i] - x_[i])

#     for i in range(n):
#         preds = []
#         succs = []
#         for key, vals in arc_dict.items():
#             if vals[0] == i+1:
#                 preds.append(key)
#             if vals[1] == i+1:
#                 succs.append(key)

#         constraint = sum(xi[p] for p in preds) - sum(xi[s] for s in succs) ==  0
#         constraints.append(constraint)


#     # Construct the problem.
#     objective = cp.Minimize(G.lam*var.T*xi**2 + 2*np.multiply(var,x_).T*xi)
#     # print(sum(np.multiply(var,x_)))
#     # print(sum(x_))
#     prob = cp.Problem(objective, constraints)

#     result = prob.solve(solver='MOSEK',verbose=False)
#     # print("status:", prob.status)

#     cvx_soln = xi.value
#     cvx_obj = objective.value
#     cvx_elapsed = time.time() - cvx_time_st

#     i = 0
#     arc_dict = {}
#     summ = 0
#     for u, v, e in G.nxg.edges(data=True):
#         G.nxg[u][v]['xi'] = cvx_soln[i]
#         summ += cvx_soln[i]

#         i += 1
    # print(summ)



def solve_var(G,lam):

    cvx_time_st = time.time()
    G.set_lambda(lam=lam)

    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    mu = np.zeros(m)
    var = np.zeros(m)
    sigma = np.zeros(m)
    cap = np.zeros(m)
    d = np.zeros(n)
    F = np.zeros((m,n))
    x = cp.Variable(m)
    x_n = cp.Variable(m)

    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        mu[i] = e.get('mu',0)
        var[i] = e.get('var',0)
        sigma[i] = np.sqrt(e.get('var',0))
        cap[i] = e.get('capacity',0)
        arc_dict[i] = (u, v)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

        arc_data = {'Start':[], 'End':[]}
    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)

    constraints = [0 <= x, x <= cap]

    for i in range(n):
        preds = []
        succs = []
        for key, vals in arc_dict.items():
            if vals[0] == i+1:
                preds.append(key)
            if vals[1] == i+1:
                succs.append(key)

        constraint = sum(x[p] for p in preds) - sum(x[s] for s in succs) ==  d[i]
        constraints.append(constraint)

    # Construct the problem.
    objective = cp.Minimize(mu.T*x + lam*var.T*x**2)

    prob = cp.Problem(objective, constraints)

    result = prob.solve(solver='MOSEK',verbose=True)
    print("status:", prob.status)

    cvx_soln = x.value
    cvx_obj = objective.value
    cvx_elapsed = time.time() - cvx_time_st
    print('elapsed time by cvx solve: ', cvx_elapsed)
    # print(cvx_soln)

    i = 0
    arc_dict = {}
    for u, v, e in G.nxg.edges(data=True):
        G.nxg[u][v]['flow'] = cvx_soln[i]
        i += 1


## toy problems
def create_toy_problem(lam):

    # Toy problem generation (With Stochastic Arc Costs) \
    # This creates the graph on pg 318 (Network Flows by Ahuja)

    GAstoch = MCF_DiGraph()

    GAstoch.nxg.add_edge(1,2, weight=1,mu=2,var=3,flow=3,capacity=10)
    GAstoch.nxg.add_edge(1,3, weight=10,mu=2,var=3,flow=3,capacity=10)
    GAstoch.nxg.add_edge(4,2, weight=0,mu=2,var=3,flow=3,capacity=10)
    GAstoch.nxg.add_edge(2,3, weight=3,mu=2,var=3,flow=3,capacity=10)
    GAstoch.nxg.add_edge(3,4, weight=2,mu=2,var=3,flow=3,capacity=10)
    GAstoch.nxg.add_edge(4,1, weight=8,mu=2,var=3,flow=3,capacity=10)

    GAstoch.nxg.node[1]['demand'] = 1
    GAstoch.nxg.node[2]['demand'] = 0
    GAstoch.nxg.node[3]['demand'] = 0
    GAstoch.nxg.node[4]['demand'] = -1

    GAstoch.nxg.node[1]['pos'] = (0,0)
    GAstoch.nxg.node[2]['pos'] = (3,3)
    GAstoch.nxg.node[3]['pos'] = (3,-3)
    GAstoch.nxg.node[4]['pos'] = (6,0)

    # Lets assign it to another to keep the original
    G = GAstoch

    # Set the lambda parameter for mean-variance tradeoff

    return G

def create_nfb_instance():
    # Toy problem generation (With Stochastic Arc Costs) \
    # This creates the graph on pg 318 (Network Flows by Ahuja)

    GAstoch = MCF_DiGraph()
    GAstoch.nxg.add_edge(1,2, weight=2, capacity=4, flow=3)
    GAstoch.nxg.add_edge(1,3, weight=2, capacity=2, flow=1)
    GAstoch.nxg.add_edge(2,3, weight=1, capacity=2, flow=0)
    GAstoch.nxg.add_edge(2,4, weight=3, capacity=3, flow=3)
    GAstoch.nxg.add_edge(3,4, weight=1, capacity=5, flow=1)

    GAstoch.nxg.node[1]['demand'] = 4
    GAstoch.nxg.node[2]['demand'] = 0
    GAstoch.nxg.node[3]['demand'] = 0
    GAstoch.nxg.node[4]['demand'] = -4

    GAstoch.nxg.node[1]['pos'] = (0,0)
    GAstoch.nxg.node[2]['pos'] = (-3,3)
    GAstoch.nxg.node[3]['pos'] = (3,-3)
    GAstoch.nxg.node[4]['pos'] = (3,3)

    G = GAstoch

    return G


 def check_nmcc(self, G):
        """Temporary solution
        """
        print('hey')
        pdb.set_trace()
        tol = 1e-10
        valid = []
        for c in nx.simple_cycles(self.nxg):
            if len(c) == 2:
                pass
            else:
                c.append(c[0])
                valid.append(c)

        if len(valid) == 0:
            return False
        else:

            path_cost = []
            for c in valid:
                cost = 0
                for i in range(len(c)-1):
                    u = c[i]
                    v = c[i+1]
                    try:
                        cost += self.nxg[u][v]['mu'] + 2*G.lam*G[u][v]['flow']*self.nxg[u][v]['var']
                    except:
                        cost += self.nxg[u][v]['mu'] + 2*G.lam*G[v][u]['flow']*self.nxg[u][v]['var']
                path_cost.append(cost)

            if  np.min(path_cost) < 0:
                W = valid[np.argmin(path_cost)]
                self.nmcc_cost = np.min(path_cost)
                self.nmcc = W

                pdb.set_trace()
                print(W)

                delta = self.find_delta()
                if delta <= 0:
                    return False

                if abs(self.nmcc_cost) < tol:
                    self.nmcc_cost = 0.0
                    return False
                else:
                    return True
            else:
                return False

    def get_nmcc(self):

        return self.nmcc


    def nmcc_exists_test(self, G):
            ## Head
            g = self.nxg

            #test#
            if not nx.is_strongly_connected(g):
                print('not strongly connected')
                scc = [c for c in sorted(nx.strongly_connected_components(g),key=len, reverse=True)]
                scc_count = len(scc)
                min_lam = []
                M_l = []
                K_l = []
                v_l =[]
                pi_l = []
                glen_l =[]
                curg_l = []
                dkv_l = []
                for i in range(scc_count):
                    nodelist = scc[i]
                    if len(nodelist) <= 1:
                        continue

                    cur_g = self.nxg.copy()
                    for node, data in self.nxg.edges():
                        if node not in nodelist:
                            cur_g.remove_node(node)

                    n = len(cur_g.nodes())
                    m = len(cur_g.edges())
                    d_kv = np.ones((n+1, n))
                    d_kv = np.inf * d_kv
                    pi = -1*np.ones((n+1, n))
                    Visit = np.ones((n, n))
                    M = np.zeros(n)
                    K = np.zeros(n)
                    v_star = None

                    for v in cur_g.nodes():
                        Visit[v-1, 0] = False
                        Visit[v-1, 1] = False

                    s = self.FindSource()
                    d_kv[0, s-1] = 0

                    pi[0, s-1] = None
                    turn = 0
                    Visit[s-1, turn] = True

                    ## Body
                    for k in range(n):
                        for v in cur_g.nodes():
                            if Visit[v-1, turn] == True:
                                Visit[v-1, turn] = False
                                for u in cur_g.neighbors(v):

                                    w_vu= self.nxg[v][u]['weight']

                                    if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                        d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                        pi[k+1, u-1] = v
                                        Visit[u-1, 1-turn] = True
                        turn = 1 - turn

                    ## Tail
                    lam = np.inf
                    for v in cur_g.nodes():
                        if Visit[v-1, turn] == True:
                            M[v-1] = -np.inf
                            for k in range(n):
                                if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k):
                                    M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)
                                    K[v-1] = k

                            if lam > M[v-1]:
                                lam = M[v-1]
                                v_star = v

                    min_lam.append(lam)
                    M_l.append(M)
                    K_l.append(K)
                    v_l.append(v_star)
                    pi_l.append(pi)
                    glen_l.append(len(cur_g.nodes()))
                    curg_l.append(cur_g)
                    dkv_l.append(d_kv)

                min_ind = np.argmin(min_lam)
                lam = min_lam[min_ind]
                M = M_l[min_ind]
                K = K_l[min_ind]
                v_star = v_l[min_ind]
                pi = pi_l[min_ind]
                glen = glen_l[min_ind]
                g = curg_l[min_ind]
                self.d_kv = dkv_l[min_ind]
                self.M = M
                self.K = K
                self.v_star = v_star
                self.pi = pi
                self.lam_star = lam
                self.g = g
                self.set_nmcc_test(glen, G)
                return lam
            else:

                n = len(g.nodes())
                m = len(g.edges())
                d_kv = np.ones((n+1, n))
                d_kv = np.inf * d_kv
                pi = -1*np.ones((n+1, n))
                Visit = np.ones((n, n))
                M = np.zeros(n)
                K = np.zeros(n)
                v_star = None

                for v in g.nodes():
                    Visit[v-1, 0] = False
                    Visit[v-1, 1] = False

                s = self.FindSource()
                d_kv[0, s-1] = 0

                pi[0, s-1] = None
                turn = 0
                Visit[s-1, turn] = True

                ## Body
                # pdb.set_trace()
                for k in range(n):
                    for v in g.nodes():
                        if Visit[v-1, turn] == True:
                            Visit[v-1, turn] = False
                            for u in g.neighbors(v):
                                w_vu = self.nxg[v][u]['weight']
                                if d_kv[k+1, u-1] > d_kv[k, v-1] + w_vu:
                                    d_kv[k+1, u-1] = d_kv[k, v-1] + w_vu
                                    pi[k+1, u-1] = v
                                    Visit[u-1, 1-turn] = True
                    turn = 1 - turn

                ## Tail
                lam = np.inf
                for v in g.nodes():
                    if Visit[v-1, turn] == True:
                        M[v-1] = -np.inf
                        for k in range(n):
                            if M[v-1] < (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k):
                                M[v-1] = (d_kv[n, v-1] - d_kv[k, v-1])/float(n-k)
                                K[v-1] = k

                        if lam > M[v-1]:
                            lam = M[v-1]
                            v_star = v

                self.M = M
                self.K = K
                self.v_star = v_star
                self.pi = pi
                self.lam_star = lam
                self.g = g
                self.d_kv = d_kv
                self.set_nmcc_test(len(self.nxg.nodes()),G)
                return lam

    def set_nmcc_test(self, glen, G):
        d_kv = self.d_kv
        cycle_len = int(glen - self.K[self.v_star - 1])
        path_to_v = []
        next_v = self.v_star
        print('cyle length: ', cycle_len)
        print('v_star: ', self.v_star)

        # pdb.set_trace()
        if self.lam_star < 0:
            j = 0
            for i in range(glen, 0, -1):
                next_v = int(self.pi[glen-j,next_v-1])
                if i == glen:
                    loopv = next_v
                path_to_v.append(next_v)
                if next_v == loopv:
                    break
                # if len(path_to_v) != len(set(path_to_v)):
                #     st = path_to_v[-1]
                #     pdb.set_trace()
                #     ind = path_to_v.index(st)
                #     path_to_v = path_to_v[ind:]
                #     break
                j += 1

            self.nmcc = path_to_v[::-1]
        else:
            self.nmcc = []


        print(self.nmcc)
        print('====')

        pdb.set_trace()

        val = self.M[self.v_star-1]
        indices = np.where(self.M==val)[0]
        c_d = {}
        for v_star in indices:
            path_to_v = []
            path_to_v = [v_star]
            next_v = v_star
            j = 0
            for i in range(glen, 0, -1):
                if i == glen:
                    loopv = next_v
                next_v = int(self.pi[glen-j,next_v-1])
                # if i == glen:
                #     loopv = next_v
                path_to_v.append(next_v)
                if i!=glen and next_v == loopv:
                    break
                if len(path_to_v) != len(set(path_to_v)):
                    st = path_to_v[-1]
                    ind = path_to_v.index(st)
                    path_to_v = path_to_v[ind:]
                    break
                if i == 1 and len(path_to_v) == len(set(path_to_v)):
                    path_to_v = [v_star] + path_to_v
                j += 1

            self.nmcc = path_to_v[::-1]
            if (len(self.nmcc)-1) != cycle_len and self.lam_star<-1e-6:
                print('hooop length is different')

            nmcc_cost = 0
            for i in range(len(self.nmcc)-1):
                u = self.nmcc[i]
                v = self.nmcc[i+1]
                try:
                    nmcc_cost += self.nxg[u][v]['mu'] + 2*G.lam*G.nxg[u][v]['flow']*self.nxg[u][v]['var']
                except:
                    nmcc_cost += self.nxg[u][v]['mu'] + 2*G.lam*G.nxg[v][u]['flow']*self.nxg[u][v]['var']

            if round(self.lam_star,5) != round(nmcc_cost/(len(self.nmcc)-1),5):
                print('hooop cost is different')

            c_d[v_star] = {'nmcc':self.nmcc, 'mean_cost':nmcc_cost/(len(self.nmcc)-1)}

        pdb.set_trace()



def get_xi():
    arc_data = {'Start':[], 'End':[], 'Flow':[], 'Capacity':[], 'Var_Cost':[], 'LowerBound':[], 'UpperBound':[]}
    node_data = {'Node':[], 'Imbalance':[]}
    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)
        arc_data['Flow'].append(e.get('flow', 0))
        arc_data['Capacity'].append(e.get('capacity', 0))
        arc_data['Var_Cost'].append(e.get('var', 0))
        arc_data['LowerBound'].append(-e.get('flow',0))
        arc_data['UpperBound'].append(e.get('capacity',0)-e.get('flow',0))

    for u, n in G.nxg.nodes(data=True):

        node_data['Node'].append(u)
        node_data['Imbalance'].append(0)

    arc_df = pandas.DataFrame(data=arc_data)
    node_df = pandas.DataFrame(data=node_data)
    sp = MinCostFlowXi(node_df, arc_df, lam=G.lam)
    results = sp.solve()
    sp.m.solutions.load_from(results)


    for u, v, e in G.nxg.edges(data=True):
        link = (u,v)
        G.nxg[u][v]['xi'] = sp.m.Y[link].value
    # pdb.set_trace()


def get_xMV(lam, G):
    arc_data = {'Start':[], 'End':[], 'Mean_Cost':[], 'Var_Cost':[], 'LowerBound':[], 'UpperBound':[]}
    node_data = {'Node':[], 'Imbalance':[]}
    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)
        arc_data['Mean_Cost'].append(e.get('mu', 0))
        arc_data['Var_Cost'].append(e.get('var',0))
        arc_data['LowerBound'].append(0)
        arc_data['UpperBound'].append(e.get('capacity',0))


    for u, n in G.nxg.nodes(data=True):
        node_data['Node'].append(u)
        node_data['Imbalance'].append(n.get('demand',0))


    arc_df = pandas.DataFrame(data=arc_data)
    node_df = pandas.DataFrame(data=node_data)
    sp = MinCostFlowMV(node_df, arc_df, lam=lam)
    results = sp.solve()
    sp.m.solutions.load_from(results)

    var_cost = 0
    mean_cost = 0
    soln = []
    for u, v, e in G.nxg.edges(data=True):
        link = (u,v)
        var_cost += sp.m.Y[link].value**2 * e.get('var',0)
        mean_cost += sp.m.Y[link].value * e.get('mu',0)
        soln.append(sp.m.Y[link].value)
    return var_cost, mean_cost, soln


import sys
sys.path.append('~/Downloads/graph-tool-2.27/src/graph_tool/')

import graph_tool as gt



def get_prop_type(value, key=None):
    """
    Performs typing and value conversion for the graph_tool PropertyMap class.
    If a key is provided, it also ensures the key is in a format that can be
    used with the PropertyMap. Returns a tuple, (type name, value, key)
    """
    # Deal with the value
    if isinstance(value, bool):
        tname = 'bool'

    elif isinstance(value, int):
        tname = 'float'
        value = float(value)

    elif isinstance(value, float):
        tname = 'float'

    elif isinstance(value, unicode):
        tname = 'string'

    elif isinstance(value, dict):
        tname = 'object'

    else:
        tname = 'string'
        value = str(value)

    return tname, value, key


def nx2gt(nxG):
    """
    Converts a networkx graph to a graph-tool graph.
    """
    # Phase 0: Create a directed or undirected graph-tool Graph
    pdb.set_trace()
    gtG = gt.Graph(directed=nxG.is_directed())

    # Add the Graph properties as "internal properties"
    for key, value in nxG.graph.items():
        # Convert the value and key into a type for graph-tool
        tname, value, key = get_prop_type(value, key)

        prop = gtG.new_graph_property(tname) # Create the PropertyMap
        gtG.graph_properties[key] = prop     # Set the PropertyMap
        gtG.graph_properties[key] = value    # Set the actual value

    # Phase 1: Add the vertex and edge property maps
    # Go through all nodes and edges and add seen properties
    # Add the node properties first
    nprops = set() # cache keys to only add properties once
    for node, data in nxG.nodes(data=True):

        # Go through all the properties if not seen and add them.
        for key, val in data.items():
            if key in nprops: continue # Skip properties already added

            # Convert the value and key into a type for graph-tool
            tname, _, key  = get_prop_type(val, key)

            prop = gtG.new_vertex_property(tname) # Create the PropertyMap
            gtG.vertex_properties[key] = prop     # Set the PropertyMap

            # Add the key to the already seen properties
            nprops.add(key)

    # Also add the node id: in NetworkX a node can be any hashable type, but
    # in graph-tool node are defined as indices. So we capture any strings
    # in a special PropertyMap called 'id' -- modify as needed!
    gtG.vertex_properties['id'] = gtG.new_vertex_property('string')

    # Add the edge properties second
    eprops = set() # cache keys to only add properties once
    for src, dst, data in nxG.edges(data=True):

        # Go through all the edge properties if not seen and add them.
        for key, val in data.items():
            if key in eprops: continue # Skip properties already added

            # Convert the value and key into a type for graph-tool
            tname, _, key = get_prop_type(val, key)

            prop = gtG.new_edge_property(tname) # Create the PropertyMap
            gtG.edge_properties[key] = prop     # Set the PropertyMap

            # Add the key to the already seen properties
            eprops.add(key)

    # Phase 2: Actually add all the nodes and vertices with their properties
    # Add the nodes
    vertices = {} # vertex mapping for tracking edges later
    for node, data in nxG.nodes(data=True):

        # Create the vertex and annotate for our edges later
        v = gtG.add_vertex()
        vertices[node] = v

        # Set the vertex properties, not forgetting the id property
        data['id'] = str(node)
        for key, value in data.items():
            gtG.vp[key][v] = value # vp is short for vertex_properties

    # Add the edges
    for src, dst, data in nxG.edges(data=True):

        # Look up the vertex structs from our vertices mapping and add edge.
        e = gtG.add_edge(vertices[src], vertices[dst])

        # Add the edge properties
        for key, value in data.items():
            gtG.ep[key][e] = value # ep is short for edge_properties

    # Done, finally!
    return gtG


 def build_res_network(self, demand='demand', capacity='capacity', weight='weight', mu='mu', var='var', flow='flow', jw=False):
        """Build a residual network and initialize a zero flow.
        """
        if sum(self.nxg.nodes[u].get(demand, 0) for u in self.nxg) != 0:
            raise nx.NetworkXUnfeasible("Sum of the demands should be 0.")

        R = MCF_DiGraph(residual=True)
        R.nxg.add_nodes_from((u, {'excess': -self.nxg.nodes[u].get(demand, 0),
                              'potential': 0}) for u in self.nxg)

        inf = float('inf')
        # Detect selfloops with infinite capacities and negative weights.
        for u, v, e in nx.selfloop_edges(self.nxg, data=True):
            if e.get(weight, 0) < 0 and e.get(capacity, inf) == inf:
                raise nx.NetworkXUnbounded(
                    'Negative cost cycle of infinite capacity found. '
                    'Min cost flow may be unbounded below.')

        # Extract edges with positive capacities. Self loops excluded.
        if self.nxg.is_multigraph():
            edge_list = [(u, v, k, e)
                         for u, v, k, e in self.nxg.edges(data=True, keys=True)
                         if u != v and e.get(capacity, inf) > 0]
        else:

            edge_list = [(u, v, e) for u, v, e in self.nxg.edges(data=True)
                         if u != v and e.get(capacity, inf) > 0]

        inf = max(sum(abs(R.nxg.nodes[u]['excess']) for u in R.nxg), 2 * sum(e[capacity] for u, v, e in edge_list if capacity in e and e[capacity] != inf)) or 1

        for u, v, e in edge_list:
            # pdb.set_trace()
            r_vu = e.get(flow, 0)
            r_uv = e.get(capacity,0) - e.get(flow,0)

            # mu_val = e.get(mu, 0)
            # var_val = e.get(var, 0)
            flow_val = e.get(flow, 0)
            cap_val = e.get(capacity, 0)
            weight_val = e.get(weight,0)



            if not jw:
                if r_uv != 0:
                    R.nxg.add_edge(u, v, capacity=cap_val, r=r_uv)

                if r_vu != 0:
                    R.nxg.add_edge(v, u, capacity=cap_val, r=r_vu)
            else:
                if r_uv != 0:
                    R.nxg.add_edge(u, v, r=r_uv)
                if r_vu != 0:
                    R.nxg.add_edge(v, u, r=r_vu)

        R.nxg.graph['inf'] = inf
        return R

def mosek_solve_var(G, lam):
    start = time.time()
    G.set_lambda(lam=lam)

    m = G.nxg.number_of_edges()
    n = G.nxg.number_of_nodes()

    mu = np.zeros(m)
    sigma = np.zeros(m)
    var = np.zeros(m)
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
        var[i] = e.get('var', 0)
        cap[i] = e.get('capacity', 0)
        arc_dict[i] = (u, v)
        rows.append(u - 1)
        cols.append(i)
        values.append(1)
        rows.append(v - 1)
        cols.append(i)
        values.append(-1)
        i += 1

    i = 0
    for node, dat in G.nxg.nodes(data=True):
        d[i] = dat['demand']
        i += 1

    a = np.asarray(mu)
    np.savetxt("mu.csv", a, delimiter=",")
    a = np.asarray(var)
    np.savetxt("var.csv", a, delimiter=",")
    a = np.asarray(cap)
    np.savetxt("cap.csv", a, delimiter=",")
    a = np.asarray(d)
    np.savetxt("d.csv", a, delimiter=",")
    a = np.asarray(rows)
    np.savetxt("Arows.csv", a, delimiter=",")
    a = np.asarray(cols)
    np.savetxt("Acols.csv", a, delimiter=",")
    a = np.asarray(values)
    np.savetxt("Avalues.csv", a, delimiter=",")
    # lam = np.array(lam)
    a = np.asarray([lam])
    np.savetxt("lam.csv", a, delimiter=",")
    with mf.Model() as model:
        pdb.set_trace()
        A = Matrix.sparse(n, m, rows, cols, values)


    #     # Create variables
    #     const_time = time.time()


    #     x = model.variable(m, Domain.greaterThan(0.0))

    #     # x_n = model.variable(m, Domain.greaterThan(0.0))

    #     # for i in range(m):
    #     #     decoy_2 = Var.vstack(x_n.index(i), x.index(i))
    #     #     model.constraint(decoy_2, Domain.inQCone())


    #     A = Matrix.sparse(n, m, rows, cols, values)
    #     model.objective("myobj",ObjectiveSense.Minimize,Expr.add(Expr.dot(mu, x), Expr.mul(Expr.dot(Expr.mulElm(x, x), var),lam)))
    #     model.constraint(x, Domain.lessThan(cap))
    #     model.constraint(Expr.sub(Expr.mul(A, x), d), Domain.equalsTo(0.0))


    #     # const_time=time.time() - const_time
    #     solve_time=time.time()

    #     model.solve()

    #     solve_time=time.time() - solve_time

    #     # print(model.primalObjValue())
    #     # print('const_time: ', const_time)
    #     # print('total_time_elapsed: ', solve_time + const_time)
    #     # print(model.getProblemStatus())

    #     cvx_obj=model.primalObjValue()
    #     i=0
    #     arc_dict={}
    #     for u, v, e in G.nxg.edges(data=True):
    #         G.nxg[u][v]['flow']=x.level()[i]
    #         i += 1
    #     # pdb.set_trace()

    #     all_elapsed = time.time()-start
    #     return all_elapsed, solve_time

def solve_var(G, lam):
    start = time.time()
    cvx_time_st=time.time()
    G.set_lambda(lam=lam)

    m=G.nxg.number_of_edges()
    n=G.nxg.number_of_nodes()

    mu=np.zeros(m)
    var=np.zeros(m)
    sigma=np.zeros(m)
    cap=np.zeros(m)
    d=np.zeros(n)
    F=np.zeros((m, n))
    x=cp.Variable(m)
    x_n=cp.Variable(m)

    i=0
    arc_dict={}
    for u, v, e in G.nxg.edges(data=True):
        mu[i]=e.get('mu', 0)
        var[i]=e.get('var', 0)
        sigma[i]=np.sqrt(e.get('var', 0))
        cap[i]=e.get('capacity', 0)
        arc_dict[i]=(u, v)
        i += 1

    i=0
    for node, dat in G.nxg.nodes(data=True):
        d[i]=dat['demand']
        i += 1

        arc_data={'Start': [], 'End': []}
    for u, v, e in G.nxg.edges(data=True):
        arc_data['Start'].append(u)
        arc_data['End'].append(v)

    constraints=[0 <= x, x <= cap]

    for i in range(n):
        preds=[]
        succs=[]
        for key, vals in arc_dict.items():
            if vals[0] == i + 1:
                preds.append(key)
            if vals[1] == i + 1:
                succs.append(key)

        constraint=sum(x[p] for p in preds) - sum(x[s]
                                                    for s in succs) == d[i]
        constraints.append(constraint)

    # Construct the problem.
    objective=cp.Minimize(mu.T * x + lam * var.T * x**2)

    prob=cp.Problem(objective, constraints)

    cvx_time_st = time.time()
    result=prob.solve(solver='MOSEK', verbose=False)
    # print("status:", prob.status)

    cvx_soln=x.value
    cvx_obj=objective.value
    cvx_elapsed=time.time() - cvx_time_st
    # print('elapsed time by cvx solve: ', cvx_elapsed)
    # print(cvx_soln)

    i=0
    arc_dict={}
    for u, v, e in G.nxg.edges(data=True):
        G.nxg[u][v]['flow']=cvx_soln[i]
        i += 1
    return time.time() - start, cvx_elapsed

def solve_mcf_sa_old(G, R, fully=False, tol=1e-2):

    prev_var_cost=-1.0
    R.DG_2_dict(G)

    while R.nmcc_exists:

        R.set_nmcc(G)
        nmcc=R.get_nmcc()

        delta=R.find_delta()
        xsi_star=R.find_xsi_star(G)

        augment_amount=min(xsi_star, delta)
        print('----Augment Amount----: ', augment_amount)

        if abs(augment_amount) < 1e-20:
            break

        R.augment_flow(augment_amount)
        prev_var_cost=G.var_cost()
        G.adjust_flow(nmcc, augment_amount)
        G.set_weights()

        R.DG_2_dict(G)

        if not fully:
            if abs(prev_var_cost - G.var_cost()) < tol:
                break

def gurobi_solve_var(G, lam, max_iters=None):

    start=time.time()
    m=G.nxg.number_of_edges()
    n=G.nxg.number_of_nodes()

    mu=np.zeros(m)
    var=np.zeros(m)
    cap=np.zeros(m)
    d=np.zeros(n)
    F=np.zeros((m, n))

    rows = []
    cols = []
    values = []
    i=0
    arc_dict={}

    for u, v, e in G.nxg.edges(data=True):
        mu[i]=e.get('mu', 0)
        var[i]=e.get('var', 0)
        cap[i]=e.get('capacity', 0)
        arc_dict[i]=(u, v)
        rows.append(u - 1)
        cols.append(i)
        values.append(1)
        rows.append(v - 1)
        cols.append(i)
        values.append(-1)
        i += 1

    i=0
    for node, dat in G.nxg.nodes(data=True):
        d[i]=dat['demand']
        i += 1

    model=Model('gurobi')

    # Create variables
    x=model.addVars(m, name="x")
    decoy=model.addVar()

    # obj
    mobj2=sum(x[k] * x[k] * var[k] for k in range(len(x)))
    # mobj=sum(x[k] * mu[k] for k in range(len(x)))
    model.setObjective((mobj2), GRB.MINIMIZE)

    # Arc capacity constraints
    model.addConstrs((x[i] <= cap[i] for i in range(m)), "upper")
    model.addConstrs((0 <= x[i] for i in range(m)), "lower")
    # model.addConstr( decoy*decoy >= (sum(x_n[k]*x_n[k] for k in range(len(x_
    # const_time=time.time()

    for i in range(n):
        preds=[]
        succs=[]
        for key, vals in arc_dict.items():
            if vals[0] == i + 1:
                preds.append(key)
            if vals[1] == i + 1:
                succs.append(key)

        # preds succs reversed
        # pdb.set_trace()
        model.addConstr(sum(x[p] for p in preds) - sum(x[s]
                                                       for s in succs) == d[i])
    # const_time=time.time() - const_time
    # Flow conservation constraints
    # model.addConstrs((flow.sum(h,'*',j) + inflow[h,j] == flow.sum(h,j,'*')
    #     for h in commodities for j in nodes), "node")

    # Compute optimal solution

    # if max_iters == None:
    #     model.Params.BarQCPConvTol = 1e-5
    #     # model.Params.BarIterLimit = 25
    #     model.Params.BarConvTol = 1e-5
    #     model.Params.FeasibilityTol = 1e-9
    # else:
    #     model.Params.BarQCPConvTol = 1
    #     model.Params.BarConvTol = 1e-5
    #     model.Params.FeasibilityTol = 1e-9
    #     model.Params.BarIterLimit = max_iters

    model.Params.OutputFlag=True
    model.update()
    cvx_time_st=time.time()
    model.optimize()

    # Print solution
    # print(model.status)

    # cvx_obj = round(model.ObjVal, cutoff)
    cvx_elapsed=time.time() - cvx_time_st
    cvx_obj=model.ObjVal


    i=0
    for u, v, e in G.nxg.edges(data=True):
        G.nxg[u][v]['flow']=model.x[i]
        i += 1
    elapsed=time.time() - start
    return elapsed, cvx_elapsed


def gurobi_solve_xi(G, lam, max_iters=None):

    start=time.time()

    m=G.nxg.number_of_edges()
    n=G.nxg.number_of_nodes()

    x_=np.zeros(m)
    mu=np.zeros(m)
    var=np.zeros(m)
    sigma=np.zeros(m)
    cap=np.zeros(m)
    d=np.zeros(n)
    F=np.zeros((m, n))
    # xi = cp.Variable(m)

    i=0
    arc_dict={}
    for u, v, e in G.nxg.edges(data=True):
        x_[i]=e.get('flow', 0)
        # if abs(x_[i]) < tol:
        #     x_[i] = 0
        mu[i]=e.get('mu', 0)
        var[i]=e.get('var', 0)
        sigma[i]=np.sqrt(e.get('var', 0))
        cap[i]=e.get('capacity', 0)
        arc_dict[i]=(u, v)
        i += 1

    i=0
    for node, dat in G.nxg.nodes(data=True):
        d[i]=dat['demand']
        i += 1

    model=Model('gurobi')
    xi=model.addVars(m, name="xi")

    mobj=sum(xi[k] * xi[k] * var[k] for k in range(len(xi)))
    mobj2=sum(var[k] * x_[k] * xi[k] for k in range(len(xi)))
    model.setObjective((lam * mobj + 2 * mobj2), GRB.MINIMIZE)

    # J = []
    K=[]
    L=[]

    for i in range(len(x_)):
        #     if abs(x_[i] - cap[i]) <= 1e-5:
        #         J.append(i)
        if abs(x_[i]) == 0:
            K.append(i)
        else:
            L.append(i)

    # pdb.set_trace()
    # model.addConstrs((-x_[i] <= xi[i] <= 0 for i in J), "J")
    model.addConstrs((0 == xi[i] for i in K))
    model.addConstrs((-x_[i] <= xi[i] <= cap[i] - x_[i] for i in L))

    # model.addConstr( decoy*decoy >= (sum(x_n[k]*x_n[k] for k in range(len(x_

    for i in range(n):
        preds=[]
        succs=[]
        for key, vals in arc_dict.items():
            if vals[0] == i + 1:
                preds.append(key)
            if vals[1] == i + 1:
                succs.append(key)

        # preds succs reversed
        # pdb.set_trace()
        model.addConstr(sum(xi[p] for p in preds) - sum(xi[s]
                                                        for s in succs) == 0)

    # for i in range(n):
    #     preds = []
    #     succs = []
    #     for key, vals in arc_dict.items():
    #         if vals[0] == i + 1:
    #             preds.append(key)
    #         if vals[1] == i + 1:
    #             succs.append(key)

    #     nodelist = [n for n in G.nxg.nodes()]
    #     node = nodelist[i]
    #     plus = 0
    #     for v in G.nxg.neighbors(node):
    #         plus += G.nxg[node][v]['flow']

    #     # pdb.set_trace()
    #     minus = 0
    #     for e in G.nxg.in_edges(node):
    #         minus += G.nxg[e[0]][e[1]]['flow']

    #     model.addConstr(sum(xi[p] for p in preds) - sum(xi[s]
    # for s in succs) == 0 + plus - minues)

        # constraint = sum(xi[p] for p in preds) - sum(xi[s]
        #                                              for s in succs) == 0 + plus - minus
        # constraints.append(constraint)

    # if max_iters == None:
    #     model.Params.BarQCPConvTol = 1e-5
    #     # model.Params.BarIterLimit = 25
    #     model.Params.BarConvTol = 1e-5
    #     model.Params.FeasibilityTol = 1e-9
    # else:
    #     model.Params.BarQCPConvTol = 1
    #     model.Params.BarConvTol = 1e-5
    #     model.Params.FeasibilityTol = 1e-9
    #     model.Params.BarIterLimit = max_iters

    model.Params.OutputFlag=True
    model.update()
    cvx_time_st=time.time()
    model.optimize()

    # Print solution
    print(model.status)
    cvx_obj=model.ObjVal
    # cvx_obj = round(model.ObjVal, cutoff)
    cvx_elapsed=time.time() - cvx_time_st

    i=0
    for u, v, e in G.nxg.edges(data=True):
        G.nxg[u][v]['xi']=model.xi[i]
        i += 1

    elapsed=time.time() - start
    return elapsed, cvx_elapsed

def solve_xi(G, lam, tol=1e-19):

    cvx_time_st=time.time()
    G.set_lambda(lam=lam)

    m=G.nxg.number_of_edges()
    n=G.nxg.number_of_nodes()

    x_=np.zeros(m)
    mu=np.zeros(m)
    var=np.zeros(m)
    sigma=np.zeros(m)
    cap=np.zeros(m)
    d=np.zeros(n)
    F=np.zeros((m, n))
    xi=cp.Variable(m)

    i=0
    arc_dict={}
    for u, v, e in G.nxg.edges(data=True):
        x_[i]=e.get('flow', 0)
        # if abs(x_[i]) < tol:
        #     x_[i] = 0
        mu[i]=e.get('mu', 0)
        var[i]=e.get('var', 0)
        sigma[i]=np.sqrt(e.get('var', 0))
        cap[i]=e.get('capacity', 0) + x_[i]
        # cap[i] = e.get('capacity', 0)

        arc_dict[i]=(u, v)
        i += 1

    xi=xi + x_

    i=0
    for node, dat in G.nxg.nodes(data=True):
        d[i]=dat['demand']
        i += 1
        # arc_data = {'Start': [], 'End': []}

    # for u, v, e in G.nxg.edges(data=True):
    #     arc_data['Start'].append(u)
    #     arc_data['End'].append(v)

    # J = []
    # K = []
    # L = []

    # for i in range(len(x_)):
    #     if abs(x_[i] - cap[i]) <= 1e-5:
    #         J.append(i)
    #     elif abs(x_[i]) <= 1e-5:
    #         K.append(i)
    #     else:
    #         L.append(i)

    # constraints = []
    # constraints = [-x_[i] <= xi[i] for i in J]
    # constraints = [xi[i] <= 0 for i in J]
    # constraints = [0 == xi[i] for i in K]
    # constraints = [-x_[i] <= xi[i] for i in L]
    # constraints = [xi[i] <= cap[i] - x_[i] for i in L]

    # model.addConstrs((-x_[i] <= xi[i] <= 0 for i in J), "J")
    # model.addConstrs((0 == xi[i] for i in K), "K")
    # model.addConstrs((-x_[i] <= xi[i] <= cap[i] - x_[i] for i in L), "L")

    constraints=[]
    constraints=[0 <= xi, xi <= cap]

    for i in range(n):
        preds=[]
        succs=[]
        for key, vals in arc_dict.items():
            if vals[0] == i + 1:
                preds.append(key)
            if vals[1] == i + 1:
                succs.append(key)

        nodelist=[n for n in G.nxg.nodes()]
        node=nodelist[i]
        plus=0
        for v in G.nxg.neighbors(node):
            plus += G.nxg[node][v]['flow']

        minus=0
        for e in G.nxg.in_edges(node):
            minus += G.nxg[e[0]][e[1]]['flow']

        # constraint = sum(xi[p] for p in preds) - sum(xi[s]
        #                                              for s in succs) == 0
        constraint=sum(xi[p] for p in preds) - sum(xi[s]
                                                     for s in succs) == 0 + plus - minus
        constraints.append(constraint)

    # Construct the problem.
    objective=cp.Minimize(G.lam * var.T * xi**2 +
                            2 * np.multiply(var, x_).T * xi)

    prob=cp.Problem(objective, constraints)

    result=prob.solve(solver='GUROBI', verbose=False)
    # pdb.set_trace()
    # print("status:", prob.status)

    cvx_soln=xi.value
    cvx_obj=objective.value
    elapsed=time.time() - cvx_time_st
    cvx_elapsed=prob.solver_stats.solve_time

    i=0
    arc_dict={}
    for u, v, e in G.nxg.edges(data=True):
        G.nxg[u][v]['xi']=cvx_soln[i]
        i += 1
    return elapsed, cvx_elapsed


def bisect_solve_alt(lam_bar, seed, num_nodes, mid, tol=1e-2):
    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()
    # pdb.set_trace()
    start=time.time()
    disc_tot=0
    nmcc_time_tot=0
    scc_time_tot=0

    # G.set_lambda(lam=1e-5)
    # discount, nmcc_time, scc_time =solve_mcf_sa(G, R, tol=1e-1)
    # disc_tot = discount
    # nmcc_time_tot = nmcc_time
    # scc_time_tot = scc_time
    # sigma_cost = sqrt(G.var_cost())
    # low = lam_bar/(2*sigma_cost)
    sigma_lam=get_sigma_cost(G, 0)
    low=lam_bar / (2 * sigma_lam) - tol

    # sigma_lam = get_sigma_cost(G,low)

    # G.set_lambda(low)
    # discount, nmcc_time, scc_time = solve_mcf_sa(G, R, fully=False, tol=1e-10)
    # disc_tot += discount
    # nmcc_time_tot += nmcc_time
    # scc_time_tot += scc_time
    # sigma_lam = sqrt(G.var_cost())

    # f_low = low - lam_bar/float(2*sigma_lam)

    found=False
    iters=0

    # G.set_lambda(lam=1e9)
    # discount, nmcc_time, scc_time =solve_mcf_sa(G, R, tol=1e-1)
    # disc_tot += discount
    # nmcc_time_tot += nmcc_time
    # scc_time_tot += scc_time
    # sigma_cost = sqrt(G.var_cost())
    # high = lam_bar/(2*sigma_cost)
    sigma_lam=get_sigma_cost(G, 1e6)
    high=lam_bar / (2 * sigma_lam) + tol

    sigma_lam=get_sigma_cost(G, high)

    # G.set_lambda(high)
    # discount, nmcc_time, scc_time = solve_mcf_sa(G, R, fully=False, tol=1e-10)
    # disc_tot += discount
    # nmcc_time_tot += nmcc_time
    # scc_time_tot += scc_time
    # sigma_lam = sqrt(G.var_cost())

    # print('low, high, mid:', low, high, mid)
    f_high=high - float(lam_bar) / float((2 * sigma_lam))
    # print(f_low, f_high)

    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()

    disc_tot=0
    nmcc_time_tot=0
    scc_time_tot=0
    mid_prev=1e6

    while not found:
        iters += 1
        mid=(high + low) / 2.0
        # print('BSlam: ', mid)

        G.set_lambda(lam=mid)
        if abs(mid_prev - mid) < tol:
            found=True
            discount, nmcc_time, scc_time=solve_mcf_sa(G, R, fully=True)
            # discount, nmcc_time, scc_time =solve_mcf_sa(G, R, tol=1e-20)

            disc_tot += discount
            nmcc_time_tot += nmcc_time
            scc_time_tot += scc_time
            break
        else:
            # discount, nmcc_time, scc_time =solve_mcf_sa(G, R, tol=1e-1)
            discount, nmcc_time, scc_time=solve_mcf_sa(G, R, fully=True)

            # discount, nmcc_time, scc_time =solve_mcf_sa(G, R)

            disc_tot += discount
            nmcc_time_tot += nmcc_time
            scc_time_tot += scc_time

        var_cost=G.var_cost()

        sigma_lam=np.sqrt(var_cost)
        f_mid=mid - float(lam_bar) / float((2 * sigma_lam))

        if np.sign(f_mid) == np.sign(f_high):
            high=mid
            f_high=f_mid
        else:
            low=mid

        mid_prev=mid

    # print('FINAL LAM: ', mid)
    end=time.time()
    var_cost=G.var_cost()
    algo_obj=lam_bar * np.sqrt(var_cost) + G.mu_cost()
    algo_elapsed=end - start - disc_tot + nmcc_time_tot + scc_time_tot
    algo_soln=G.get_soln()
    return algo_obj, algo_elapsed, algo_soln, iters

def NR_solve_alt(lam_bar, seed, num_nodes, cvx_obj=None, tol=1e-1):
    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()

    start=time.time()
    disc_tot=0
    nmcc_time_tot=0
    scc_time_tot=0

    found=False
    iters=0

    sigma_lam=get_sigma_cost(G, 1e6)
    high=lam_bar / (2 * sigma_lam) + tol

    # sigma_lam = get_sigma_cost(G,high)

    discount, nmcc_time, scc_time=solve_mcf_sa(G, R, fully=False, tol=1e-1)
    var_cost=G.var_cost()
    sigma_lam=np.sqrt(var_cost)

    # var_cost = sigma_lam**2
    # G.set_lambda(high)
    # discount, nmcc_time, scc_time = solve_mcf_sa(G, R, fully=False,
    # tol=1e-10)
    disc_tot += discount
    nmcc_time_tot += nmcc_time
    scc_time_tot += scc_time

    iters=0
    deduct=0
    curios_time=0
    found=False
    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()

    lam_prev=1e6
    disc_me_tot=0
    lam=high
    while not found:

        # print('NRlam: ', lam)
        iters += 1
        G.set_lambda(lam=lam)
        # print('----')
        # print(lam)

        if iters != 1:
            if abs(lam_prev - lam) < tol:
                found=True
                # discount, nmcc_time, scc_time = solve_mcf_sa(G, R,
                # fully=True)
                discount, nmcc_time, scc_time=solve_mcf_sa(G, R, tol=1e-20)

                disc_tot += discount
                nmcc_time_tot += nmcc_time
                scc_time_tot += scc_time
                break
            else:
                discount, nmcc_time, scc_time=solve_mcf_sa(
                    G, R, fully=False, tol=1e-1)
                # discount, nmcc_time, scc_time = solve_mcf_sa(G, R,
                # fully=True)

                # discount, nmcc_time, scc_time = solve_mcf_sa(G, R)

                disc_tot += discount
                nmcc_time_tot += nmcc_time
                scc_time_tot += scc_time
                var_cost=G.var_cost()
                sigma_lam=np.sqrt(var_cost)

        f_lam=lam - lam_bar / (2 * sigma_lam)

        # start_curios = time.time()
        # solve_xi(G, lam)
        # var_cost_der = G.xi_cost()
        # curios_time +=time.time() - start_curios
        #######################
        disc_me=time.time()
        G.find_feasible_flow_xi()
        R_xi=G.build_res_network_xi()
        disc_me_tot += time.time() - disc_me
        # discount, nmcc_time, scc_time = solve_xi_mmc(G, R_xi, fully=False,
        # tol=1e-1)
        discount, nmcc_time, scc_time=solve_xi_mmc(G, R_xi)

        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time
        var_cost_der=G.xi_cost()
        #######################
        f_lam_der=1 + (lam_bar / 2) * (var_cost**(-3 / 2)) * var_cost_der

        if np.sign(f_lam_der) == -1:
            print('hopppala')
            pdb.set_trace()

        lam_prev=lam
        lam=lam_prev - f_lam / f_lam_der

        # algo_obj = lam_bar*np.sqrt(var_cost) + G.mu_cost()
    # print('final lam: ', lam)
    end=time.time()
    var_cost=G.var_cost()
    algo_obj=lam_bar * np.sqrt(var_cost) + G.mu_cost()
    algo_elapsed=end - start - disc_tot - \
        disc_me_tot + nmcc_time_tot + scc_time_tot
    algo_soln=G.get_soln()
    return algo_obj, algo_elapsed, algo_soln, iters, curios_time

def cheatbsbigint(lam_bar, seed, num_nodes, tol=1e-3):
    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()

    start=time.time()

    low=0
    high=lam_bar

    found=False
    iters=0
    # pdb.set_trace()
    # start = time.time()
    G.set_lambda(lam=high)
    solve_var(G, high)
    # print(time.time()-start)
    # G = create_random_graph(seed=seed, num_nodes=num_nodes)
    # R = G.build_res_network()
    # G.set_lambda(lam=high)
    # start = time.time()
    # discount, nmcc_time, scc_time =solve_mcf_sa(G, R, fully=True, tol=1e-20)
    # print(time.time() - start)
    # pdb.set_trace()
    var_cost=G.var_cost()
    sigma_lam=np.sqrt(var_cost)
    f_high=high - float(lam_bar) / float((2 * sigma_lam))
    mid_prev=1e6

    while not found:
        iters += 1
        mid=(high + low) / 2.0

        G.set_lambda(lam=mid)
        solve_var(G, mid)

        var_cost=G.var_cost()
        sigma_lam=np.sqrt(var_cost)
        f_mid=mid - float(lam_bar) / float((2 * sigma_lam))

        if np.sign(f_mid) == np.sign(f_high):
            high=mid
            f_high=f_mid
        else:
            low=mid

        mid_prev=mid
        if abs(f_mid) == 0:
            break
        print(f_mid)

    # print('FINAL LAM: ', mid)
    end=time.time()
    var_cost=G.var_cost()
    algo_obj=lam_bar * np.sqrt(var_cost) + G.mu_cost()
    algo_elapsed=end - start
    algo_soln=G.get_soln()
    return algo_obj, algo_elapsed, iters

def cheatbssmallint(lam_bar, seed, num_nodes, tol=1e-3):
    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()

    start=time.time()

    sigma_lam=get_sigma_cost(G, 0)
    low=lam_bar / (2 * sigma_lam) - tol

    found=False
    iters=0

    tol=1e-2

    sigma_lam=get_sigma_cost(G, 1e6)
    high=lam_bar / (2 * sigma_lam) + tol

    sigma_lam=get_sigma_cost(G, high)

    f_high=high - float(lam_bar) / float((2 * sigma_lam))

    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()

    mid_prev=1e6

    while not found:
        iters += 1
        mid=(high + low) / 2.0
        # print('BSlam: ', mid)

        G.set_lambda(lam=mid)
        solve_var(G, mid)
        if abs(mid_prev - mid) < tol:
            found=True

        var_cost=G.var_cost()

        sigma_lam=np.sqrt(var_cost)
        f_mid=mid - float(lam_bar) / float((2 * sigma_lam))

        if np.sign(f_mid) == np.sign(f_high):
            high=mid
            f_high=f_mid
        else:
            low=mid

        mid_prev=mid

    # print('FINAL LAM: ', mid)
    end=time.time()
    var_cost=G.var_cost()
    algo_obj=lam_bar * np.sqrt(var_cost) + G.mu_cost()
    algo_elapsed=end - start
    algo_soln=G.get_soln()
    return algo_obj, algo_elapsed, iters

def cheat_solve_nrrrrr(lam_bar, seed, num_nodes, tol=1e-5):
    G=create_random_graph(seed=seed, num_nodes=num_nodes)
    R=G.build_res_network()
    start=time.time()
    diff=1000
    first_guess=lam_bar / 2
    lam=first_guess
    iters=0

    while abs(diff) > tol:
        G.set_lambda(lam=lam)
        solve_var(G, lam)
        var_cost=G.var_cost()

        f_lam=lam - lam_bar / (2 * np.sqrt(var_cost))

        solve_xi(G, lam)
        var_cost_der=G.xi_cost()

        f_lam_der=1 + (lam_bar / 2) * (var_cost**(-3 / 2)) * var_cost_der
        print(np.sign(f_lam_der))
        lam_prev=lam
        lam=lam_prev - f_lam / f_lam_der

        iters += 1
        diff=f_lam

    # G.set_lambda(lam=lam)
    # solve_var(G, lam)
    print(lam)
    end=time.time()
    var_cost=G.var_cost()
    algo_obj=lam_bar * np.sqrt(var_cost) + G.mu_cost()
    algo_elapsed=end - start
    # algo_soln = nx.get_edge_attributes(G.nxg,'flow')

    return algo_obj, algo_elapsed, iters



def cheat_solve_nr(**kwargs):
    if kwargs['G'] != None:
        G=copy.deepcopy(kwargs['G'])
    else:
        G=create_random_graph(**kwargs)
    if kwargs['R'] != None:
        R=copy.deepcopy(kwargs['R'])
    else:
        R=G.build_res_network()

    lam_bar=kwargs['lam_bar']
    stop_tol=kwargs['stop_tol']
    subproblem_tol=kwargs['subproblem_tol']
    precision_start_tol=kwargs['precision_start_tol']

    disc_me_tot=0

    # lam_high, f_lam, lam_bound_time = get_upper_bound(G, R, lam_bar)
    # varcost = G.var_cost()
    lam=lam_bar
    start=time.time()

    # lam = lam_high
    iters=0
    deduct=0
    found=False
    disc_tot=0
    nmcc_time_tot=0
    scc_time_tot=0
    lam_prev=1e6
    subproblem_times=[]
    gap=[]
    gap_perc=[]
    times=[]
    elapsed=0
    cvx_elapsed=0
    xicost=1e6
    while not found:

        iters += 1
        G.set_lambda(lam=lam)

        if iters != 1:
            sub_start=time.time()

            if abs(f_lam) < precision_start_tol:
                elapsed_, cvx_elapsed_=gurobi_solve_var(G, lam=lam)
            else:
                elapsed_, cvx_elapsed_=gurobi_solve_var(G, lam=lam)

            elapsed += elapsed_
            cvx_elapsed += cvx_elapsed_

        varcost=G.var_cost()
        sigmacost=np.sqrt(varcost)
        f_lam=lam - lam_bar / (2 * sigmacost)

        disc_me=time.time()
        G.set_lambda(lam=lam, xi=True)
        if iters == 1:
            G.find_feasible_flow_xi()
            R_xi=G.build_res_network_xi()
        disc_tot += time.time() - disc_me
        discount, nmcc_time, scc_time, xicost=solve_xi_mmc(
            G, R_xi, xicost, fully=False)
        disc_tot += discount
        nmcc_time_tot += nmcc_time
        scc_time_tot += scc_time

        # G.set_lambda(lam=lam, xi=True)
        # elapsed_, cvx_elaped_ = solve_xi(G, lam=lam)
        # elapsed += elapsed_
        # cvx_elapsed += cvx_elaped_
        # xicost = G.xi_cost()
        f_lam_der=1 + (lam_bar / 2) * (varcost**(-3 / 2)) * xicost

        lam_prev=lam
        lam=lam_prev - f_lam / f_lam_der

        disc_me=time.time()
        G.set_lambda(lam=lam)
        varcost=G.var_cost()
        sigmacost=np.sqrt(varcost)

        algo_obj=round((lam_bar * sigmacost + G.mu_cost() *
                          kwargs['mu_scalar']) * kwargs['obj_scalar'], kwargs['cutoff'])
        disc_tot += disc_me - time.time()
        gap.append(abs(algo_obj - params['cvx_obj']))
        gap_perc.append(abs(1 - algo_obj / params['cvx_obj']))
        times.append(time.time() - start - elapsed + cvx_elapsed -
                     disc_tot + nmcc_time_tot - scc_time_tot)
        if algo_obj - params['cvx_obj'] == 0:
            found=True

    end=time.time()
    # var_cost = G.var_cost()
    # algo_obj = round((lam_bar * np.sqrt(var_cost) + G.mu_cost() *
    #           kwargs['mu_scalar']) * kwargs['obj_scalar'],kwargs['cutoff'])
    # algo_elapsed = end - start - disc_tot - disc_me_tot + \
    #     nmcc_time_tot - scc_time_tot + lam_bound_time
    algo_elapsed=end - start - elapsed + cvx_elapsed - \
        disc_tot + nmcc_time_tot - scc_time_tot
    # algo_soln = G.get_soln()
    return algo_obj, algo_elapsed, iters, np.array(subproblem_times), gap, gap_perc, times


def cheat_solve_nrr(**kwargs):
    lam_bar=kwargs['lam_bar']
    G=create_random_graph(**kwargs)
    R=G.build_res_network()
    diff=1000
    lam_high, f_lam, lam_bound_time=get_upper_bound(G, R, lam_bar)
    start=time.time()

    lam=lam_high
    iters=0
    elapsed=0
    real_elapsed=0
    while abs(diff) > kwargs['stop_tol']:
        G.set_lambda(lam=lam)
        elapsed_, real_elapsed_=gurobi_solve_var(G, lam=mid, max_iters=None)
        elapsed += elapsed_
        real_elapsed += real_elapsed_
        var_cost=G.var_cost()

        f_lam=lam - lam_bar / (2 * np.sqrt(var_cost))

        solve_xi(G, lam)
        var_cost_der=G.xi_cost()

        f_lam_der=1 + (lam_bar / 2) * (var_cost**(-3 / 2)) * var_cost_der
        lam_prev=lam
        lam=lam_prev - f_lam / f_lam_der

        iters += 1
        diff=f_lam
        print(diff)

    end=time.time()

    var_cost=G.var_cost()
    algo_obj=(lam_bar * np.sqrt(var_cost) + G.mu_cost() *
                kwargs['mu_scalar']) * kwargs['obj_scalar']
    algo_elapsed=end - start - elapsed + real_elapsed

    return algo_obj, algo_elapsed, iters


def cheat_solve_bs(**kwargs):
    if kwargs['G'] != None:
        G=copy.deepcopy(kwargs['G'])
    else:
        G=create_random_graph(**kwargs)
    if kwargs['R'] != None:
        R=copy.deepcopy(kwargs['R'])
    else:
        R=G.build_res_network()
    lam_bar=kwargs['lam_bar']
    stop_tol=kwargs['stop_tol']
    subproblem_tol=kwargs['subproblem_tol']
    precision_start_tol=kwargs['precision_start_tol']

    low=0
    found=False
    iters=0

    ####

    # tol = 1e-5

    # G.set_lambda(lam=1e-5)
    # elapsed, real_elapsed = gurobi_solve_var(G, lam=1e-5, max_iters=None)

    # # sigma_lam = get_sigma_cost(G, 1e-5)
    # sigma_lam = np.sqrt(G.var_cost())
    # low = lam_bar / (2 * sigma_lam) - tol

    # G.set_lambda(lam=1e4)
    # elapsed, real_elapsed = gurobi_solve_var(G, lam=1e4, max_iters=None)

    # # sigma_lam = get_sigma_cost(G, 1e4)
    # sigma_lam = np.sqrt(G.var_cost())
    # high = lam_bar / (2 * sigma_lam) + tol

    # G.set_lambda(lam=high)
    # elapsed, real_elapsed = gurobi_solve_var(G, lam=high, max_iters=None)

    # f_high = high - float(lam_bar) / float((2 * sigma_lam))

    #####

    # lam_high, f_high, lam_bound_time = get_upper_bound(G, R, lam_bar)
    # high = lam_high

    #####

    high=lam_bar / 2

    #####

    start=time.time()
    f_mid=1000
    disc_tot=0
    nmcc_time_tot=0
    scc_time_tot=0

    subproblem_times=[]
    mid_prev=1e6
    gap=[]
    times=[]
    gap_perc=[]
    algotimes=[]
    elapsed=0
    real_elapsed=0
    discount_me = 0
    while not found:
        iters += 1
        mid=(high + low) / 2.0

        sub_start=time.time()
        G.set_lambda(lam=mid)
        elapsed_, real_elapsed_=mosek_solve_var(G, lam=mid)
        elapsed += elapsed_
        real_elapsed += real_elapsed_

        subproblem_times.append(
            time.time() - sub_start - elapsed + real_elapsed)

        var_cost=G.var_cost()
        sigma_lam=np.sqrt(var_cost)
        f_mid=mid - float(lam_bar) / float((2 * sigma_lam))

        # if abs(f_mid) < stop_tol:
        #     found = True

        if np.sign(f_mid) == np.sign(1):
            high=mid
            # f_high = f_mid
        else:
            low=mid

        disc_me=time.time()
        mid_prev=mid
        # pdb.set_trace()
        algo_obj=round((lam_bar * sigma_lam + G.mu_cost() *kwargs['mu_scalar']) * kwargs['obj_scalar'], kwargs['cutoff'])
        print(algo_obj)
        print(mid)

        gap.append(abs(algo_obj - params['cvx_obj']))
        gap_perc.append(abs(1 - algo_obj / params['cvx_obj']))
        discount_me += time.time() - disc_me
        times.append(time.time() - start - elapsed + real_elapsed - discount_me)

        if algo_obj - params['cvx_obj'] == 0:
            found=True

    algo_elapsed=times[-1]
    return algo_obj, algo_elapsed, iters, np.array(subproblem_times), gap, gap_perc, times


# import progressbar
# def find_stop_tol():
#     stopTolList = [1e-6, 1e-12, 1e-21]
#     muTopList = [10, 30]
#     varTopList = [10, 30]
#     dTopList = [10, 200]
#     reallamBarList = np.array([1e0, 1e2])
#     numNodesList = [100, 500]
#     scale = reallamBarList.max() * 1e1
#     seeds = 99*np.arange(5)
#     precStartTolList = [1e-3, 1e-1, 1e1]
#     spTolList = [1e-1, 1e-12, 1e-21]

#     paramPerms = list(itertools.product(muTopList, varTopList,
# dTopList, reallamBarList, numNodesList, stopTolList, precStartTolList,
# spTolList))

#     column_list = ['mu_top', 'var_top', 'd_top', 'lam_bar', 'seed', 'num_nodes', 'stopTol', 'subproblem_tol', 'precision_start_tol',
#                    'cvx_obj', 'nr_obj', 'bs_obj', 'nr_gap', 'bs_gap', 'nr_iter', 'bs_iter', 'nr_elapsed', 'bs_elapsed', 'cvx_elapsed', 'nrSPtimes', 'bsSPtimes', 'nrSPtime_mean', 'bsSPtime_mean', 'nrSPtime_var', 'bsSPtime_var']
#     exp_df = pandas.DataFrame(columns=column_list)

#     bar = progressbar.ProgressBar(widgets=widgets, max_value=len(paramPerms))

#     bar.start()
#     barcount = 0

#     for perm in paramPerms:
#         barcount += 1
#         params = {}
#         params['mu_top'] = perm[0]
#         params['var_top'] = perm[1]
#         params['d_top'] = perm[2]
#         params['lam_bar'] = perm[3] / scale
#         params['num_nodes'] = perm[4]
#         params['stop_tol'] = perm[5]
#         params['precision_start_tol'] = perm[6]
#         params['subproblem_tol'] = perm[7]
#         params['cvx_obj'] = None
#         params['mu_scalar'] = 1 / scale
#         params['obj_scalar'] = scale

#         # print('==========')
#         # print(perm)
#         # print('==========')

#         for i in range(5):
#             seed = seeds[i]
#             params['seed'] = seed

#             cvx_obj, cvx_elapsed, cvx_soln = cvx_solve(**params)
#             params['cvx_obj'] = cvx_obj

#             # print(tabulate([[params['num_nodes'], params['seed'], params['lam_bar'], params[
#             #       'stop_tol']]], headers=['num_nodes', 'seed', 'lam_bar', 'stop_tol']))

#             nr_obj, nr_elapsed, nr_soln, nr_iters, curios_time, nr_subproblem_times = NR_solve(
#                 **params)
#             bs_obj, bs_elapsed, bs_soln, bs_iters, bs_subproblem_times = bisect_solve(
#                 **params)
#             nr_gap = abs(1 - nr_obj / cvx_obj)
#             bs_gap = abs(1 - bs_obj / cvx_obj)
#             # print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters, bs_gap], [
#                   # 'CVX', cvx_obj, cvx_elapsed, 1]], headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))
#             row = pandas.DataFrame([[params['mu_top'], params['var_top'], params['d_top'], perm[3], params['seed'], params['num_nodes'], params['stop_tol'], params['subproblem_tol'], params['precision_start_tol'],
#                                      cvx_obj, nr_obj, bs_obj, nr_gap, bs_gap, nr_iters, bs_iters, nr_elapsed, bs_elapsed, cvx_elapsed, nr_subproblem_times, bs_subproblem_times, np.mean(nr_subproblem_times), np.mean(bs_subproblem_times), np.var(nr_subproblem_times), np.var(bs_subproblem_times)]], columns=column_list)
#             exp_df = exp_df.append(row)
#             exp_df.to_pickle("stop_tol_experiment.pkl")
#         bar.update(barcount)
#     bar.finish()
#     return exp_df


def print_experiment():
    # print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS',
    # bs_obj, bs_elapsed, bs_iters,bs_gap], ['CVX',cvx_obj, cvx_elapsed, 1]],
    # headers=['Method','Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))
    pass

# experiments = find_stop_tol()
# experiments.to_pickle("stop_tol_experiment.pkl")
# pdb.set_trace()
# print('finished')




  # nr_gap = abs(1 - nr_obj / cvx_obj)
        # bs_gap = abs(1 - bs_obj / cvx_obj)

        # print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters, bs_gap], [
        #     'CVX', cvx_obj, cvx_elapsed, 1]], headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))

        # # print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters, bs_gap]], headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))

        # nr_iters_l.append(nr_iters)
        # bs_iters_l.append(bs_iters)
        # nr_elapsed_l.append(nr_elapsed)
        # bs_elapsed_l.append(bs_elapsed)
        # nr_gaps.append(nr_gap)
        # bs_gaps.append(bs_gap)


# pdb.set_trace()
# try:
#   times_cvx = np.mean(times_cvx_list)
#   times_bs = np.mean(times_bs_list)
#   times_nr = np.mean(times_nr_list)
#   gap_cvx = np.mean(gap_cvx_list)
#   gap_bs = np.mean(gap_bs_list)
#   gap_nr = np.mean(gap_nr_list)
# except:
#   pdb.set_trace()

# pdb.set_trace()
# plt.plot(times_cvx, gap_cvx, label='GUROBI', marker='D')
# plt.plot(times_bs, gap_bs, label='Bisection', marker='P')
# plt.plot(times_nr, gap_nr, label='Newton', marker='X')
# plt.yscale('log')
# plt.title('log')
# plt.legend(loc='upper right')
# plt.xlabel('Time (s) - Wall Clock')
# plt.ylabel('Log(Gap percentage)')
# plt.show()


#   G = create_random_graph(**params)
#   nnodes = len(G.nxg.nodes())
#   nedges = len(G.nxg.edges())
#   print('=========================================')
#   print('| #Nodes : ', nnodes, '   |   ', ' #Edges : ', nedges, ' |')
#   print('=========================================')
#   print(tabulate([['NR', nr_obj, np.mean(nr_elapsed_l), np.mean(nr_iters_l), np.mean(nr_gaps)], ['BS', bs_obj, np.mean(bs_elapsed_l), np.mean(
# bs_iters_l), np.mean(bs_gaps)], ['CVX', cvx_obj, cvx_elapsed, 1, 0]],
# headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))


# nr_gap = abs(1 - nr_obj / cvx_obj)
# bs_gap = abs(1 - bs_obj / cvx_obj)

# print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters, bs_gap], [
#     'CVX', cvx_obj, cvx_elapsed, 1]], headers=['Method', 'Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))


# for i  in range(15):
#     num_nodes = random.randint(5000, 7500)
#     seed = random.randint(1,10000)
#     lam_bar = np.random.uniform(0.1, 1)
#     print(num_nodes, seed, lam_bar)
#     print('c')
#     cvx_obj, cvx_elapsed, cvx_soln = cvx_solve(lam_bar, seed, num_nodes)
#     print('n')
#     nr_obj, nr_elapsed, nr_soln, nr_iters, curios_time, nr_subproblem_times= NR_solve(lam_bar, seed, num_nodes, cvx_obj=cvx_obj)
#     print('b')
#     bs_obj, bs_elapsed, bs_soln, bs_iters, bs_subproblem_times = bisect_solve(lam_bar, seed, num_nodes, cvx_obj=cvx_obj)
#     nr_gap=abs(1-nr_obj/cvx_obj)
#     bs_gap=abs(1-bs_obj/cvx_obj)
#     print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters,bs_gap], ['CVX',cvx_obj, cvx_elapsed, 1]], headers=['Method','Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))
#     pdb.set_trace()


# nodelist = [100, 500, 1000, 5000, 10000]
# for n in nodelist:
#     num_nodes = int(n)
#     nr_iters_l = []
#     bs_iters_l = []
#     nr_elapsed_l = []
#     bs_elapsed_l = []
#     nr_gaps=[]
#     bs_gaps=[]
#     curious_time_l = []
#     lam_ll = []
#     print(n)
#     print('=============')
#     # seed_l = np.arange(10)*10
#     # lam_bar_l = [0.14, 0.66, 0.59, 0.24, 0.16, 0.012, 0.35, 0.55, 0.445, 0.894]
#     # lam_bar_l = [20, 100, 250, 10, 25, 50, 80, 560, 1000, 2]
#     for j in range(5):
#         seed = random.randint(1,10000)
#         lam_bar = np.random.uniform(0.001, 10)
#         # lam_bar = 10
#         # seed = seed_l[j]
#         # lam_bar = lam_bar_l[j]
#         print('seed: ' + str(seed) + ' | lam_bar: ' + str(lam_bar))
#         lam_l = []
#         lam_l.append(lam_bar)
#         # nr_obj, nr_elapsed, nr_soln, nr_iters, curios_time = NR_solve(lam_bar, seed, num_nodes, cvx_obj=None)
#         # robj, relapsed, rsoln, riters, mid = bisect_solve(lam_bar, seed, num_nodes)
#         mid = 0
#         # bs_obj, bs_elapsed, bs_soln, bs_iters = bisect_solve_alt(lam_bar, seed, num_nodes, mid)
# nr_obj, nr_elapsed, nr_soln, nr_iters, curios_time =
# NR_solve_alt(lam_bar, seed, num_nodes)

#         cvx_obj, cvx_elapsed, cvx_soln = cvx_solve(lam_bar, seed, num_nodes)
#         nr_gap=abs(1-nr_obj/cvx_obj)
#         # bs_gap=abs(1-bs_obj/cvx_obj)
#         # r_gap=abs(1-robj/cvx_obj)
#         # print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters, bs_gap], ['CVX',cvx_obj, cvx_elapsed, 1]], headers=['Method','Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))
#         print(tabulate([['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['CVX',cvx_obj, cvx_elapsed, 1]], headers=['Method','Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))
#         # pdb.set_trace()
#         # print(tabulate([['RB', robj, relapsed, riters, r_gap], ['NR', nr_obj, nr_elapsed, nr_iters, nr_gap], ['BS', bs_obj, bs_elapsed, bs_iters, bs_gap], ['CVX',cvx_obj, cvx_elapsed, 1]], headers=['Method','Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))

#         nr_iters_l.append(nr_iters)
#         # bs_iters_l.append(bs_iters)
#         nr_elapsed_l.append(nr_elapsed)
#         # bs_elapsed_l.append(bs_elapsed)
#         nr_gaps.append(nr_gap)
#         # bs_gaps.append(bs_gap)

#     lam_ll.append(lam_l)
#     G = create_random_graph(seed=seed, num_nodes=num_nodes)
#     nnodes = len(G.nxg.nodes())
#     nedges = len(G.nxg.edges())
#     print('=========================================')
#     print('| #Nodes : ',nnodes ,'   |   ', ' #Edges : ',nedges, ' |')
#     print('=========================================')
#     # print(tabulate([['NR', nr_obj, np.mean(nr_elapsed_l), np.mean(nr_iters_l), np.mean(nr_gaps)], ['BS', bs_obj, np.mean(bs_elapsed_l), np.mean(bs_iters_l), np.mean(bs_gaps)], ['CVX',cvx_obj, cvx_elapsed, 1, 0]], headers=['Method','Obj_Val', 'Solve_Time', 'Iterations', 'Obj_Gap']))
# print(tabulate([['NR', nr_obj, np.mean(nr_elapsed_l),
# np.mean(nr_iters_l), np.mean(nr_gaps)], ['CVX',cvx_obj, cvx_elapsed, 1,
# 0]], headers=['Method','Obj_Val', 'Solve_Time', 'Iterations',
# 'Obj_Gap']))

        # def find_tind(t, tlist, gap):
        #   for x in tlist:
        #       if x == tlist[0] and x-t >0:
        #           # val = min(tlist, key=lambda x: abs(x-t))
        #           # tlist = list(tlist)
        #           # idx = tlist.index(val)
        #           return 1e5
        #       elif x-t > 0:
        #           break
        #       else:
        #           val = x
        #   if t > tlist[-1]:
        #       return None

        #   tlist = list(tlist)
        #   idx = tlist.index(val)
        #   return gap[idx]

        # try:
        #   taxis = np.linspace(0.25, times_cvx[-1], num=20)
        #   nr = []
        #   bs = []
        #   cvx = []
        #   for t in taxis:
        #       nr.append(find_tind(t, times_nr, gap_nr))
        #       bs.append(find_tind(t, times_bs, gap_bs))
        #       cvx.append(find_tind(t, times_cvx, gap_cvx))

        # except:
        #   pdb.set_trace()



               ##
        # params['max_iters'] = 5
        # obj, elapsed = gurobi_solve(**params)
        # gap_cvx.append(abs(obj - cvx_obj))
        # gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        # times_cvx.append(elapsed)
        # params['max_iters'] = 10
        # obj, elapsed = gurobi_solve(**params)
        # gap_cvx.append(abs(obj - cvx_obj))
        # gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        # times_cvx.append(elapsed)
        # params['max_iters'] = 15
        # obj, elapsed = gurobi_solve(**params)
        # gap_cvx.append(abs(obj - cvx_obj))
        # gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        # times_cvx.append(elapsed)
        # if gap_cvx[-1] != 0:
        #     params['max_iters'] = 25
        #     obj, elapsed = gurobi_solve(**params)
        #     gap_cvx.append(abs(obj - cvx_obj))
        #     gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        #     times_cvx.append(elapsed)
        # if gap_cvx[-1] != 0:
        #     params['max_iters'] = 50
        #     obj, elapsed = gurobi_solve(**params)
        #     gap_cvx.append(abs(obj - cvx_obj))
        #     gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        #     times_cvx.append(elapsed)
        # if gap_cvx[-1] != 0:
        #     params['max_iters'] = 100
        #     obj, elapsed = gurobi_solve(**params)
        #     gap_cvx.append(abs(obj - cvx_obj))
        #     gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        #     times_cvx.append(elapsed)
        # if gap_cvx[-1] != 0:
        #     params['max_iters'] = 200
        #     obj, elapsed = gurobi_solve(**params)
        #     gap_cvx.append(abs(obj - cvx_obj))
        #     gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        #     times_cvx.append(elapsed)
        # if gap_cvx[-1] != 0:
        #     params['max_iters'] = 250
        #     obj, elapsed = gurobi_solve(**params)
        #     gap_cvx.append(abs(obj - cvx_obj))
        #     gap_perc_cvx.append(abs(1 - obj / cvx_obj))
        #     times_cvx.append(elapsed)
        # pdb.set_trace()
        # times_cvx, gap_cvx, gap_perc_cvx = process_gurobi_log(cvx_obj, time_elapsed)
        # del_gurobi_log()

        # save(times_cvx, 'times_cvx' + str(i), n)
        # save(gap_perc_cvx, 'gap_alg3' + str(i), n)
        # save(gap_cvx, 'gap_cvx' + str(i), n)


def bisect_solve_exp(**kwargs):
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

    lam_bar = kwargs['lam_bar']
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
    # algo_obj = (lam_bar * sigmacost + G.mu_cost() *
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
        algo_obj= (lam_bar * sigmacost + G.mu_cost())

        # algo_obj = (mid * varcost + G.mu_cost() *
        #             kwargs['mu_scalar']) * kwargs['obj_scalar']
        times.append(time.time() - start - disc_tot +
                     nmcc_time_tot - scc_time_tot + nmcc_time_start)

        objs.append(algo_obj)
        lams.append(mid)
# 
        # print('algo_obj: {}, f: {}, lam: {}'.format(algo_obj, f, mid))

        disc_tot += time.time() - disc_me
        f = mid - float(lam_bar) / float((2 * sigmacost))

        if abs(f) <= kwargs['stop_tol']:
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

    return algo_obj, iters, np.array(subproblem_times), times, objs, lams, 5, G.mu_cost(), np.sqrt(G.var_cost())
# def read_log(logfilename):
#     f = open(logfilename + '_processed.txt', 'w')
#     ff = open(logfilename + '.txt', 'r+')
#     w_to_processed = False

#     for line in ff.readlines():
#         if line.find('ITE PFEAS') >= 0:
#             w_to_processed = True
#         if line.find('Optimizer terminated') >= 0:
#             w_to_processed = False
#         if w_to_processed:
#             f.write(line)
#     f.close()

#     log = pandas.read_table(logfilename + '_processed.txt',
#                             delim_whitespace=True, sep='\t', lineterminator='\n')
#     return log
# def update_weight(R, G, lam, delta=0):
#     for u, v, e in R.edges(data=True):
#         try:
#             w = (e['flow'] + delta) + 2 * lam * \
#                  G[u][v]['var'] * (e['flow'] + delta)
#         except:
#             w = (e['flow'] + delta) + 2 * lam * \
#                  G[v][u]['var'] * (e['flow'] + delta)

#         e['weight'] = w
# # __all__ = ['capacity_scaling']
# # from itertools import chain
# # from math import log
# # import networkx as nx
# # from networkx.utils import *
# def _detect_unboundedness(R):
#     """Detect infinite-capacity negative cycles.
#     """
#     s = generate_unique_node()
#     G = nx.DiGraph()
#     G.add_nodes_from(R)

#     # Value simulating infinity.
#     inf = R.graph['inf']
#     # True infinity.
#     f_inf = float('inf')
#     for u in R:
#         for v, e in R[u].items():
#             # Compute the minimum weight of infinite-capacity (u, v) edges.
#             w = f_inf
#             for k, e in e.items():
#                 if e['capacity'] == inf:
#                     w = min(w, e['weight'])
#             if w != f_inf:
#                 G.add_edge(u, v, weight=w)

#     if nx.negative_edge_cycle(G):
#         raise nx.NetworkXUnbounded(
#             'Negative cost cycle of infinite capacity found. '
#             'Min cost flow may be unbounded below.')
# def _build_residual_network(G, demand, capacity, weight):
#     """Build a residual network and initialize a zero flow.
#     """
#     if sum(G.node[u].get(demand, 0) for u in G) != 0:
#         raise nx.NetworkXUnfeasible("Sum of the demands should be 0.")

#     R = nx.MultiDiGraph()
#     R.add_nodes_from((u, {'excess': G.node[u].get(demand, 0),
#                           'potential': 0}) for u in G)

#     inf = float('inf')
#     # Detect selfloops with infinite capacities and negative weights.
#     for u, v, e in G.selfloop_edges(data=True):
#         if e.get(weight, 0) < 0 and e.get(capacity, inf) == inf:
#             raise nx.NetworkXUnbounded(
#                 'Negative cost cycle of infinite capacity found. '
#                 'Min cost flow may be unbounded below.')

#     # Extract edges with positive capacities. Self loops excluded.
#     if G.is_multigraph():
#         edge_list = [(u, v, k, e)
#                      for u, v, k, e in G.edges(data=True, keys=True)
#                      if u != v and e.get(capacity, inf) > 0]
#     else:
#         edge_list = [(u, v, 0, e) for u, v, e in G.edges(data=True)
#                      if u != v and e.get(capacity, inf) > 0]
#     # Simulate infinity with the larger of the sum of absolute node imbalances
#     # the sum of finite edge capacities or any positive value if both sums are
#     # zero. This allows the infinite-capacity edges to be distinguished for
#     # unboundedness detection and directly participate in residual capacity
#     # calculation.
#     inf = max(sum(abs(R.node[u]['excess']) for u in R),
#               2 * sum(e[capacity] for u, v, k, e in edge_list
#                       if capacity in e and e[capacity] != inf)) or 1
#     for u, v, k, e in edge_list:
#         r = min(e.get(capacity, inf), inf)
#         w = e.get(weight, 0)
#         # Add both (u, v) and (v, u) into the residual network marked with the
#         # original key. (key[1] == True) indicates the (u, v) is in the
#         # original network.
#         R.add_edge(u, v, key=(k, True), capacity=r, weight=w, flow=0)
#         R.add_edge(v, u, key=(k, False), capacity=0, weight=-w, flow=0)

#     # Record the value simulating infinity.
#     R.graph['inf'] = inf

#     _detect_unboundedness(R)

#     return R
# def _build_flow_dict(G, R, capacity, weight):
#     """Build a flow dictionary from a residual network.
#     """
#     inf = float('inf')
#     flow_dict = {}
#     if G.is_multigraph():
#         for u in G:
#             flow_dict[u] = {}
#             for v, es in G[u].items():
#                 flow_dict[u][v] = dict(
#                     # Always saturate negative selfloops.
#                     (k, (0 if (u != v or e.get(capacity, inf) <= 0 or
#                                e.get(weight, 0) >= 0) else e[capacity]))
#                     for k, e in es.items())
#             for v, es in R[u].items():
#                 if v in flow_dict[u]:
#                     flow_dict[u][v].update((k[0], e['flow'])
#                                            for k, e in es.items()
#                                            if e['flow'] > 0)
#     else:
#         for u in G:
#             flow_dict[u] = dict(
#                 # Always saturate negative selfloops.
#                 (v, (0 if (u != v or e.get(capacity, inf) <= 0 or
#                            e.get(weight, 0) >= 0) else e[capacity]))
#                 for v, e in G[u].items())
#             flow_dict[u].update((v, e['flow']) for v, es in R[u].items()
#                                 for e in es.values() if e['flow'] > 0)
#     return flow_dict
# def capacity_scaling(G, demand='demand', capacity='capacity', weight='weight',heap=BinaryHeap, lam=0.7):

#     R = _build_residual_network(G, demand, capacity, weight)

#     inf = float('inf')
#     # Account cost of negative selfloops.
#     flow_cost = sum(
#         0 if e.get(capacity, inf) <= 0 or e.get(weight, 0) >= 0
#         else e[capacity] * e[weight]
#         for u, v, e in G.selfloop_edges(data=True))

#     # Determine the maxmimum edge capacity.
#     wmax = max(chain([-inf],
#                      (e['capacity'] for u, v, e in R.edges(data=True))))
#     if wmax == -inf:
#         # Residual network has no edges.
#         return flow_cost, _build_flow_dict(G, R, capacity, weight)

#     R_node = R.node
#     R_succ = R.succ

#     delta = 2 ** int(log(wmax, 2))
#     update_weight(R, G, lam)
#     while delta >= 0.0025:
#         # Determine the Δ-active nodes.
#         S = set()
#         T = set()
#         S_add = S.add
#         S_remove = S.remove
#         T_add = T.add
#         T_remove = T.remove
#         for u in R:
#             excess = R_node[u]['excess']
#             if excess >= delta:
#                 S_add(u)
#             elif excess <= -delta:
#                 T_add(u)
#         # Repeatedly augment flow from S to T along shortest paths until
#         # Δ-feasibility is achieved.
#         print(delta)
#         while S and T:
#             s = next(iter(S))
#             t = None
#             # Search for a shortest path in terms of reduce costs from s to
#             # any t in T in the Δ-residual network.
#             d = {}
#             pred = {s: None}
#             h = heap()
#             h_insert = h.insert
#             h_get = h.get
#             h_insert(s, 0)
#             while h:
#                 u, d_u = h.pop()
#                 d[u] = d_u
#                 if u in T:
#                     # Path found.
#                     t = u
#                     break
#                 p_u = R_node[u]['potential']
#                 for v, es in R_succ[u].items():
#                     if v in d:
#                         continue
#                     wmin = inf
#                     # Find the minimum-weighted (u, v) Δ-residual edge.
#                     for k, e in es.items():
#                         if e['capacity'] - e['flow'] >= delta:
#                             w = e['weight']
#                             if w < wmin:
#                                 wmin = w
#                                 kmin = k
#                                 emin = e
#                     if wmin == inf:
#                         continue
#                     # Update the distance label of v.
#                     d_v = d_u + wmin - p_u + R_node[v]['potential']
#                     if h_insert(v, d_v):
#                         pred[v] = (u, kmin, emin)
#             # pdb.set_trace()
#             if t is not None:
#                 # Augment Δ units of flow from s to t.
#                 while u != s:
#                     v = u
#                     u, k, e = pred[v]
#                     e['flow'] += delta
#                     R_succ[v][u][(k[0], not k[1])]['flow'] -= delta
#                 # Account node excess and deficit.
#                 R_node[s]['excess'] -= delta
#                 R_node[t]['excess'] += delta
#                 if R_node[s]['excess'] < delta:
#                     S_remove(s)
#                 if R_node[t]['excess'] > -delta:
#                     T_remove(t)
#                 # Update node potentials.
#                 d_t = d[t]
#                 for u, d_u in d.items():
#                     R_node[u]['potential'] -= d_u - d_t
#             else:
#                 # Path not found.
#                 S_remove(s)
#             update_weight(R, G, lam)

#         # Saturate Δ-residual edges with negative reduced costs to achieve
#         # Δ-optimality.
#         for u in R:
#             p_u = R_node[u]['potential']
#             for v, es in R_succ[u].items():
#                 for k, e in es.items():
#                     flow = e['capacity'] - e['flow']
#                     if e['weight'] - p_u + R_node[v]['potential'] < 0:
#                         flow = e['capacity'] - e['flow']
#                         if flow >= delta:
#                             e['flow'] += flow
#                             R_succ[v][u][(k[0], not k[1])]['flow'] -= flow
#                             R_node[u]['excess'] -= flow
#                             R_node[v]['excess'] += flow
#         update_weight(R, G, lam)

#         delta = delta / 2.0

#     # pdb.set_trace()
#     if any(R.node[u]['excess'] > 2 * delta for u in R):
#         raise nx.NetworkXUnfeasible('No flow satisfying all demands.')

#     # Calculate the flow cost.
#     for u in R:
#         for v, es in R_succ[u].items():
#             for e in es.values():
#                 flow = e['flow']
#                 if flow > 0:
#                     flow_cost += flow * e['weight']

#     return flow_cost, _build_flow_dict(G, R, capacity, weight)
# def push_flow(G, P, delta, g):
#     for i in range(len(P)-1):
#         u = P[i]
#         v = P[i+1]

#         if G.nxg.has_edge(u, v):
#             flow =  G.nxg[u][v]['flow']
#             G.nxg[u][v]['flow'] = flow + delta
#         else:
#             flow =  G.nxg[v][u]['flow']
#             G.nxg[v][u]['flow'] = flow - delta
#         g[u-1] -= delta
#         g[v-1] += delta

#     return g
# def ASSP(G, lam):
#     nx.set_node_attributes(G.nxg, 0, 'potential')
#     nx.set_edge_attributes(G.nxg, 0, 'flow')

#     theta = 0.5
#     rho = 1.0 / 3.0

#     m = len(G.nxg.edges())

#     mu = np.zeros(m)
#     var = np.zeros(m)
#     sigma = np.zeros(m)
#     cap = np.zeros(m)
#     d = np.zeros(len(G.nxg.nodes()))
#     rows = []
#     cols = []
#     values = []
#     i = 0
#     arc_dict = {}
#     for u, v, e in G.nxg.edges(data=True):
#         pdb.set_trace()
#         mu[i] = e.get('mu', 0)
#         sigma[i] = np.sqrt(e.get('var', 0))
#         var[i] = e.get('var', 0)
#         cap[i] = e.get('capacity', 0)
#         arc_dict[i] = (u, v)
#         rows.append(u - 1)
#         cols.append(i)
#         values.append(1)
#         rows.append(v - 1)
#         cols.append(i)
#         values.append(-1)
#         i += 1

#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         d[i] = dat['demand']
#         i += 1


#     max_d = max(d)
#     max_d_vec = np.ones(len(cap)) * max_d
#     x_carrot = np.minimum(cap, max_d_vec)

#     eps_max = np.inf
#     eps = min(eps_max, max(1.0 / (len(G.nxg.nodes()) + 1),
#               rho * np.max(mu + var * x_carrot)))

#     # x = np.zeros(len(cap))
#     # for i,j,e in G.nxg.edges(data=True):
#     #     G.nxg[i][j]['flow'] = max(0, min(G.nxg[i][j]['capacity'], - (G.nxg[i][j]['mu'] - (G.nxg._node[i]['potential'] - G.nxg._node[j]['potential'] )/(2*lam*G.nxg[u][v]['var']))))

#     g = np.zeros(len(G.nxg.nodes()))
#     for i in G.nxg.nodes():
#         outgoing = 0
#         for j, es in G.nxg._succ[i].items():
#             outgoing += G.nxg[i][j]['flow']
#         incoming = 0
#         for j, es in G.nxg._pred[i].items():
#             incoming += G.nxg[j][i]['flow']

#         g[i-1] = outgoing - incoming + G.nxg._node[i]['demand']

#     surp = g[g > 0]
#     pos_surp = [ list(g).index(i) for i in list(surp)]
#     s = pos_surp[0] + 1
#     P = [s]
#     while len(g[g > 0]) > 0:

#         print('-------')
#         print(eps)
#         print(P)
#         print(g)
#         '''step 1'''
#         i = P[-1]
#         step2time = True

#         for j, es in G.nxg._succ[i].items():
#             if G.nxg._node[i]['potential'] - G.nxg._node[j]['potential'] > G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow'] + theta * eps:
#                 step2time = False
#                 status = 'active'
#                 addme = j
#                 break

#         if step2time:
#             for j, es in G.nxg._pred[i].items():
#                 if G.nxg._node[j]['potential'] - G.nxg._node[i]['potential'] < G.nxg[j][i]['mu'] + 2*lam*G.nxg[j][i]['var']*G.nxg[j][i]['flow'] - theta * eps:
#                     step2time = False
#                     status = 'inactive'
#                     addme = j
#                     break

#         if not step2time:
#             P.append(addme)
#             if len(P) == len(list(set(P))):
#                 '''step 3a: push flow'''
#                 if g[addme-1] >= 0:
#                     continue
#                 else:

#                     fmarginP = np.inf
#                     for k in range(len(P)):
#                         if k == len(P)-1:
#                             break
#                         u = P[k]
#                         v = P[k+1]

#                         try:
#                             var = 2*lam*G.nxg[u][v]['var']
#                             deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] - (G.nxg[u][v]['mu'] + var*G.nxg[u][v]['flow']))/(var)
#                             # deltaij = min(deltaij, G.nxg[u][v]['capacity'] - G.nxg[u][v]['flow'])
                            
#                         except:
#                             var = 2*lam*G.nxg[v][u]['var']
#                             deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] + (G.nxg[v][u]['mu'] + var*G.nxg[v][u]['flow']))/(var)
#                             # deltaij = min(deltaij, G.nxg[v][u]['flow'])
                      
#                         if deltaij <= fmarginP:
#                             fmarginP = deltaij

#                     delta = min(fmarginP, g[s-1], -g[addme-1])
#                     g = push_flow(G, P, delta, g)
#                     surp = g[g > 0]
#                     pos_surp = [ list(g).index(o) for o in list(surp)]
#                     s = pos_surp[0] + 1
#                     P = [s]

#             else:
#                 '''step 3b: augment the flow on the cycle'''
#                 cs = set([x for x in P if P.count(x) > 1])
#                 C = P[P.index(cs.pop()):]
#                 fmarginC = np.inf

#                 for k in range(len(C)):
#                     if k == len(C)-1:
#                         break
#                     u = C[k]
#                     v = C[k+1]

#                     try:
#                         var = 2*lam*G.nxg[u][v]['var']
#                         deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] - (G.nxg[u][v]['mu'] + var*G.nxg[u][v]['flow']))/(var)
#                         deltaij = min(deltaij, G.nxg[u][v]['capacity'] - G.nxg[u][v]['flow'])
                        
#                     except:
#                         var = 2*lam*G.nxg[v][u]['var']
#                         deltaij = (G.nxg._node[u]['potential'] - G.nxg._node[v]['potential'] + (G.nxg[v][u]['mu'] + var*G.nxg[v][u]['flow']))/(var)
#                         deltaij = min(deltaij, G.nxg[v][u]['flow'])

#                     if deltaij <= fmarginC:
#                         fmarginC = deltaij

#                 g = push_flow(G, P, fmarginC, g)
#                 surp = g[g > 0]
#                 pos_surp = [ list(g).index(o) for o in list(surp)]
#                 s = pos_surp[0] + 1
#                 P = [s]

#         else:
#             '''step 2: price rise on i'''
#             minpi = np.inf
#             for j, es in G.nxg._succ[i].items():
#                 curpi = G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow'] + eps + G.nxg._node[j]['potential']
#                 if curpi <= minpi:
#                     minpi = curpi
#             for j, es in G.nxg._pred[i].items():
#                 curpi = -G.nxg[j][i]['mu'] - 2*lam*G.nxg[j][i]['var']*G.nxg[j][i]['flow'] + eps + G.nxg._node[j]['potential']
#                 if curpi <= minpi:
#                     minpi = curpi

#             # delta_1 = np.inf
#             # for j, es in G.nxg._succ[i].items():
#             #     curpi = G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow'] + eps - (G.nxg._node[i]['potential'] - G.nxg._node[j]['potential'])
#             #     if curpi <= delta_1:
#             #         delta_1 = curpi
#             # delta_2 = np.inf
#             # for j, es in G.nxg._succ[i].items():
#             #     curpi = -(G.nxg[i][j]['mu'] + 2*lam*G.nxg[i][j]['var']*G.nxg[i][j]['flow']) + eps + (G.nxg._node[j]['potential'] - G.nxg._node[i]['potential'])
#             #     if curpi <= delta_2:
#             #         delta_2 = curpi

#             # delta_3 = np.inf
#             # minpi = min(delta_1, delta_2, delta_3)

#             G.nxg._node[i]['potential'] = minpi
#             if i != s:
#                 P.remove(i)

#     pdb.set_trace()
#     print('wait here')
# def dcinv(arc, t, lam=0.7, rev=False):
#     if not rev:
#         return (t - arc['mu'])/(2.0*lam*arc['var'])
# def tseng_ASSG(G):
#     nx.set_node_attributes(G.nxg, 0, 'potential')
#     nx.set_edge_attributes(G.nxg, 0, 'flow')

#     theta = 0.5
#     rho = 1.0 / 3.0

#     startn = {}
#     endn = {}
#     fin = {}
#     fout = {}
#     nxtin = {}
#     nxtou = {}
#     finalin = {}
#     finalou = {}

#     m = len(G.nxg.edges())
#     n = len(G.nxg.nodes())

#     for node in range(1,n+1):
#         fin[node] = 0
#         fout[node] = 0
#         finalin[node] = 0
#         finalou[node] = 0

#     for s,e,arc in G.nxg.edges(data=True):
#         if fout[s] != 0:
#             nxtou[finalou[s]] = arc
#         else:
#             fout[s] = (arc, e)
#         if fin[e] != 0:
#             nxtin[finalin[e]] = arc
#         else:
#             fin[e]=arc
#         finalou[s]=arc
#         finalin[e]=arc
#         nxtin[arc]=0
#         nxtou[arc]=0

#     for u,v,e in G.nxg.edges(data=True):
#         qfpushf[e] = (0, u, v, e)
#         qfpushb[e] = (0, u, v, e)
#         p[e] = 0

#     mu = np.zeros(m)
#     var = np.zeros(m)
#     sigma = np.zeros(m)
#     cap = np.zeros(m)
#     d = np.zeros(len(G.nxg.nodes()))
#     rows = []
#     cols = []
#     values = []
#     i = 0
#     arc_dict = {}
#     for u, v, e in G.nxg.edges(data=True):
#         mu[i] = e.get('mu', 0)
#         sigma[i] = np.sqrt(e.get('var', 0))
#         var[i] = e.get('var', 0)
#         cap[i] = e.get('capacity', 0)
#         arc_dict[i] = (u, v)
#         rows.append(u - 1)
#         cols.append(i)
#         values.append(1)
#         rows.append(v - 1)
#         cols.append(i)
#         values.append(-1)
#         i += 1

#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         d[i] = dat['demand']
#         i += 1

#     max_d = max(d)
#     max_d_vec = np.ones(len(cap)) * max_d
#     x_carrot = np.minimum(cap, max_d_vec)

#     eps_max = np.inf
#     eps = min(eps_max, max(1.0 / (len(G.nxg.nodes()) + 1),
#               rho * np.max(mu + var * x_carrot)))

#     DFCT1 = {}
#     PRDCSR = {}
#     g = np.zeros(len(G.nxg.nodes()))
#     for i in G.nxg.nodes():
#         outgoing = 0
#         for j, es in G.nxg._succ[i].items():
#             outgoing += G.nxg[i][j]['flow']
#         incoming = 0
#         for j, es in G.nxg._pred[i].items():
#             incoming += G.nxg[j][i]['flow']

#         g[i-1] = outgoing - incoming + G.nxg._node[i]['demand']
#         DFCT1[i] = -g[i-1]
#         PRDCSR[i] = 0
#         SURPLUS[u] = g[i-1]

#     LARGE=2000000000
#     REF_THR = e-10
#     MAXPRICE = LARGE
#     TOTALITER = 0
#     niter = 0


#     #reduce arc capacities
#     for node in G.nxg.nodes():
#         CAPOUT = 0
#         arc = fout[node]
#         if arc > 0:
#             CAPOUT = min(LARGE - CAPOUT - arc['capacity'])
#             arc = nxtou[arc]
#         CAPOUT = min(LARGE - CAPOUT - SURPLUS[i])
#         if CAPOUT < 0:
#             print('infeasible')
#         CAPIN = 0
#         arc = fin[node]
#         if arc > 0:
#             if arc['capacity'] > CAPOUT:
#                 arc['capacity'] = CAPOUT
#             CAPIN = min(LARGE, CAPIN+arc['capacity'])
#             arc = ntxtin[arc]
#         CAPIN = min(LARGE, CAPIN + SURPLUS[node])
#         if CAPIN < 0:
#             print('infeasible')
#         arc = fout[node]
#         if arc > 0:
#             if arc['capacity'] > CAPIN:
#                 arc['capacity'] = CAPIN
#             arc = nxtou[arc]

#     #set arc flows to satisfy e-cs
#     eps2 = eps/2.0
#     eps3 = eps2/2.0

#     # for u,v,arc in G.nxg.edges(data=True):
#     #     start = u
#     #     end = v
#     #     PSTART = P[start]
#     #     PEND = P[end]
#     ###########

#     for i in G.nxg.nodes():
#         PRDCSR[i] = 0

#     NXTQUEUE = {}
#     for i in range(1, n):
#         NXTQUEUE[i] = i+1

#     NXTQUEUE[n]=1
#     NODE = 1
#     PREVNODE = NODE
#     LASTQUEUE = n

#     TOTALITER += TOTALITER
#     PRDCSR[NODE] = m + 1
#     i = NODE

#     qarcf = qfpushf[i]
#     qarcb = qfpushb[i]
#     price = p[i]

#     ###Forward arcs
#     if qarcf > 0:
#         j = qarcf[2]
#         if price - p[j] - (qarcf[3]['mu'] + 2*lam*qarcf[3]['var']*qarcf[3]['flow']) > eps2:
#             resid = min(qarcf['capacity'], dcinv(quarcf, price-p[j], lam=lam)) - qarcf[3]['flow']
#             if resid > REF_THR:
#                 qfpushf[i] = qarcf
#                 if -SURPLUS[j] > REF_THR:
#                     path_push1(j, qarcf[0])
#                     # go to 620
#                 else:
#                     if PRDCSR[j] == 0:
#                         PRDCSR[j] = qarcf
#                         i = j
#                         #go to 110
#                     else:
#                         cycle_push1(i,j,qarcf[0])
#                         #go to 620
#             else:
#                 qarcf = qnxtpushf[qarcf]
#         else:
#             qarcf = qnxtpushf[qarcf]
#         #go to 122
#     qfpushf[i] = qarcf

#     ###backward arcs
#     if qarcb > 0:
#         j = qarcb[2]
#         if price - p[j] - (qarcb[3]['mu'] + 2*lam*qarcb[3]['var']*qarcb[3]['flow']) > eps2:
#             resid = qarcb[3]['flow']  - max(0, dcinv(qarcb, p[j] - price, lam=lam)) 
#             if resid > REF_THR:
#                 qfpushb[i] = qarcb
#                 if -SURPLUS[j] > REF_THR:
#                     path_push1(j, -qarcb[0])
#                     # go to 620
#                 else:
#                     if PRDCSR[j] == 0:
#                         PRDCSR[j] = -qarcb
#                         i = j
#                         #go to 110
#                     else:
#                         cycle_push1(i,j,-qarcb[0])
#                         #go to 620
#             else:
#                 qarcb = qnxtpushf[qarcb]
#         else:
#             qarcb = qnxtpushf[qarcb]
#         #go to 123
#     qfpushb[i] = qarcb

#     if SURPLUS[i] > 0:
#         print(i, qarcf, qarcb, price)

#     if qarcb[0] + qarcf[0] == 0:
#         refprice = price
#         price = LARGE
#         niter = niter + 1
#         arc = fout[i][0]
#         if arc > 0:
#             up = arc['capacity']
#             if (arc['flow'] - up) < -REF_THR:
#                 xp = p[fout[i][1]] + (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -eps:
#                     price = xp + eps

#                     qarcf = arc
#                     qnxtpushf[arc] = 0
                    
#                 else:
#                     if (xp - price) < -eps2:
#                         if qarcf == 0:
#                             qarcf = arc
#                             qnxtpush[arc] = 0
#                         else:
#                             qnxtpushf[arc] = qarcf
#                             qarcf = arc
#             arc = nxtou[arc]
#             #goto 111

#         arc = fin[i][0]
#         if arc > 0:
#             if arc['flow'] > REF_THR:
#                 xp = p[fin[i][1]] - (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -eps:
#                     price = xp + eps
#                     qarcf = 0
#                     qarcb = arc
#                     qnxtpushb[arc] = 0
#                 else:
#                     if (xp-price) < -eps2:
#                         if qarcb == 0:
#                             qarcb = arc
#                             qnxtpushb[qarcb] = 0
#                         else:
#                             qnxtpushb[arc] = qarcb
#                             qarcb = arc

#             arc = nxtin[arc]
#             # goto 112


#         if price > LARGE:
#             qarcf = 0
#             qarcb = 0
#             arc = fout[i][0]
#             if arc > 0:
#                 xp = p[fout[i][1]] + (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -REF_THR:
#                     price = xp
#                 else:
#                     if xp > LARGE:
#                         print('price exceeds the range')
#                 arc = nxtou[arc]
#                 #goto 211
#             arc = fin[i]

#             if arc > 0:
#                 xp = p[fout[i][1]] - (arc['mu'] + lam*arc['var']*arc['flow']*arc['flow'])
#                 if (xp-price) < -REF_THR:
#                     price = xp
#                 else:
#                     if xp > LARGE:
#                         print('price exceeded the range')
#                 arc = nxtin[arc]
#                 #go to 212

#             if price > int(LARGE/10)
#                 price = int(LARGE/10)

#         p[i] = price
#         qfpushf[i] = qarcf
#         qfpushb[i] = qarcb


#         arc = PRDCSR[i]
#         if arc != m+1:
#             PRDCSR[i] = 0
#             if arc > 0:
#                 i = 0#startnode of arc
#             else:
#                 i = 0#endnode of -arc

#         # goto 110


#     if -SURPLUS[NODE] > THRESHSURPLUS:
#         TOTALITER +=1
#         PRDCSR[NODE] = m+1
#         i = NODE
# def subgrad(G, lam):
#     n = len(G.nxg.nodes())
#     m = len(G.nxg.edges())


#     var = np.zeros(m)
#     cap = np.zeros(m)
#     mu = np.zeros(m)

#     arc_dict = {}
#     A = np.zeros((n,m))
#     i = 0
#     for u, v, e in G.nxg.edges(data=True):
#         var[i] = e.get('var', 0)
#         cap[i] = e.get('capacity', 0)
#         mu[i] = e.get('mu', 0)
#         A[u-1, i] = 1
#         A[v-1, i] = -1
#         arc_dict[i] = (u,v)
#         i += 1

#     b = np.zeros(n)
#     i = 0
#     for node, dat in G.nxg.nodes(data=True):
#         b[i] = dat['demand']
#         i += 1

#     sumvar = 2*lam*var*np.eye(m)
#     sumvarinv = np.linalg.inv(sumvar)
#     H = np.dot(np.dot(A, sumvarinv),A.T)
#     Hinv = np.linalg.inv(H)

#     print('starting')
#     st = time.time()
#     p = np.zeros(n)
#     alpha = 1
#     g_p = np.zeros(n)
#     g_p[0] = 100
#     g_p[-1] = -100
#     close = False
#     iters = 1
#     while not close:
#         iters += 1

#         for i in G.nxg.nodes():
#             ij_term = 0
#             ji_term = 0
#             for j in np.where(A[i-1,:]==-1)[0]:
#                 idx = j
#                 u = arc_dict[j][0]
#             # for j, es in G.nxg._pred[i].items():
#                 # x_ji = max(0, min(G.nxg[j][i]['capacity'],(p[j-1] - p[i-1] - G.nxg[j][i]['mu'])/(2*lam*G.nxg[j][i]['var'])))
#                 x_ji = max(0, min(cap[idx],(p[u-1] - p[i-1] - mu[idx])/(2*lam*var[idx])))
#                 ji_term += x_ji
            
#             # for j, es in G.nxg._succ[i].items():
#             for j in np.where(A[i-1,:]==1)[0]:
#                 idx = j
#                 v = arc_dict[j][1]
#                 x_ij = max(0, min(cap[idx],(p[i-1] - p[v-1] - mu[idx])/(2*lam*var[idx])))
#                 ij_term += x_ij
            
#             g_p[i-1] = ji_term - ij_term + b[i-1]
#         p = p -alpha*np.dot(-Hinv,g_p)
#         if np.linalg.norm(g_p) <= 0.5:
#             close = True

#     print('time: ', time.time()-st)
#     for i,j,e in G.nxg.edges(data=True):
#             G.nxg[i][j]['flow'] = max(0, min(G.nxg[i][j]['capacity'],(p[i-1] - p[j-1] - G.nxg[i][j]['mu'])/(2*lam*G.nxg[i][j]['var'])))
#     tot_cost = G.mu_cost() + lam*G.var_cost()
#     print('obj: ', tot_cost)
#     pdb.set_trace()

