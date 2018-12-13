def run_GP(module, function, int_vars = [], values = [],direction = [], verbosity = 0):
    exec("from %s import %s as Model"%(module, function))
    M = Model()
    for i,var  in enumerate(int_vars):
        if values[i] != None:
            s = 'M.'+var
            M.substitutions[eval(s)]=values[i]
    try:
        sol = (M.solve(verbosity=verbosity,solver="mosek"))
        #print sol.summary()
    except:
        print "Likely infeasible model at %s\n %s"%(int_vars, values)
        return (None, None,None)

    var_l = []
    vals = []
    for var in sol["freevariables"]:
        var_l.append(var.__str__())
        vals.append(sol["freevariables"][var])

    return (var_l, vals, sol["cost"])

def con_GP(module, function, int_vars = [], values = [],direction = [], verbosity = 0):
    #Like run_GP, adds inequality constraints instead of substitutions.
    #Direction of the constraint is controlled by the direction value associated with it
    #">=" is _geq_ constraint
    #"<=" is _leq_ constraint
    #"==" is an equality constraint
    exec("from %s import %s as Model"%(module, function))
    M = Model()
    for i,var  in enumerate(int_vars):
        if values[i] != None:
            if len(direction) > 0:
                s = 'M.'+var+" " + direction[i]+ " "+str(values[i])
                M.append([eval(s)])
            else:
                s = 'M.'+var
                M.substitutions[eval(s)]=values[i]
    M.cost =  M.aircraft.Wtotal
    try:
        #print M
        sol = (M.solve(verbosity=verbosity,solver="mosek"))
        #print sol.summary()
    except:
        print "Likely infeasible model at %s\n %s"%(int_vars, values)
        return (None, None,None)

    var_l = []
    vals = []
    for var in sol["freevariables"]:
        var_l.append(var.__str__())
        vals.append(sol["freevariables"][var])

    return (var_l, vals, sol["cost"])

if __name__ == "__main__":
    con_GP("evtol", "Mission", ["INT_N_lift_motors", "INT_N_thrust_motors"], [2, 4], [">=", "=="])
