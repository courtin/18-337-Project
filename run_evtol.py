from evtol import Aircraft, Mission

def run_opt(int_vars=[], values=[]):
        Vehicle = Aircraft()
        M = Mission(Vehicle)
        M.substitutions[M.R_total] = 60
        for i,var  in enumerate(int_vars):
            if values[i] != None:
                s = 'M.'+var
                M.substitutions[eval(s)]=values[i]
        M.cost =  M.aircraft.Wtotal
        try:
            sol = (M.solve(verbosity=0,solver="mosek"))
            #print sol.summary()
        except RuntimeWarning:
            print "Likely infeasible model at %s\n %s"%(int_vars, values)
            return (None, None,None)

        var_l = []
        vals = []
        for var in sol["freevariables"]:
            var_l.append(var.__str__())
            vals.append(sol["freevariables"][var])

        return (var_l, vals, sol["cost"])
def none_test():
    return None
if __name__ == "__main__":
    (va, vv, cost) = run_opt()
    print cost
