from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled, Variable, SignomialEquality, units

class Test1(Model):
    """ Test equation

    Variables
    ---------
    INT_x1                  [-]     x1
    INT_x2                  [-]     x2
    INT_x3                  [-]     x3
    Z                   [-]     Z

    """
    def setup(self):
        exec parse_variables(Test1.__doc__)
        ZERO = 1e-30
        x1 = INT_x1
        x2 = INT_x2
        x3 = INT_x3
        constraints = [Z >= x1+x2,
            6*x2 >= 10+ 2*x3,
            4*x1 >= 3*x3,
            2*x1 >= x2,
            x1 >= ZERO,
            x2 >= ZERO,
            x3 >= ZERO

        ]

        self.cost = Z

        return constraints


def run_Test1(int_vars=[], values=[], verbosity = 0):
    #Comment
    M = Test1()
    for i,var  in enumerate(int_vars):
        if values[i] != None:
            s = 'M.'+var
            M.substitutions[eval(s)]=values[i]
    
    try:
        sol = (M.solve(verbosity=verbosity,solver="mosek"))
        #print sol.table()
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
    #(M1, M2, M3) = run_Test1(verbosity=1)
    #(M1, M2, M3) = run_Test1(["INT_x3"], [1e-30],verbosity=1)
    (M1, M2, M3) = run_Test1(["INT_x3","INT_x2"], [2, 1e-30],verbosity=1)
    #(M1, M2, M3) = run_Test1(["INT_x2","INT_x1", "INT_x3"], [2,1, 1e-30],verbosity=1)


