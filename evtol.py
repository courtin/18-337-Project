" eVTOL Aircraft Model "
#pylint: disable=invalid-name, too-many-instance-attributes, too-many-locals
#pylint: disable=redefined-variable-type, too-many-statements, not-callable
from os.path import abspath, dirname
from os import sep
from numpy import hstack
import pandas as pd
from gpkit import Model, parse_variables, Vectorize, SignomialsEnabled, Variable, SignomialEquality, units
from gpkitmodels.GP.aircraft.wing.wing_core import WingCore
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar as BoxSparGP
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSecondStruct
from gpkitmodels.GP.aircraft.tail.horizontal_tail import HorizontalTail
from gpkitmodels.GP.aircraft.tail.vertical_tail import VerticalTail
from gpkitmodels.GP.materials import cfrpud, cfrpfabric
from gpkitmodels.GP.aircraft.fuselage.elliptical_fuselage import Fuselage
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp2
from gpkitmodels.GP.aircraft.motor.motor import Motor
from gpkitmodels import g
from my_fit_constraintset import FitCS as FCS
from numpy import pi
from gpkit.constraints.tight import Tight as TCS
#import  matplotlib.pyplot as plt


class AircraftPerf(Model):
    """ Aircaft Performance

    Variables
    ---------
    CD                  [-]     aircraft drag coefficient
    cda                 [-]     non-wing drag coefficient
    mfac        1.0     [-]     drag margin factor
    mpower      1.0     [-]     power margin
    T                   [lbf]   thrust
    L                   [lbf]   Total lift
    L_wing              [lbf]   Wing lift
    L_rotors            [lbf]   Lift from dedicated rotors
    L_tilt              [lbf]   Lift from tilt rotors
    T_tilt              [lbf]   Thrust from tilt rotors
    T_rotors            [lbf]   Thrust from dedicated rotors
    W                   [lbf]   Current Weight
    Poper               [kW]    Operating power
    P_rotors            [kW]    Power going to rotors
    Pavn        500     [W]     Avioncs power
    Ppay        500     [W]     Payload power
    D                   [lbf]   Drag
    """
    def setup(self, static, state):
        exec parse_variables(AircraftPerf.__doc__)

        fd = dirname(abspath(__file__)) + sep + "dai1336a.csv"

        self.wing = static.wing.flight_model(static.wing, state, fitdata=fd)
        self.htail = static.htail.flight_model(static.htail, state)
        self.vtail = static.vtail.flight_model(static.vtail, state)
        self.thrust_prop = static.thrust_prop.flight_model(static.thrust_prop, state)
        self.lift_prop = static.lift_prop.flight_model(static.lift_prop, state)
        self.thrust_motor = static.thrust_motor.flight_model(static.thrust_motor, state)
        self.lift_motor = static.lift_motor.flight_model(static.lift_motor, state)
        self.boom = static.boom.flight_model(static.boom, state)
        self.flight_models = [self.wing, self.htail, self.vtail,self.boom, 
                                self.thrust_prop, self.lift_prop,
                                self.thrust_motor, self.lift_motor
                                ]

        e = self.e = self.wing.e
        cdht = self.cdht = self.htail.Cd
        cdvt = self.cdvt = self.vtail.Cd
        Sh = self.Sh = static.htail.planform.S
        Sv = self.Sv = static.vtail.planform.S
        Sw = self.Sw = static.wing.planform.S
        cdw = self.cdw = self.wing.Cd
        CL = self.CL = self.wing.CL


        Qprop_thrust = self.thrust_prop.Q
        RPMprop_thrust = self.thrust_prop.omega
        Qmotor_thrust = self.thrust_motor.Q
        RPMmotor_thrust = self.thrust_motor.omega

        Qprop_lift = self.lift_prop.Q
        RPMprop_lift = self.lift_prop.omega
        Qmotor_lift = self.lift_motor.Q
        RPMmotor_lift = self.lift_motor.omega
        

        self.wing.substitutions[e] = 0.85
        self.wing.substitutions[self.wing.CLstall] = 1.4
        dvars = [cdht*Sh/Sw, cdvt*Sv/Sw]
        self.fuse = static.fuselage.flight_model(static.fuselage,
                                                         state)
        self.flight_models.extend([self.fuse])
        cdfuse = self.fuse.Cd
        Sfuse = static.fuselage.S
        dvars.extend([cdfuse*Sfuse/Sw])
        self.fuse.substitutions[self.fuse.mfac] = 1.1
        
        cdboom = self.boom.Cd
        Sboom = static.boom.S
        dvars.extend([cdboom*Sboom/Sw*static.N_booms])
        Wzf = static.Wzf
        rho = state.rho
        V = state.V
        M = state.M
 
        #Sh = static.htail.planform.S
        #CLhmax = static.htail.planform.CLmax
        #Sv = static.vtail.planform.S
        #CLvmax = static.vtail.planform.CLmax




        constraints = [cda >= sum(dvars),
                       CD/mfac >= cda + cdw,
                       
                       Poper >= Pavn + Ppay + P_rotors, 

                       Qprop_thrust == Qmotor_thrust,
                       RPMprop_thrust == RPMmotor_thrust,
                       Qprop_lift == Qmotor_lift,
                       RPMprop_lift == RPMmotor_lift,
                       #cdht >= .005,
                       #cdvt >= .005
                      ]



        return self.flight_models, constraints

class Aircraft(Model):
    """ Aircraft Model

    Variables
    ---------
    Wpay        400     [lbf]    payload weight
    Wavn        40.     [lbf]    avionics weight
    Wtotal              [lbf]   aircraft weight
    Wzf                 [lbf]   zero fuel weight
    mfac        1.0      [-]     total weight margin
    Wh                  [lbf]     htail weight
    Wv                  [lbf]     vtail weight
    bmax        40      [ft]    span constraint
    minvttau    0.09    [-]     minimum vertical tail tau ratio
    minhttau    0.06    [-]     minimum horizontal tail tau ratio
    Wbatt               [kgf]   Battery weight
    Ebatt               [W*hr]  Battery Energy
    Vbatt               [cm^3]  Battery Volume
    Pavn        80      [W]     avionics power draw
    Ppay        60      [W]     payload power draw
    Poper               [W]     operating power
    T_flight            [min]     flight time
    Estar_batt  200     [W*hr/kgf] battery specific energy
    Vstar_batt  .047    [W*hr/cm^3] battery specific volume 
    N_lift_motors       [-]     Number of dedicated lifting motors
    N_thrust_motors     [-]     Number of dedicated thrust motors
    N_tilt_motors       [-]     Number of tilting motors
    N_booms             [-]     Number of mounting booms for motors
    W_f                 [lbf]   fuse weight
    f_h                 [-]     htail weight fraction
    f_v                 [-]     vtail weight fraction
    f_f         .2      [-]     fuselage weight fraction   
    R_fuse      5       [ft]     fuselage radius  
    L_fuse      25      [ft]     fuselage length
    R_boom      .15      [ft]    boom radius
    N_booms             [-]     number of booms

    SKIP VERIFICATION

    """
    fuseModel = None
    flight_model = AircraftPerf

    def setup(self):
        exec parse_variables(Aircraft.__doc__)

        self.fuselage = fuselage = Fuselage()
        self.boom = boom = Fuselage()

        HorizontalTail.sparModel = BoxSparGP
        HorizontalTail.fillModel = None
        HorizontalTail.skinModel = WingSecondStruct
        VerticalTail.sparModel = BoxSparGP
        VerticalTail.fillModel = None
        VerticalTail.skinModel = WingSecondStruct
        self.htail = HorizontalTail()
        self.vtail = VerticalTail()

        
        
        
        WingGP.sparModel = BoxSparGP
        WingGP.fillModel = None
        WingGP.skinModel = WingSecondStruct
        self.wing = WingGP(N=3)
        #self.wing.substitutions["Abar"] = .0492 #Firefly FFw3_g.048_root airfoil
        #self.wing.substitutions["tau"]  = .1095 #Constant thickness airfoil


        Sw      = self.Sw       = self.wing.planform.S
        cmac    = self.cmac     = self.wing.planform.cmac
        tau     = self.tau      = self.wing.planform.tau
        croot   = self.croot    = self.wing.planform.croot
        b       = self.b        = self.wing.planform.b
        Vh      = self.Vh       = self.htail.Vh
        lh      = self.lh       = self.htail.lh
        Sh      = self.Sh       = self.htail.planform.S
        Vv      = self.Vv       = self.vtail.Vv
        Sv      = self.Sv       = self.vtail.planform.S
        lv      = self.lv       = self.vtail.lv

        self.vtail.substitutions[Vv] = 0.04
        self.htail.substitutions[Vh] = 0.75

        vttau = self.vtail.planform.tau
        httau = self.htail.planform.tau
        
        self.wing.substitutions[self.wing.mfac] = 1.0
        
        #Volfuse = self.Volfuse = self.fuselage.Vol
        
        Wfuse = self.fuselage.W
        Wwing = self.wing.W
        Wboom = self.boom.W
        
        self.htail.substitutions[self.htail.mh] = 0.1
        self.htail.substitutions[self.htail.skin.rhoA] = 0.4
        self.vtail.substitutions[self.vtail.skin.rhoA] = 0.4
        self.htail.substitutions[self.htail.planform.AR] = 4
        self.vtail.substitutions[self.vtail.planform.AR] = 4
        self.htail.substitutions[self.htail.planform.CLmax] = 1.5
        self.vtail.substitutions[self.vtail.planform.CLmax] = 1.5

        self.fuselage.substitutions[self.fuselage.nply] = 5
        self.boom.substitutions[self.boom.nply] = 2

        self.lift_motor = Motor()
        self.thrust_motor = Motor()
        #self.tilt_motor = Motor()

        W_l_m = self.lift_motor.W
        W_thr_m = self.thrust_motor.W
        #W_tlt_m = self.tilt_motor.W

        Propeller.flight_model = ActuatorProp2
        self.lift_prop = Propeller()
        self.thrust_prop = Propeller()
        #self.tilt_prop = Propeller()

        W_l_p = self.lift_prop.W
        W_thr_p = self.thrust_prop.W

        self.propulsors = [self.lift_motor, self.thrust_motor, #self.tilt_motor,
                            self.lift_prop, self.thrust_prop,  #self.tilt_prop
                            ]

        self.components = [self.wing, self.fuselage, 
                            self.htail, self.vtail,
                            ]

        
 
        #self.loading = [self.htailg, vtailg]
        constraints = [Vh <= Sh*lh/Sw/cmac,
                       Vv <= Sv*lv/Sw/b,
                       fuselage.R == R_fuse,
                       fuselage.l == L_fuse,
                       lv ==  L_fuse,
                       lh == L_fuse,
                       vttau >= minvttau,
                       httau >= minhttau,
                       Wh == f_h*Wtotal,
                       Wv == f_v*Wtotal,
                       W_f == f_f*Wtotal,
                        TCS([Wtotal/mfac >= (Wpay + Wavn + Wbatt +
                        N_lift_motors*W_l_m + 
                        N_thrust_motors*W_thr_m +
                        N_lift_motors*W_l_p + 
                        N_thrust_motors*W_thr_p +
                        N_booms*Wboom +
                        #N_tilt_motors*W_tlt_m +
                            sum([c.W for c in self.components])
                        +W_f)]),
                        b <= bmax,

                        self.htail.W == Wh,
                        self.vtail.W == Wv,
                        Wbatt == Ebatt/Estar_batt, 

                        bmax**2 >= N_lift_motors*self.lift_prop.R**2*pi,
                        bmax >= N_thrust_motors*self.thrust_prop.R,

                        boom.R == R_boom,
                        #boom.l == b,

                        N_booms >= N_lift_motors/(boom.l/self.lift_prop.R)

                        #self.lift_prop.R >= 1*units("ft"),
                        #self.N_lift_motors <= 10,
                      ]

        return constraints, self.components, self.propulsors, self.boom#, self.loading


class FlightState(Model):
    """State variable container, with atmosphere model

    Variables
    ---------
    V                                   [kts]           flight speed
    M                                   [-]             Flight Mach number
    nu                                  [m^2/s]         Kinematic Viscosity
    rho                                 [kg/m^3]        Density
    mu                                  [kg/(m*s)]      Dynamic Viscosity
    a                                   [m/s]           Speed of sound
    T                                   [K]             Temperature
    p                                   [kPa]           Pressure
    alt                                 [ft]            current altitude
    p_sl             101325             [Pa]            Pressure at sea level
    T_sl             288.15             [K]             Temperature at sea level
    L_atm           .0065               [K/m]           Temperature lapse rate
    M_atm           .0065               [kg/mol]        Molar mass of dry air
    R_atm           287                 [m^2/s^2/K]     Gas Constant
    T_s             110.4               [K]             Sutherland Temperature
    C_1             1.458e-6            [kg/(m*s*K^.5)] Sutherland coefficient
    gam_atm         1.4                 [-]             Ratio of specific heats
    x                                   [km]            Downrange distance

    """
    
    def setup(self):
        exec parse_variables(FlightState.__doc__)

    
        with SignomialsEnabled():
            constraints = [
                nu == mu/rho,
                M == V/a,
                (p/p_sl)**(1/5.257) == T/T_sl,
                rho == p/(R_atm*T),
                SignomialEquality(T_sl, T + L_atm*alt),
                SignomialEquality((T + T_s) * mu, C_1*T**1.5),
                a == (gam_atm*R_atm*T)**(.5)
            ]
        return constraints 

class SeaLevelFlightState(Model):
    """State variable container, with atmosphere model

    Variables
    ---------
    V                                   [kts]           flight speed
    M                                   [-]             Flight Mach number
    nu                                  [m^2/s]         Kinematic Viscosity
    rho             1.225               [kg/m^3]        Density
    mu              1.789e-5            [kg/(m*s)]      Dynamic Viscosity
    a               344                 [m/s]           Speed of sound
    T               288.15              [K]             Temperature
    p               101.325             [kPa]           Pressure
    alt             .1                  [ft]            current altitude
    p_sl             101325             [Pa]            Pressure at sea level
    T_sl             288.15             [K]             Temperature at sea level
    L_atm           .0065               [K/m]           Temperature lapse rate
    M_atm           .0065               [kg/mol]        Molar mass of dry air
    R_atm           287                 [m^2/s^2/K]     Gas Constant
    T_s             110.4               [K]             Sutherland Temperature
    C_1             1.458e-6            [kg/(m*s*K^.5)] Sutherland coefficient
    gam_atm         1.4                 [-]             Ratio of specific heats
    x                                   [km]            Downrange distance

    """
    
    def setup(self):
        exec parse_variables(SeaLevelFlightState.__doc__)

    
        with SignomialsEnabled():
            constraints = [
                nu == mu/rho,
                M == V/a,
            ]
        return constraints
                
class Hover(Model):
    """Hover segment performance integration

    Variables
    ---------
    t               .5      [min]           hover time
    E                       [W*hr]             segment energy
    control_mrgn    1.3     [-]             control margin
    """
    def setup(self, aircraft):
        exec parse_variables(Hover.__doc__)

        state = self.state = SeaLevelFlightState()
        perf = self.perf = aircraft.flight_model(aircraft,self.state)
        dedt = perf.Poper

        
        constraints = [perf.lift_prop.T*aircraft.N_lift_motors >= aircraft.Wtotal*control_mrgn,
                       perf.P_rotors >= perf.lift_motor.Pelec*aircraft.N_lift_motors,
                        E == dedt*t,
                        perf.thrust_prop.T*aircraft.N_thrust_motors >= 1e-30*units("lbf"),
                        perf.CD <= 1,
                            ]
            
        return [constraints,
                self.perf,
                self.state,
                ]
class LevelFlight(Model):
    """Level flight segment performance integration

    Variables
    ---------
    R                       [nmi]           range
    V                       [kts]           segment speed
    V_min           100     [kts]           minimum speed
    t                       [min]           segment time
    E                       [W*hr]             segment energy
    """
    def setup(self, aircraft):
        exec parse_variables(LevelFlight.__doc__)

        state = self.state = SeaLevelFlightState()
        perf = self.perf = aircraft.flight_model(aircraft,self.state)
        
        # integrate energy burn
        dedt = perf.Poper
        rho = state.rho

        self.wingg = aircraft.wing.spar.loading(
            aircraft.wing, state, out=False)
        self.wingg.substitutions[self.wingg.Nmax] = 5
        self.wingg.substitutions[self.wingg.Nsafety] = 1.5
        Sh = aircraft.htail.planform.S
        CLhmax = aircraft.htail.planform.CLmax
        Sv = aircraft.vtail.planform.S
        CLvmax = aircraft.vtail.planform.CLmax

        self.htailg = aircraft.htail.spar.loading(
            aircraft.htail, state)
        self.vtailg = aircraft.vtail.spar.loading(
            aircraft.vtail, state)
        #self.fuseload = aircraft.fuselage.loading(aircraft.fuselage,Wtotal)
        
        self.loading = [self.wingg,
                        self.htailg, self.vtailg,
                        #self.fuseload
                        ]

        constraints = [state.V == V,
                        perf.thrust_prop.T*aircraft.N_thrust_motors >= .5*rho*V**2*perf.Sw*perf.CD,
                        perf.lift_prop.T*aircraft.N_lift_motors >= 1e-30*units("lbf"),
                        0.5*rho*V**2*perf.CL*perf.Sw == perf.W,
                        R == t * V,
                        E >= dedt*t,
                        V >= V_min,

                        self.wingg.W >= aircraft.Wtotal,
                        self.htailg.W >= .5*rho*V**2*Sh*CLhmax,
                        self.vtailg.W >= .5*rho*V**2*Sv*CLvmax,

                        perf.P_rotors >= perf.thrust_motor.Pelec*aircraft.N_thrust_motors,
                ]

        return [self.perf,
                self.state,
                constraints,
                self.loading
                ]

class Mission(Model):
    """define mission for aircraft
    Variables
    ---------
    R_total               60             [nmi]           total range
    INT_N_thrust_motors                  [-]             integer thrust motors      
    INT_N_lift_motors                    [-]             integer lift motors
    INT_N_booms                          [-]             number of lift motor support booms
    """
    def setup(self):
        exec parse_variables(Mission.__doc__)
        self.aircraft = aircraft = Aircraft()
        self.cruise = cruise = LevelFlight(self.aircraft)
        self.hover = hover = Hover(self.aircraft)

        self.mission = [self.cruise, 
                        self.hover,
                        ]
        
        constraints = [cruise.perf.W == aircraft.Wtotal,
                        R_total == cruise.R,
                        aircraft.Ebatt >= cruise.E + hover.E,
                        aircraft.N_thrust_motors == INT_N_thrust_motors,
                        aircraft.N_lift_motors == INT_N_lift_motors,
                        aircraft.N_booms == INT_N_booms,
                        ]

        self.cost = self.aircraft.Wtotal
        return self.mission, self.aircraft, constraints

if __name__ == "__main__":
    M = Mission()
    print M
    M.substitutions[M.R_total] = 20
    sol = (M.solve("mosek", iteration_limit = 200, verbosity = 0))
    print sol.summary()