"""Now that you know everything about data and visualization, let's get started with optimization!
Optimeed provides high-level interface to perform optimization with visualization and data storage.
The Wiki gives more details about the optimization. To get started, you need the following key ingredients:

    - A device that contains the variables to be optimized ("Device") and other parameters you would like to save
    - A list of optimization variables ("OptimizationVariable")
    - An evaluation function ("Characterization")
    - One or more objective functions ("Objectives")
    - (optional) Eventual constraints ("Constraints")
    - An optimization algorithm ("Optimization Algorithm")
    - Something that will fill the "Device" object with the optimization variables coming from the optimization algorithm. ("MathsToPhysics")
        Don't get scared with this one, if you do not know how it can be useful, the proposition by default works perfectly fine.
    - Something that will link all the blocks together ("Optimizer")
"""


# These are what we need for the optimization
from optimeed.optimize.optiAlgorithms import MultiObjective_GA as OptimizationAlgorithm
from optimeed.optimize import Real_OptimizationVariable, InterfaceObjCons, InterfaceCharacterization, OptiHistoric
from optimeed.optimize.optimizer import OptimizerSettings, run_optimization

import csv


# These are the high-level visualization tools
from optimeed.visualize.displayOptimization import OptimizationDisplayer
from optimeed.visualize import Onclick_representDevice, Represent_brut_attributes, start_qt_mainloop
from optimeed.core.evaluators import PermanentMultiprocessEvaluator, MultiprocessEvaluator, Evaluator
import time
import os

from LinearMotor import *

class Device:
    """Define the Device to optimize."""
    tau_k: float  # Type hinted -> will be automatically saved
    e: float
    hm: float
    ha: float
    Lz: float
    lq: float
    


    def __init__(self):
        self.tau_k = 0
        self.e = 0
        self.hm = 0
        self.ha = 0
        self.Lz = 0
        self.lq = 0

class Characterization(InterfaceCharacterization):
    """Define the Characterization scheme. In this case nothing is performed,
     but this is typically where model code will be executed and results saved inside 'theDevice'.
     The arguments time_initialization and time_evaluation are there to mimic the evaluation time of the model.
     They will be use to compare different evaluators
     """
    def __init__(self, time_initialization=1, time_evaluation=0.01):
        self.initialized = False
        self.time_initialization = time_initialization
        self.time_evaluation = time_evaluation

    def compute(self, thedevice):
        if not self.initialized:
            print("Characterization is initializing -> {} second penalty".format(self.time_initialization))
            time.sleep(self.time_initialization)
            self.initialized = True
        print("Process {} is active. The id of the characterization: {}".format(hex(os.getpid()), hex(id(self))))
        time.sleep(self.time_evaluation)




class MyObjective1(InterfaceObjCons):
    """First objective function (to be minimized)"""
    def compute(self, thedevice):
        p2 = LinearMotor(thedevice.tau_k, thedevice.e, thedevice.hm, thedevice.ha, thedevice.Lz, thedevice.lq)
        THD = p2.get_THD()
        return THD
    
class MyObjective2(InterfaceObjCons):
    """Second objective function (to be minimized)"""
    def compute(self, thedevice):
        p1 = LinearMotor(thedevice.tau_k, thedevice.e, thedevice.hm, thedevice.ha, thedevice.Lz, thedevice.lq)
        F_active = p1.get_F_active()
        return -F_active

class MyObjective3(InterfaceObjCons):
    """Third objective function (to be minimized)"""
    def compute(self, thedevice):
        p1 = LinearMotor(thedevice.tau_k, thedevice.e, thedevice.hm, thedevice.ha, thedevice.Lz, thedevice.lq)
        Time = p1.get_time()
        return Time


class MyConstraint1(InterfaceObjCons):
    """Constraints, that needs to be <= 0"""
    def compute(self, thedevice):
        lq = thedevice.lq
        e = thedevice.e
        return lq - e/2
    
class MyConstraint2(InterfaceObjCons):
    """Constraints, that needs to be <= 0"""
    def compute(self, thedevice):
        tau_k = thedevice.tau_k
        e = thedevice.e
        return e - tau_k


if __name__ == "__main__":  # This line is necessary to spawn new processes
    """Start the main code. Instantiate previously defined classes."""
    theDevice = Device()
    theAlgo = OptimizationAlgorithm()
    theAlgo.set_option(theAlgo.OPTI_ALGORITHM, "NSGAII")  # You can change the algorithm if you need ;)

    theCharacterization = Characterization(time_initialization=1, time_evaluation=0.01)

    """Variable to be optimized"""
    optimizationVariables = list()
    """
    optimizationVariables.append(Real_OptimizationVariable('tau_k', 1e-3, 30e-3))
    optimizationVariables.append(Real_OptimizationVariable('e', 1e-3, 10e-3))
    optimizationVariables.append(Real_OptimizationVariable('hm', 1e-3, 20e-3))
    optimizationVariables.append(Real_OptimizationVariable('ha', 1e-3, 15e-3))
    optimizationVariables.append(Real_OptimizationVariable('Lz', 1e-3, 50e-3))
    optimizationVariables.append(Real_OptimizationVariable('lq', 1e-3, 10e-3))
    """
    
    optimizationVariables.append(Real_OptimizationVariable('tau_k', 1e-3, 25e-3))
    optimizationVariables.append(Real_OptimizationVariable('e', 1e-3, 20e-3))
    optimizationVariables.append(Real_OptimizationVariable('hm', 1e-3, 10e-3))
    optimizationVariables.append(Real_OptimizationVariable('ha', 1e-3, 5e-3))
    optimizationVariables.append(Real_OptimizationVariable('Lz', 1e-3, 50e-3))
    optimizationVariables.append(Real_OptimizationVariable('lq', 1e-3, 5e-3))

    """Objective and constraints"""
    listOfObjectives = [MyObjective1(), MyObjective3()]
    listOfConstraints = [MyConstraint1(), MyConstraint2()]

    """Set the optimizer"""
    theOptiSettings = OptimizerSettings(theDevice, listOfObjectives, listOfConstraints, optimizationVariables,
                                        theOptimizationAlgorithm=theAlgo, theCharacterization=theCharacterization)

    """The logger (to automatically save the points)"""
    theOptiHistoric = OptiHistoric(optiname="opti", autosave_timer=10, autosave=True, create_new_directory=True)

    """The evaluator. It is used to manage the evaluations of the parameters from the optimization algorithm, and the behaviour in parallel run.
    The first evaluator does not allow parallel run. You can see it with the constant process and id characterization. It is the recommended evaluator if you do not need parallelism. 
    The second evaluator allow parallel run. It is the recommended evaluator for parallelism. Notice each process getting its own id, but the characterization is forked at each call.
    In other words, if initialization of the models is long to execute (for instance, opening matlab), then you will pay the full price each time.
    The third evaluator also allows parallel run, but the characterizations are forked only once -> no extra penalty on initialization afterwards.
    """
    theEvaluator = Evaluator(theOptiSettings)  # First evaluator -> One initialization, then proceeds in same thread
    # theEvaluator = MultiprocessEvaluator(theOptiSettings, number_of_cores=2)  # Second evaluator -> Parallel run, initializes Charac() at each run.
    # theEvaluator = PermanentMultiprocessEvaluator(theOptiSettings, number_of_cores=2)  # Third evaluator -> Parallel run, initializes Charac() at startup.

    """Start the optimization"""
    max_opti_time_sec = 300

    display_opti = True

    if display_opti:  # Display real-time graphs
        optiDisplayer = OptimizationDisplayer(theOptiSettings, theOptiHistoric, light_background=True)
        _, theDataLink, _ = optiDisplayer.generate_optimizationGraphs()

        # Here we set the actions on click.
        theActionsOnClick = list()
        theActionsOnClick.append(Onclick_representDevice(theDataLink, [Represent_brut_attributes()]))
        optiDisplayer.set_actionsOnClick(theActionsOnClick)

        start_time = time.time()  # Mark the starting time of the optimization


        resultsOpti, convergence = optiDisplayer.launch_optimization([theOptiSettings, theOptiHistoric],
                                                                     {"max_opti_time_sec": max_opti_time_sec, "evaluator": theEvaluator},
                                                                     refresh_time=0.1, max_nb_points_convergence=None)  # Refresh the graphs each nth seconds

        # Calculate remaining time
        elapsed_time = time.time() - start_time
        remaining_time = max_opti_time_sec - elapsed_time
        print(f"Remaining time: {remaining_time} seconds")

    else:  # Otherwise just focus on results ... That can be helpful if you are confident the optimizations will converge and you need to launch several optimizations.
        resultsOpti, convergence = run_optimization(theOptiSettings, theOptiHistoric, max_opti_time_sec=max_opti_time_sec, evaluator=theEvaluator)

    """Gather results"""
    # Pro hint: you would probably never work with these next few lines of code, instead you would move to the next tutorial
    # to retrieve the results from the automatically saved files.
    print("Best individuals :")
    for device in resultsOpti:
        print("tau_k : {} \t e : {} \t hm : {} \t ha : {} \t Lz : {} \t lq : {}". format(device.tau_k, device.e, device.hm, device.ha, device.Lz, device.lq))
        p = LinearMotor(device.tau_k, device.e, device.hm, device.ha, device.Lz, device.lq)
        print("F_active = {}".format(p.get_F_active()))
        print("THD = {}".format(p.get_THD()))
            

    # Open a file named OptimalMachines.csv in write mode
    with open("OptimalMachines.csv", "w", newline='') as csvfile:
        fieldnames = ['Machine', 'tau_k', 'e', 'hm', 'ha', 'Lz', 'lq', 'F_active', 'THD', 'mass', 'time']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for idx, device in enumerate(resultsOpti, start=1):
            p = LinearMotor(device.tau_k, device.e, device.hm, device.ha, device.Lz, device.lq)
            writer.writerow({'Machine': idx,
                            'tau_k': device.tau_k,
                            'e': device.e,
                            'hm': device.hm,
                            'ha': device.ha,
                            'Lz': device.Lz,
                            'lq': device.lq,
                            'F_active': p.get_F_active(),
                            'THD': p.get_THD(),
                            'mass': p.get_totalMass(),
                            'time': p.get_time()})


    if display_opti:
        start_qt_mainloop()  # To keep windows alive

    """Note that the results are automatically saved if KWARGS_OPTIHISTO autosaved=True.
    In this case, optimization folder is automatically generated in Workspace/optiX. It contains five files:
    -> autosaved: contains all the devices evaluated during the optimization
    -> logopti: contains all the information relating to the optimization itself: objectives, constraints, evaluation time.
    -> opticonvergence: contains all the information relative to the convergence of the optimization (saved only at the end)
    -> results: all the best devices as decided by the optimization algorithm
    -> optimization_parameters: the class OptimizationParameters that can be reloaded using SingleObjectSaveLoad.load 
    -> summary.html: a summary of the optimization problem
    See other tutorials on how to save/load these information.
    """
