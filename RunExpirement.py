# IMPORTS
## Std Lib
import time
import csv
from statistics import fmean as mean
## Local files
import GlobalConstants
import DebugConstants as db
from MulticastPackingInstance import MulticastPackingInstance
from Solvers.Impls.PureColGenMcpSolver import PureColGenMcpSolver
from Solvers.Impls.JansenZhangMinMaxer import JansenZhangMinMaxer
from Solvers.Impls.PerturbedColGenMcpSolver import PerturbedColGenMcpSolver
from Solvers.Impls.ColGenIPSolver import ColGenIPSolver
from Solvers.Impls.WarmStartColGenMcpSolver import WarmStartColGenMcpSolver

def runExpirement(datetime_str, numReps, paramList, approxLevels, solverTypes, solverList):
    #tableRows = list()
    #tableRows2 = list()
    
    for label,n,m,k,s,d in paramList:
        if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_0:
            print("Parameter Loop: |V|={}, |E|={}, k={}, |S_i|<={}, D = {}".format(n,m,k,s,d))
        
        totalIter = {(apx, solver_id): list() 
                     for apx in approxLevels
                     for solver_id in solverTypes}
        secsPerIt = {(apx, solver_id): list() 
                     for apx in approxLevels
                     for solver_id in solverTypes}
        totalSecs = {(apx, solver_id): list() 
                     for apx in approxLevels
                     for solver_id in solverTypes}
        FinObjVal = {(apx, solver_id): list() 
                     for apx in approxLevels
                     for solver_id in solverTypes}
        FinPotVal = {(apx, solver_id): list() 
                     for apx in approxLevels
                     for solver_id in solverTypes}
        StopFlags = {(apx, solver_id): 0b0000
                     for apx in approxLevels
                     for solver_id in solverTypes}

        for i in range(numReps):
            if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_0:
                print("\t Instance Repitition: {}".format(i))
            instance =  MulticastPackingInstance(n,m,k,s,d)

            for apx in approxLevels:
                for solver_id in solverTypes:
                    if solver_id == GlobalConstants.COLGEN_ID:
                        solver = PureColGenMcpSolver(instance=instance, block_approx=apx)
                    elif solver_id == GlobalConstants.JZ2008_ID:
                        solver = JansenZhangMinMaxer(instance=instance, block_approx=apx)
                    elif solver_id == GlobalConstants.PERTCG_ID:
                        solver = PerturbedColGenMcpSolver(instance=instance, block_approx=apx)
                    elif solver_id == GlobalConstants.FULLIP_ID:
                        solver = ColGenIPSolver(instance=instance, block_approx=apx)
                    elif solver_id == GlobalConstants.WARMST_ID:
                        solver = WarmStartColGenMcpSolver(instance=instance, block_approx=apx)

                    if db.DEBUG_LEVEL >= db.DEBUG_LEVEL_THEORY_0:
                        print("\t\t Algo: {} \t Block Approx: {}".format(solver_id, apx))
                    try:
                        timer = list()
                        total_time = 0
                        while(not solver.stop_flag):
                            prev_time = time.perf_counter()
                            solver.perform_iteration()
                            new_time = time.perf_counter()
                            timer.append(new_time - prev_time)
                            total_time += new_time - prev_time
                            if total_time > GlobalConstants.MAX_TIME:
                                solver.stop_flag |= GlobalConstants.STOP_FLAG_TIMEOUT

                        x = solver.solution[solver.iteration]
                        totalIter[apx, solver_id].append(solver.iteration)
                        secsPerIt[apx, solver_id].append(mean(timer))
                        totalSecs[apx, solver_id].append(sum(timer))
                        FinObjVal[apx, solver_id].append(solver.lamb(x))
                        FinPotVal[apx, solver_id].append(solver.phi(x, solver.t))
                        StopFlags[apx, solver_id] |= solver.stop_flag

                        newRow = [
                            label,
                            solver_id, apx,
                            n,m,k,s,d,
                            solver.lamb(x),
                            solver.phi(x, solver.t),
                            solver.iteration,
                            sum(timer),
                            mean(timer),
                            "{:05b}".format(solver.stop_flag)
                        ]
                        #solverList.append(solver) # I Think this causes the program to use too much RAM
                        with open("outputs/{}.csv".format(datetime_str),'a') as allFile:
                            allWriter = csv.writer(allFile)
                            allWriter.writerow(newRow)
                        
                    except Exception as e:
                        print(repr(e))

#         for apx in approxLevels:
#             for solver_id in solverTypes:
#                 try:
#                     avgTotalIter = mean(totalIter[apx, solver_id])
#                     avgSecsPerIt = mean(secsPerIt[apx, solver_id])
#                     avgtotalSecs = mean(totalSecs[apx, solver_id])
#                     avgFinObjVal = mean(FinObjVal[apx, solver_id])
#                     avgFinPotVal = mean(FinPotVal[apx, solver_id])
#                     newRow = [label, solver_id, apx, n, m, k, s, d,
#                           avgFinObjVal, avgFinPotVal, avgTotalIter, 
#                           avgtotalSecs, avgSecsPerIt, 
#                           "{:05b}".format(StopFlags[apx, solver_id])]


#     with open("outputs/summaries/summary_{}.csv".format(datetime_str),'a') as sumFile:
#         sumWriter = csv.writer(sumFile)
#         sumWriter.writerow(newRow)
        
    return solverList