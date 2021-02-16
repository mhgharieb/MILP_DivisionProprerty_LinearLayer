


from gurobipy import *
import time
from AESLike import *
from SolveModel import *
from constants import Matrix
from ShurnkMatrixModel import *
import multiprocessing as mp
import shutil
import os


machine = os.uname()[1]

_logDir = '/tmp/aes_4_round/'

if not os.path.exists(_logDir):
    try:
        os.mkdir(_logDir)
    except OSError:
        print ("Creation of the directory %s failed" % _logDir)
        exit()
    else:
        print ("Successfully created the directory %s " % _logDir)
    




# Gurobi general parapmeters
ISTIMELIMIT = 0
TIMELIMIT = 3600
ISThREADS = 1
ISNodefileDir = 1
NodefileStart = 0.5
UseOptimize = 0

def printResult(inputDP, balanceBits):
    activebits = inputDP.count('1')
    fileobj = open(filename_result, "a", buffering=1)
    print('# of rounds: %i' % Round + " , activebits:%i " % activebits + ", Input DP: " + inputDP)
    print("balanced bits: " + str(len(balanceBits)) + ' ' + str(balanceBits))
    fileobj.write('# of rounds: %i' % Round + " , activebits:%i " % activebits + ", Input DP: " + inputDP + "\n")
    fileobj.write("balanced bits: " + str(len(balanceBits)) + ' ' + str(balanceBits) + "\n")
    fileobj.close()


def WriteModel(constraints, variables, objective, filename_model):
    file_co = open(filename_model, 'w+')
    if UseOptimize:
        file_co.write('Minimize\n')
        file_co.write(' + '.join(objective) + '\n')
    file_co.write('Subject To\n')
    for enq in constraints:
        file_co.write(enq + '\n')
    file_co.write('Binary\n')
    for v in variables:
        file_co.write(v + '\n')
    file_co.close()


def createModel(Rounds, inputDP_column, X_Output, filename_model):
    variables = []
    constraints = []
    InputOutputVariablesLinearLayer = []
    InputVariables = []
    OutputVariables = []
    

    # round 0 MixColumn
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    Z = ShiftRow(CreateStateVariable('y', blocksize//wordSize, wordSize, 0))
    Xout = CreateStateVariable('x', blocksize//wordSize, wordSize, 1)
    Z_column = [Z[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    Xout_column = [Xout[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    for ColIndex in range(len(inputDP_column)):
        if len(inputDP_column[ColIndex]) != ColSize:
            currentColZ = Z_column[ColIndex][:-1]
            tmpVariables , tmpConstraints, tmpInputOutputVariablesLinearLayer = ConstraintsByShrunkMixColumnNew(currentColZ, Xout_column[ColIndex], 0, ColIndex)    
            variables += words2bits(currentColZ)
            variables += tmpVariables
            constraints += tmpConstraints
            InputOutputVariablesLinearLayer += [tmpInputOutputVariablesLinearLayer]
            InputVariables += words2bits(currentColZ)   
        else:
            InputVariables += words2bits(Xout_column[ColIndex])   
    
    constraints = ['%s = %s' %(a, b) for a, b in zip(InputVariables, ''.join(sum(inputDP_column,[])))] + constraints
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    for r in range(1, Rounds - 1):
        tmpVariables, tmpConstraints, tmpInputOutputVariablesLinearLayer = ConstraintByRound(r)
        variables += tmpVariables
        constraints += tmpConstraints
        InputOutputVariablesLinearLayer += tmpInputOutputVariablesLinearLayer
        
    
    # round r - 1
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    Xin = CreateStateVariable('x', blocksize//wordSize, wordSize, Rounds - 1)
    Y = CreateStateVariable('y', blocksize//wordSize, wordSize, Rounds - 1)
    Xout = CreateStateVariable('x', blocksize//wordSize, wordSize, Rounds)

    variables += words2bits(Xin)
    variables += words2bits(Y)
    
    ## Sbox layer
    constraints += ConstraintsBySboxLayer(Xin, Y)
    ## Shift Row
    Z = ShiftRow(Y)
    ## MixColumn 
    Z_column = [Z[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    Xout_column = [Xout[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    XorOutput = []
    ZeroOutput = []
    for ColIndex in range(blocksize//(wordSize * ColSize)):
        currentZcol = words2bits(Z_column[ColIndex])
        currentXoutcol = words2bits(Xout_column[ColIndex]) 
        if X_Output in currentXoutcol:
            outputCheckBitIndx = currentXoutcol.index(X_Output)
            outputBit = X_Output
            XorOutput = list(compress(currentZcol, Matrix[outputCheckBitIndx]))
            ZeroOutput += [a for a in currentZcol if a not in XorOutput]
        else:
            ZeroOutput += currentZcol   
    
    variables += [outputBit]
    constraints.append( '%s - %s = 0' %(' + '.join(XorOutput), outputBit))
    constraints.append('%s = 1' %outputBit)
    constraints += ['%s = 0' %a for a in ZeroOutput]

    OutputVariables = [outputBit]
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------

    WriteModel(constraints, variables, OutputVariables, filename_model)

    model = read(filename_model)
    model._InputVariables = getVarsByName(model, InputVariables)
    model._OutputVariables = getVarsByName(model, OutputVariables)
    model._InputOutputVariablesLinearLayer = [(getVarsByName(model, X), getVarsByName(model, Y)) for X, Y in InputOutputVariablesLinearLayer]
    model._XorOutput = getVarsByName(model, OutputVariables)

    return model






def ConstraintByRound(r):
    NBytesState = blocksize//wordSize
    Xin = CreateStateVariable('x', NBytesState, wordSize, r)
    Y = CreateStateVariable('y', NBytesState, wordSize, r)
    Xout = CreateStateVariable('x', NBytesState, wordSize, r + 1)

    Variables = []
    Constraints = []

    Variables += sum(Xin, [])  # convert 2d to 1d
    Variables += sum(Y, [])

    tmpConstraints = ConstraintsBySboxLayer(Xin, Y)
    Constraints += tmpConstraints

    Z = ShiftRow(Y)

    tmpVariables, tmpConstraints, InputOutputVariablesLinearLayer = ConstraintsByLinearLayer(Z, Xout, ColSize, r)
    Variables += tmpVariables
    Constraints += tmpConstraints

    return (Variables, Constraints, InputOutputVariablesLinearLayer)



def findAllUNKNOWN(m):
    Z_UNKNOWN = m._XorOutput[getVarsValue(m._XorOutput).index(1)].getAttr('VarName')
    Y = CreateStateVariable('y', blocksize//wordSize, wordSize, m._Rounds - 1)
    Xout = CreateStateVariable('x', blocksize//wordSize, wordSize, m._Rounds)
    Z = ShiftRow(Y)
    Z_column = [Z[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    Xout_column = [Xout[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    for ColIndex in range(len(Z_column)):
        currentZcol = words2bits(Z_column[ColIndex])
        currentXoutcol = words2bits(Xout_column[ColIndex]) 
        if Z_UNKNOWN in currentZcol:
            Z_UNKNOWN_Indx = currentZcol.index(Z_UNKNOWN)
            allUNKNOWN = list(compress(currentXoutcol, MatrixTranspose[Z_UNKNOWN_Indx])) # column in Matrix
            Xout_serialized = words2bits(Xout)
            _UNKNOWN_mask = [0] * blocksize
            for v in allUNKNOWN:
                _UNKNOWN_mask[Xout_serialized.index(v)] = 1
            return ([Z_UNKNOWN] + allUNKNOWN, _UNKNOWN_mask)



def deleteNodefileDir():
    shutil.rmtree(NodefileDir, True)

def SubRun(Round, shrunkColIndex, outputBit, Y_Output):


    ''''
        It is AES encryption using parameters

    Round: AK -> SB -> SR -> MC
    where:
        SB: use SboxInqs in constants.py
        SR: use ShiftRowIndex in constants.py
        MC: use Matrix in constants.py
    '''
    
    
    
    filename_model = "tmp/shrunkColIndex_%i_outputCheckBit_%i_%s_" %(shrunkColIndex, outputCheckBit, outputBit) + Blockcipher + "_" + str(Round) + "r" + "_%s" %(str(time.time())) +  ".lp"
    filename_result = "tmp/shrunkColIndex_%i_outputCheckBit_%i_%s_" %(shrunkColIndex, outputCheckBit, outputBit) + Blockcipher + "_Result_%s_" %("optimize" if UseOptimize else "bitbybit") + str(Round) + "r" + "_%s_" %machine.split('.')[0] + "_%s" %(str(time.time())) + ".txt"


    fileobj = open(filename_result, "a", buffering=1)
    
    
    inputDP_column = []
    for i in range(blocksize//(wordSize*ColSize)):
        if i ==  shrunkColIndex:
            inputDP_column.append(['1' * wordSize] * (ColSize - 1))
        else:
            inputDP_column.append(['0' * wordSize] * ColSize)
    
    
    m = createModel(Round, inputDP_column, outputBit, filename_model)
    m._Rounds = Round
    m._blocksize = blocksize
    m._Blockcipher = Blockcipher
    m._filenameResult = filename_result
    m._InputCheckBit = shrunkColIndex
    m._DetailedOutput = 1
    m.setParam("LogToConsole", 0)
    if ISTIMELIMIT:
        m.setParam("TIME_LIMIT", TIMELIMIT)
        print("Time limit: %i seconds" %TIMELIMIT, flush=True)
    if ISThREADS:
        m.setParam("Threads", Threads)
    if ISNodefileDir:
        m.setParam("NodefileDir", NodefileDir[:-1])
        m.setParam("NodefileStart", NodefileStart)

    fileobj.write("shrunkColIndex_%i_outputCheckBit_%i_%s_" %(shrunkColIndex, outputCheckBit, Y_Output) + "\n")
        
    logFile = logDir +"shrunkColIndex_%i_outputCheckBit_%i_%s_" %(shrunkColIndex, outputCheckBit, Y_Output) +  machine.split('.')[0] + "_GurobiLog_" + m._Blockcipher + "_Intgral_%02i_%s" % (Round, str(time.time()))
    m.setParam("LogFile", logFile)
    fileobj.write("Log: %s\n" % logFile)
    print("Log: %s" % logFile, flush=True)
    fileobj.write("shrunkColIndex: %i\n" %shrunkColIndex)
    balanceBits = SolveModelBitbyBitSingleBit(m)
    print(str(balanceBits), flush=True)
    ##### read log #########
    a_file = open(logFile, "r")
    lines = a_file.readlines()
    fileobj.write("log part \n")
    fileobj.write("".join(lines[-20:]))

    if len(balanceBits) == 0:
        allUNKNOWN, UNKNOWN_mask_currentRun = findAllUNKNOWN(m)
        for v in allUNKNOWN:
            fileobj.write("%03i" %(0) + " (" + v + ") " + " -- UNKNOWN, run time: %f Seconds \n" %(0))
    else:
        UNKNOWN_mask_currentRun = []
    fileobj.close()
    return (balanceBits, UNKNOWN_mask_currentRun, filename_result, shrunkColIndex)


def SubRunCallback(Result):
    print("Call SubRun:%s" %str(Result), flush=True)
    balanceBits, UNKNOWN_mask_part, filename_result, shrunkColIndex = Result
    fileResults.append(filename_result)
    if len(balanceBits) == 0:
        for i in range(blocksize):
            UNKNOWN_mask[i] |= UNKNOWN_mask_part[i]


        file_UNKNOWN_mask = open("txt/UNKNOWN_mask_print_%03i" %shrunkColIndex, "a", buffering=1)
        file_UNKNOWN_mask.write(''.join(list(map(str,UNKNOWN_mask))) + "\n")
        file_UNKNOWN_mask.close()


def NumberOfSubRun(Rounds, outputCheckBit):
    # round r - 1
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    Y = CreateStateVariable('y', blocksize//wordSize, wordSize, Rounds - 1)
    Xout = CreateStateVariable('x', blocksize//wordSize, wordSize, Rounds)

    ## Shift Row
    Z = ShiftRow(Y)
    
    Z_column = [Z[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    Xout_column = [Xout[i:i+ColSize] for i in range(0, blocksize//wordSize, ColSize)]
    
    ColIndex = outputCheckBit // (wordSize * ColSize)
    outputCheckBit_in_Col = outputCheckBit % (wordSize * ColSize)
    
    currentZcol = words2bits(Z_column[ColIndex])
    currentXoutcol = words2bits(Xout_column[ColIndex]) 
    
    outputBit = currentXoutcol[outputCheckBit_in_Col]
    XorOutput = list(compress(currentZcol, Matrix[outputCheckBit_in_Col]))

    return((outputBit, XorOutput))
        





if __name__ == "__main__":
    

    

    outputCheckBit_start = 0
    outputCheckBit_end = 128
    InputCheckBit_list = [0]
    
    
    

    global blocksize
    global wordSize
    global ColSize
    global UNKNOWN_mask
    global fileResults
    global logDir
    global NodefileDir
    global Blockcipher
    global Threads

    Round = 4
    blocksize = 128
    wordSize = 8
    ColSize = 4
    
    Blockcipher = "AESInverse"

    for shrunkColIndex in InputCheckBit_list:
        time_start_InputCheckBit = time.time()

        ######## Create log dir ###############

        logDir = _logDir + "%s_%ir_shrunkColIndex_%03i/" %(Blockcipher, Round, shrunkColIndex)
        if not os.path.exists(logDir):
            try:
                os.mkdir(logDir)
            except OSError:
                print ("Creation of the directory %s failed" % logDir)
                exit()
            else:
                print ("Successfully created the directory %s " % logDir)
        
        ######## Create log dir ###############
        
        print("=================================================================================================")
        print("shrunkColIndex: %i\n" %shrunkColIndex)
        

        
        filename_TotalResult = "txt/Total_%i.txt" %shrunkColIndex 
        fileTotal = open(filename_TotalResult, "a", buffering=1)
        fileTotal.write("shrunkColIndex: %i\n" %shrunkColIndex)



        try:
            file_UNKNOWN_mask = open("txt/UNKNOWN_mask_print_%03i" %shrunkColIndex, "r")
            last = file_UNKNOWN_mask.read().splitlines()[-1]
            file_UNKNOWN_mask.close()
            UNKNOWN_mask = list(map(int, list(last)))
            if len(UNKNOWN_mask) != blocksize:
                UNKNOWN_mask = [0] * blocksize
        except:
            UNKNOWN_mask = [0] * blocksize
        





        
        for outputCheckBit in range(outputCheckBit_start, outputCheckBit_end):
            time_start = time.time()


            outputBit, XorOutput = NumberOfSubRun(Round, outputCheckBit)
            numSubRun = len(XorOutput)
            numCPU = mp.cpu_count()
            Threads = numCPU

            fileTotal.write("**********************************\n")
            print("**********************************")
            fileTotal.write("Run %i\n" %outputCheckBit)
            print("Run %i" %outputCheckBit)
            if UNKNOWN_mask[outputCheckBit] == 1:
                fileTotal.write("UNKNOWN\n")  
                print("UNKNOWN")  
                continue  



            NodefileDir = logDir + "NodefileDir/"

            deleteNodefileDir()

            if not os.path.exists(NodefileDir):
                try:
                    os.mkdir(NodefileDir)
                except OSError:
                    print ("Creation of the directory %s failed" % NodefileDir)
                    exit()
                else:
                    print ("Successfully created the directory %s " % NodefileDir)


            fileResults = []
            UNKNOWN_mask_before = list(UNKNOWN_mask)
            
            
            print(SubRun(Round, shrunkColIndex, outputBit, XorOutput[0]))
            
            if UNKNOWN_mask_before == UNKNOWN_mask:
                # no division trail
                Xout = words2bits(CreateStateVariable('x', blocksize//wordSize, wordSize, Round))
                fileTotal.write("%03i" %outputCheckBit + " (" +  Xout[outputCheckBit] + ") " + " -- Balanced\n")
                print("%03i" %outputCheckBit + " (" +  Xout[outputCheckBit] + ") " + " -- Balanced")

            for filename_result in fileResults:
                f = open(filename_result, 'r')
                fileTotal.write(f.read())
                f.close()



            time_end = time.time()
            fileTotal.write("run time: %f Seconds \n" %(time_end - time_start) )
            print("run time: %f Seconds " %(time_end - time_start) )

        deleteNodefileDir()
            
        fileTotal.write("=================================================================================================")
        fileTotal.write("Total run time: %f Seconds \n" %(time_end - time_start_InputCheckBit) )
        print("Total run time: %f Seconds " %(time_end - time_start_InputCheckBit) )
        fileTotal.close()