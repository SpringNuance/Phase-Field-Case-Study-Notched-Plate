from driverConstants import *
from driverStandard import StandardAnalysis
import driverUtils, sys
options = {
    'SIMExt':'.sim',
    'aba_timing':1,
    'adaptiveMesh':OFF,
    'ams':OFF,
    'analysisType':STANDARD,
    'applicationName':'analysis',
    'aqua':OFF,
    'beamSectGen':OFF,
    'biorid':OFF,
    'cavityTypes':[],
    'cavparallel':OFF,
    'complexFrequency':OFF,
    'contact':OFF,
    'cosimulation':OFF,
    'coupledProcedure':ON,
    'cpus':1,
    'cpusValidated':True,
    'cse':OFF,
    'cyclicSymmetryModel':OFF,
    'directCyclic':OFF,
    'direct_solver':DMP,
    'direct_solver_use_default':'1',
    'dsa':OFF,
    'dynStepSenseAdj':OFF,
    'dynamic':OFF,
    'exciteBodies':[],
    'externalField':OFF,
    'externalFieldCSEAux':OFF,
    'externalFieldExtList':['.sim', '.SMAManifest'],
    'externalFieldFiles':[],
    'externalFieldSimReader':None,
    'fieldImport':OFF,
    'filPrt':[],
    'fils':[],
    'finitesliding':OFF,
    'flexibleBody':OFF,
    'foundation':OFF,
    'freqSimReq':OFF,
    'geostatic':OFF,
    'heatTransfer':OFF,
    'impJobExpVars':{},
    'importJobList':[],
    'importSim':OFF,
    'importer':OFF,
    'importerParts':OFF,
    'includes':[],
    'inertiarelief':OFF,
    'initialConditionsFile':OFF,
    'input':'phase_field_3D_millimeter_UEL',
    'inputFormat':INP,
    'interactive':None,
    'interpolExtList':['.odb', '.sim', '.SMAManifest'],
    'jleJobExcluded':False,
    'job':'phase_field_3D_millimeter_UEL',
    'keyword_licenses':[],
    'lanczos':OFF,
    'libs':[],
    'magnetostatic':OFF,
    'massDiffusion':OFF,
    'materialresponse':OFF,
    'modifiedTet':OFF,
    'moldflowFiles':[],
    'moldflowMaterial':OFF,
    'mp_mode':THREADS,
    'mp_mode_requested':MPI,
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'noMultiHostElemLoop':[],
    'no_domain_check':1,
    'numLoadCases':0,
    'onestepinverse':OFF,
    'outputKeywords':ON,
    'parameterized':OFF,
    'partsAndAssemblies':ON,
    'parval':OFF,
    'pgdHeatTransfer':OFF,
    'postOutput':OFF,
    'preDecomposition':ON,
    'randomResponse':OFF,
    'restart':OFF,
    'restartEndStep':OFF,
    'restartIncrement':0,
    'restartStep':0,
    'restartWrite':OFF,
    'rezone':OFF,
    'runCalculator':OFF,
    'simpackBodies':[],
    'soils':OFF,
    'soliter':OFF,
    'solverTypes':['DIRECT'],
    'ssd':OFF,
    'ssdCheckOK':False,
    'standard_parallel':ALL,
    'staticNonlinear':OFF,
    'steadyStateTransport':OFF,
    'step':ON,
    'stepSenseAdj':OFF,
    'stressExtList':['.odb', '.sim', '.SMAManifest'],
    'subGen':OFF,
    'subGenResidual':OFF,
    'submodel':OFF,
    'substrLibDefs':OFF,
    'substructure':OFF,
    'symModGenModel':OFF,
    'symmetricModelGeneration':OFF,
    'symmetricResultsTransfer':OFF,
    'tempNoInterpolExtList':['.fil', '.odb', '.sim', '.SMAManifest'],
    'thermal':OFF,
    'tmpdir':'C:\\Users\\nguyenb5\\AppData\\Local\\Temp',
    'tracer':OFF,
    'transientSensitivity':OFF,
    'unfold_param':OFF,
    'unsymm':OFF,
    'user':'phase_field_3D.f90',
    'verbose':1,
    'visco':OFF,
    'xplSelect':OFF,
}
analysis = StandardAnalysis(options)
status = analysis.run()
sys.exit(status)
