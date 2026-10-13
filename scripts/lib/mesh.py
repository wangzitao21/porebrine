from pathlib import Path
import numpy as np


def head(name, cls='dictionary'):
    return f'FoamFile {{ version 2.0; format ascii; class {cls}; object {name}; }}\n'


def make_case(path, epsilon, spacing=5e-6, end=20., dt=.5, ports=25, *, thickness=None, solver_profile='strict', component_scheme='upwind'):
    epsilon=np.asarray(epsilon,dtype=float)
    thickness=spacing if thickness is None else thickness
    if epsilon.ndim!=2 or min(epsilon.shape)<2 or not np.isfinite(epsilon).all() or np.any((epsilon<=0)|(epsilon>1)):
        raise ValueError('Porosity must be a finite 2D array in (0, 1], with at least two cells per direction')
    if not all(np.isfinite(v) and v>0 for v in [spacing,thickness,end,dt,ports]):
        raise ValueError('Mesh dimensions, duration, step and port size must be positive and finite')
    if solver_profile not in ('strict','balanced'):raise ValueError('Unknown linear solver profile')
    component_schemes={'upwind':'Gauss upwind','limitedLinear':'Gauss limitedLinear 1'}
    if component_scheme not in component_schemes:raise ValueError('Unknown component convection scheme')
    path=Path(path)
    for folder in ['0','system','constant']:
        (path/folder).mkdir(parents=True,exist_ok=True)
    nz,nx=epsilon.shape;dx=dz=spacing;length=nx*dx;height=nz*dz

    def put(name,text):
        (path/name).write_text(text)

    def field(name,data,bc=None,dimensions='[0 0 0 0 0 0 0]',cls='volScalarField'):
        if bc is None:
            bc='inlet {type zeroGradient;} outlet {type zeroGradient;} walls {type zeroGradient;}'
        if isinstance(data,np.ndarray):
            data=f'nonuniform List<scalar>\n{data.size}\n(\n'+'\n'.join(f'{x:.17g}' for x in data.ravel())+'\n)'
        elif not isinstance(data,str):
            data=f'uniform {data:.17g}'
        put('0/'+name,head(name,cls)+f'dimensions {dimensions};\ninternalField {data};\nboundaryField {{ {bc} frontAndBack {{type empty;}} }}\n')

    field('epsilon',epsilon)
    field('p',0.,'inlet {type zeroGradient;} outlet {type fixedValue;value uniform 0;} walls {type zeroGradient;}','[1 -1 -2 0 0 0 0]')
    field('U','uniform (0 0 0)','inlet {type fixedValue;value uniform (0 0 0);} outlet {type zeroGradient;} walls {type fixedValue;value uniform (0 0 0);}','[0 1 -1 0 0 0 0]','volVectorField')
    put('system/blockMeshDict',head('blockMeshDict')+f'''convertToMeters 1;
vertices ((0 0 0) ({length} 0 0) ({length} {height} 0) (0 {height} 0) (0 0 {thickness}) ({length} 0 {thickness}) ({length} {height} {thickness}) (0 {height} {thickness}));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {nz} 1) simpleGrading (1 1 1)); edges ();
boundary (inlet {{type patch; faces ((0 4 7 3));}} outlet {{type patch; faces ((1 2 6 5));}}
walls {{type wall; faces ((0 1 5 4) (3 7 6 2));}} frontAndBack {{type empty; faces ((0 3 2 1) (4 5 6 7));}}); mergePatchPairs ();
''')
    ports=min(ports,nz)
    if ports<nz:
        put('system/topoSetDict',head('topoSetDict')+f'''actions (
{{name seal;type faceSet;action new;source patchToFace;patches (inlet);}}
{{name seal;type faceSet;action subtract;source boxToFace;box (-1 -1 -1) ({dx*.1} {(ports-.1)*dz} 1);}}
{{name sealRight;type faceSet;action new;source patchToFace;patches (outlet);}}
{{name sealRight;type faceSet;action subtract;source boxToFace;box ({length-dx*.1} {height-(ports-.1)*dz} -1) (1 1 1);}}
{{name seal;type faceSet;action add;source faceToFace;set sealRight;}});
''')
        put('system/createPatchDict',head('createPatchDict')+'pointSync false;patches ({name walls;patchInfo {type wall;} constructFrom set;set seal;});\n')
    put('system/controlDict',head('controlDict')+f'''application saltLakeFoam;startFrom startTime;startTime 0;stopAt endTime;
endTime {end};deltaT {dt};writeControl runTime;writeInterval {end};purgeWrite 0;
writeFormat binary;writePrecision 16;writeCompression off;timeFormat general;timePrecision 10;runTimeModifiable false;
''')
    component_div='\n'.join(f'div(massPhi,Y{j}) {component_schemes[component_scheme]};' for j in range(7))
    put('system/fvSchemes',head('fvSchemes')+'''ddtSchemes {default Euler;} gradSchemes {default Gauss linear;}
divSchemes {default Gauss upwind;
'''+component_div+'''
} laplacianSchemes {default Gauss linear orthogonal;}
interpolationSchemes {default linear;} snGradSchemes {default orthogonal;} fluxRequired {default no;p;}
''')
    put('system/fvSolution',head('fvSolution')+'''solvers {
p {solver GAMG;tolerance 1e-20;relTol .0001;smoother GaussSeidel;}
U {solver smoothSolver;smoother symGaussSeidel;tolerance 1e-22;relTol 1e-5;}
"Y.*" {solver PBiCGStab;preconditioner DILU;tolerance 1e-12;relTol 0;}}
SIMPLE {nNonOrthogonalCorrectors 0;}
relaxationFactors {fields {p .7;} equations {U .9;}}
''')
    if solver_profile=='balanced':
        solution=path/'system/fvSolution'
        # These are inner linear solves. The nonlinear and conservation
        # acceptance criteria in reactiveProperties remain unchanged.
        solution.write_text(solution.read_text().replace('tolerance 1e-20;relTol .0001;', 'tolerance 1e-12;relTol .01;')
                            .replace('tolerance 1e-22;relTol 1e-5;', 'tolerance 1e-14;relTol .01;'))
    put('system/decomposeParDict',head('decomposeParDict')+'numberOfSubdomains 8;method scotch;\n')
    return field,put
