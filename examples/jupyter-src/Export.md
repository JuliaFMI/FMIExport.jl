# Create a FMU
at the example of the well-known *bouncing ball*

Tutorial by Tobias Thummerer, Simon Exner | Last edit: June 2026

## License


```julia
#
# Copyright (c) 2021 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#
```

## Introduction

This Julia Package FMIExport.jl enables the export of simulation models from Julia code. The FMI (Functional Mock-up Interface) is a free standard ([fmi-standard.org](https://fmi-standard.org)) that defines a *model container* and an interface to exchange simulation models using a combination of XML files, binaries and C code zipped into a single file. The magic on the Julia side happens by compiling Julia code to a DLL (Windows) or LIB (Linux) by using the *PackageCompiler.jl*, therefore a *package* is defined to be compiled - instead of just a single file. This is very important to keep in mind.

## Getting started 
To allow for editing the example, you can make a copy of the examples folder. In case you just want to lean back and enjoy the magic, you can proceed by just *opening* the bouncing ball example. 

In the creation of the FMU, an entire package is involved. However, the things that need to be defined by you can happen within only a single file. In case of the bouncing ball, this file is [BouncingBall.jl](https://github.com/JuliaFMI/FMIExport.jl/blob/main/examples/FMI2/BouncingBall/src/BouncingBall.jl). If we have a look inside, we find many interesting things, that we investigate in detail now. The code in the following only investigates *snippets*, for full running code see the [examples folder](https://github.com/JuliaFMI/FMIExport.jl/blob/main/examples/).

We start by loading the library.


```julia
using FMIExport
using FMIImport
using FMIExport.FMIBase.FMICore: fmi2True, fmi2False, fmi2Integer
```

Further, we define a start state and initial parameters.


```julia
# a minimum height to reset the ball after event
const EPS = 1e-8

# ball position, velocity (initial)
const DEFAULT_X0 = [1.0, 0.0]

# ball mass, ball radius, ball collision damping, ball minimum velocity, gravity constant 
const DEFAULT_PARAMS = [1.0, 0.1, 0.9, 1e-3, 9.81]
```




    5-element Vector{Float64}:
     1.0
     0.1
     0.9
     0.001
     9.81



### Initialization function
For solving IVPs (initial value problems), we require an *inital value* - for ODEs, this is the initial state. Therefore, we define a function, that computes (or just returns) a start state to start the algorithm for numerical integration to obtain a solution.


```julia
FMU_FCT_INIT = function()
   
    sticking = fmi2False
    counter = fmi2Integer(0)

    s = DEFAULT_X0[1]         # ball position
    v = DEFAULT_X0[2]         # ball velocity
    a = 0.0                   # ball acceleration

    t = 0.0        
    x_c = [s, v]      
    ẋ_c = [v, a]
    x_d = [sticking, counter]
    u = []
    p = DEFAULT_PARAMS

    return (t, x_c, ẋ_c, x_d, u, p)
end
```




    #1 (generic function with 1 method)



### Function for the right-hand side
After we defined start values for the integration, we need to define what is actually integrated. For ODEs, this is refered to as the *right-hand side*, a function, that defines derivatives that are numerically integrated by the ODE solver.


```julia
FMU_FCT_EVALUATE = function(t, x_c, ẋ_c, x_d, u, p, eventMode)
    m, r, d, v_min, g = p
    s, v = x_c
    sticking, counter = x_d
    _, a = ẋ_c

    if sticking == fmi2True
        a = 0.0
    elseif sticking == fmi2False

        if eventMode
            h = s-r
            if h <= 0 && v < 0
                s = r + EPS # so that indicator is not triggered again
                v = -v*d 
                counter = fmi2Integer(counter+1)

                # stop bouncing to prevent high frequency bouncing (and maybe tunneling the floor)
                if abs(v) < v_min
                    sticking = fmi2True
                    v = 0.0
                end
            end
        end

        a = (m * -g) / m     # the system's physical equation (a little longer than necessary)
    else
        @error "Unknown value for `sticking` == $(sticking)."
        return (x_c, ẋ_c, x_d, p)
    end

    # Todo: Remove these allocations. Make it inplace.
    x_c = [s, v]
    ẋ_c = [v, a]
    x_d = [sticking, counter]
    p = [m, r, d, v_min, g]

    return (x_c, ẋ_c, x_d, p)
end
```




    #3 (generic function with 1 method)



### Event indicator function
We defined what happens if we collide with the ground, but we need a way to tell the solver *when exactly* this is the case - so a mathematical definition of the collision. This happens within the *event indicator*: Whenever we have a real zero crossing within the event indicator, the solver (or the event finding routine within) searches for the exact event location and pauses the integration process here. After the event is handled (and the ball velocity changed), the numerical integration resumes with the new state.


```julia
FMU_FCT_EVENT = function(t, x_c, ẋ_c, x_d, u, p)
    m, r, d, v_min, g = p
    s, v = x_c
    _, a = ẋ_c
    sticking, counter = x_d

    # helpers
    z1 = 0.0 # first event indicator
    h = s-r # ball height

    if sticking == fmi2True
        z1 = 1.0            # event 1: ball stay-on-ground
    else
        if h > -EPS && h <= 0 && v > 0
            z1 = -EPS
        else 
            z1 = h
        end
    end

    z = [z1]

    return z
end
```




    #5 (generic function with 1 method)



### Output function

Finally, we can define *outputs* that are exposed to the user. In FMI, in theory, we could investigate all quantities we are interested in (if not explicitly hidden by the FMU), however it makes sense to define or calculate application specific outputs. In this case, we just return the ball position and velocity.


```julia
FMU_FCT_OUTPUT = function(t, x_c, ẋ_c, x_d, u, p)
    m, r, d, v_min, g = p
    s, v = x_c
    _, a = ẋ_c
    sticking, counter = x_d

    y = [s, v]

    return y
end
```




    #7 (generic function with 1 method)



Because we want to export ME as well as CS, we need to specify an ODE solver for the CS FMU. This is not necessary, if you are only using ME.


```julia
import OrdinaryDiffEqTsit5: Tsit5
FMU_FCT_SOLVER = function ()
    return Tsit5()
end
```




    #9 (generic function with 1 method)



### FMU constructor
After defining the mathematical behavior for our FMU, we need to define a constructor - so a function that is called if the FMU is loaded. Within, we define all the variables we want to be part of the *model description*, technically an XML that exposes the model structure and variables. This XML is parsed by importing tools to simulate FMUs appropriately and to provide helpful information.


```julia
FMIBUILD_CONSTRUCTOR = function(resPath="")
    fmu = createFMU2Simple("BouncingBall";
                        initializationFct=FMU_FCT_INIT,
                        evaluationFct=FMU_FCT_EVALUATE,
                        outputFct=FMU_FCT_OUTPUT,
                        eventFct=FMU_FCT_EVENT,
                        solverFct=FMU_FCT_SOLVER)

    # modes 
    addModelExchange(fmu)
    addCoSimulation(fmu)

    # states [2]
    addStateAndDerivative(fmu, "ball.s"; stateStart=DEFAULT_X0[1], stateDescr="Absolute position of ball center of mass", derivativeDescr="Absolute velocity of ball center of mass")
    addStateAndDerivative(fmu, "ball.v"; stateStart=DEFAULT_X0[2], stateDescr="Absolute velocity of ball center of mass", derivativeDescr="Absolute acceleration of ball center of mass")

    # discrete state [2]
    addIntegerDiscreteState(fmu, "sticking"; description="Indicator (boolean) if the mass is sticking on the ground, as soon as abs(v) < v_min")
    addIntegerDiscreteState(fmu, "counter"; description="Number of collision with the floor.")

    # outputs [2]
    addRealOutput(fmu, "ball.s_out"; description="Absolute position of ball center of mass")
    addRealOutput(fmu, "ball.v_out"; description="Absolute velocity of ball center of mass")

    # parameters [5]
    addRealParameter(fmu, "m";     start=DEFAULT_PARAMS[1], description="Mass of ball")
    addRealParameter(fmu, "r";     start=DEFAULT_PARAMS[2], description="Radius of ball")
    addRealParameter(fmu, "d";     start=DEFAULT_PARAMS[3], description="Collision damping constant (velocity fraction after hitting the ground)")
    addRealParameter(fmu, "v_min"; start=DEFAULT_PARAMS[4], description="Minimal ball velocity to enter on-ground-state")
    addRealParameter(fmu, "g";     start=DEFAULT_PARAMS[5], description="Gravity constant")

    addEventIndicator(fmu)

    return fmu
end
```




    #11 (generic function with 2 methods)



## Building the FMU
The compilation process for the FMU is triggered by the following lines. Keep in mind, that the entire package is compiled to an FMU. If we trigger the building process within this package, we would create a FMU that compiles FMU... this sounds cool at first glance, but is not what we want in most cases. Therefore, the build command itself `saveFMU` is within a formatted block, starting with `### FMIBUILD_NO_EXPORT_BEGIN ###` and ending with `### FMIBUILD_NO_EXPORT_END ###`. All text between these delimiters is automatically removed *before* building of the package happens. This way, the actual building package `FMIBuild` and the build command are excluded from the final FMU. Alternatively, one could trigger the build process from outside the package, e.g. by working in the REPL directly. The command `saveFMU` is commented out, to not start a building process within this notebook - this must be included in your application of course.


```julia
### FMIBUILD_NO_EXPORT_BEGIN ###
# The line above is a start-marker for excluded code for the FMU compilation process!

tmpDir = mktempdir(; prefix="fmibuildjl_test_", cleanup=false) 
@info "Saving example files at: $(tmpDir)"
fmu_save_path = joinpath(tmpDir, "BouncingBall.fmu")  

fmu = FMIBUILD_CONSTRUCTOR()
using FMIBuild: saveFMU                    # <= this must be excluded during export, because FMIBuild cannot execute itself (but it is able to build)
# saveFMU(fmu, fmu_save_path; debug=true, compress=false)    # <= this must be excluded during export, because saveFMU would start an infinite build loop with itself (debug=true allows debug messages, but is slow during execution!)

# The following line is a end-marker for excluded code for the FMU compilation process!
### FMIBUILD_NO_EXPORT_END ###
```

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39mSaving example files at: C:\Users\RUNNER~1\AppData\Local\Temp\fmibuildjl_test_qB9HOT
    

## Simulate the FMU 
FMUs exported and compiled with Julia **can not** be executed within Julia. This sounds wired, but comes from the fact that only a single Julia instance (sysimage) is allowed to run for each Julia process. However, you can use FMUs from Julia within any other tool that supports FMI, of course.

Interestingly, there is actually no need to *compile* a FMU if you want to use it within Julia anyway - you can just use the current copy within the memory with *FMI.jl*, see the following lines.


```julia
using FMIImport, Plots
fmu.executionConfig.loggingOn = true
```




    true



Now we can simulate the FMU in ME ...


```julia
solution = simulateME(fmu, (0.0, 3.0); 
    recordValues=["ball.s", "counter"], 
    saveat=0.0:0.01:3.0, 
    solver=Tsit5(), # instead of specifying a solver by hand, you can import DifferentialEquations.jl and use auto-picking a solver
)
plot(solution; states=false) # don't plot states
```

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mfmi2Instantiate!(...): This component was already registered. This may be because you created the FMU by yourself with FMIExport.jl.
    [33m[1m└ [22m[39m[90m@ FMIBase C:\Users\runneradmin\.julia\packages\FMIBase\OYsGD\src\printing.jl:30[39m
    

    [34mSimulating ME-FMU ...   0%|█                            |  ETA: N/A[39m

    [34mSimulating ME-FMU ... 100%|█████████████████████████████| Time: 0:00:05[39m
    

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mUnknown instance at Ptr{Nothing} @0x00000270cc0d2ad0.
    [33m[1m└ [22m[39m[90m@ FMIExport D:\a\FMIExport.jl\FMIExport.jl\src\FMI2\simple.jl:51[39m
    




<img src="data:image/png;base64,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" />



... as well as in CS ...


```julia
solution = simulateCS(fmu, (0.0, 3.0); recordValues=["ball.s", "counter"], saveat=0.0:0.01:3.0)
plot(solution; states=false)
```

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mfmi2Instantiate!(...): This component was already registered. This may be because you created the FMU by yourself with FMIExport.jl.
    [33m[1m└ [22m[39m[90m@ FMIBase C:\Users\runneradmin\.julia\packages\FMIBase\OYsGD\src\printing.jl:30[39m
    

    [34mSim. CS-FMU ...   0%|█                                  |  ETA: 0:07:00[39m

    [34mSim. CS-FMU ... 100%|███████████████████████████████████| Time: 0:00:01[39m
    

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mUnknown instance at Ptr{Nothing} @0x00000270cadd4490.
    [33m[1m└ [22m[39m[90m@ FMIExport D:\a\FMIExport.jl\FMIExport.jl\src\FMI2\simple.jl:51[39m
    




<img src="data:image/png;base64,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" />



## Getting started with custom project
The best way to start with your own project, is to make a *copy* of one of the example *folders*, for example the bouncing ball and edit it to fit you requirements. Because the entire package is compiled, you can control the size of the binary by adding or removing packages to the *Project.toml*. You can use *Pkg.jl* of course to manage dependencies. Then you can start overwriting the functions, that define the right-hand side of you FMU, events and event indicators and finally the output variables.
