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

    [36m[1m[ [22m[39m[36m[1mInfo: [22m[39mSaving example files at: C:\Users\RUNNER~1\AppData\Local\Temp\fmibuildjl_test_2WYMPq
    

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

    [34mSimulating ME-FMU ... 100%|█████████████████████████████| Time: 0:00:06[39m
    

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mUnknown instance at Ptr{Nothing} @0x00000214bfa1a890.
    [33m[1m└ [22m[39m[90m@ FMIExport D:\a\FMIExport.jl\FMIExport.jl\src\FMI2\simple.jl:56[39m
    




<img src="data:image/png;base64,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" />



... as well as in CS ...


```julia
solution = simulateCS(fmu, (0.0, 3.0); recordValues=["ball.s", "counter"], saveat=0.0:0.01:3.0)
plot(solution; states=false)
```

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mfmi2Instantiate!(...): This component was already registered. This may be because you created the FMU by yourself with FMIExport.jl.
    [33m[1m└ [22m[39m[90m@ FMIBase C:\Users\runneradmin\.julia\packages\FMIBase\OYsGD\src\printing.jl:30[39m
    [34mSim. CS-FMU ...   0%|█                                  |  ETA: 0:09:14[39m

    [34mSim. CS-FMU ... 100%|███████████████████████████████████| Time: 0:00:01[39m
    

    [33m[1m┌ [22m[39m[33m[1mWarning: [22m[39mUnknown instance at Ptr{Nothing} @0x00000214bf84d450.
    [33m[1m└ [22m[39m[90m@ FMIExport D:\a\FMIExport.jl\FMIExport.jl\src\FMI2\simple.jl:56[39m
    




<img src="data:image/png;base64,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" />



## Getting started with custom project
The best way to start with your own project, is to make a *copy* of one of the example *folders*, for example the bouncing ball and edit it to fit you requirements. Because the entire package is compiled, you can control the size of the binary by adding or removing packages to the *Project.toml*. You can use *Pkg.jl* of course to manage dependencies. Then you can start overwriting the functions, that define the right-hand side of you FMU, events and event indicators and finally the output variables.
