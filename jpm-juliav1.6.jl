# module JuliaPhotochem
using Measures
using DifferentialEquations
using Catalyst
using Random
using Plots



export ChemSys, jpm
press=1.0
tempk=298.0
# Set the value of hv
hv = 1

# seed the random number generator
Random.seed!(42)


arrhenius(A,EdivR,tempk) = A * exp(-EdivR/tempk)
molcm3frompressure(tempk,press) = 6.022e23 * press / (.0821 * tempk) / 1e3 # 6.022e23 is Avogadro's number, .0821 is R in L atm K-1 mol-1, and there are 1e3 cm3 in a L
function termolecular(k₀298,n, k∞298, m, tempk,press)
    # calculate rate constant in the low pressure limit
    k₀ = k₀298 * (tempk/298.0)^(-n) # cm6 molec-2 s-1
    # calculate rate constant in the high pressure limit
    k∞ = k∞298 * (tempk/298.0)^(-m) # cm3 molec-1 s-1
    # calculate M, the total gas Concentration
    M = molcm3frompressure(tempk,press) # molec/cm3
    # calculate the effective second-order rate constant
    kf = ( k∞ * k₀ * M ) / (k∞ + k₀ * M ) * 0.6^(1.0 / ( 1.0 + (log10(k₀*M/k∞))^2.0 ))
    return kf 
end


function newmodel(press=1.0::Float64, tempk=298.0::Float64)
    # function definition to convert from cm3 molec-1 s-1 to (ppm min)-1
	# 7.34e15 is Avogadro's number divided by R in L atm K-1 mol-1 divided by 1000 to convert to cm3 from L and 1e6 to convert to (ppm)-1
	# These units are valid for bimolecular reactions
    # rkppmmin(rk,press,tempk) = rk * 60.0 *  7.34e+15 * press / tempk 

    # make a function arrhenius(A,(E/R),tempk) to calculate the rate constant

    # Reaction rate constants

    # Old reaction 2 was O + O₂ → O₃
    # removed July 2024, JPL 19-5 page 434 Table 2-1
    # original was rkppmmin(6.0e-34 * ( tempk/300. )^(-2.3) * 7.34e+21 * press / tempk, press, tempk)
    # rk2(tempk,press) = molcm3frompressure(tempk,press)*6.1e-34*(tempk/298)^(-2.4)

    # Reaction 1 is NO₂ → NO + O₃
    # photolysis reaction, rate constant from original mechanism below
    # rk1(hv) = 0.5 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk1(hv) = 9.258186912960969E-003 * hv

    # Reaction 2 is O₃ + NO → NO₂ + O₂
    # updated July 2024, JPL 19-5 page 84 Table 1C 
    # original was rkppmmin(2.e-12 * exp(-1400/tempk), press, tempk)
    rk2(tempk) = arrhenius(3.0e-12,1500,tempk) 

    # Reaction 3 is HCHO + O₂ → HO₂ + CO
    # photolysis reaction, rate constant taken from original mechanism
    # rk3(hv) = 0.015 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk3(hv) = 3.442262702424394E-005 * hv

    # Reaction 4 is HCHO → H₂ + CO
    # photolysis reaction, rate constant from original mechanism
    # rk4(hv) = 0.022 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk4(hv) = 5.132233465281091E-005 * hv


    # Reaction 5 is HCHO + OH + O₂ → HO₂ + CO + H₂O
    # updated July 2024, JPL 19-5 page 106 Table 1D
    # OH + H₂CO (not HCHO or CH₂O, lol) → H₂O + HCO 
    # (formyl radical in oxygen makes HO₂ and CO, GEOS-Chem also parameterizes this as instantaneous)
    # original rate was rkppmmin(1.2e-14 * tempk^1 * exp(+287/tempk), press, tempk)
    rk5(tempk) = arrhenius(5.5e-12,-125.0,tempk)

    # Reaction 6 is HO₂ + NO → OH + NO₂
    # updated July 2024, JPL 19-5 page 83 Table 1C
    # original was rk7(press, tempk) = rkppmmin(3.7e-12 * exp(+250/tempk), press, tempk)
    rk6(tempk) = arrhenius(3.44e-12,-260,tempk)

    # Reaction 7 is OH + NO₂ → HNO₃
    # updated July 2024, JPL 19-5 page 434 Table 2-1
    # the original rate is below
    # function rk8(press, tempk)
    #     rk0 = 2.0e-30 * ( tempk/300 )^(-3.0)
    #     rki = 2.5e-11 * ( tempk/300 )^0
    #     rm = 7.34e+21 * press / tempk
    #     rb = 0.6^(1.0 / ( 1.0 + (log10(rk0*rm/rki))^2.0 ))
    #     rk8 = rk0 * rm * rb / (1.0 + rk0*rm/rki)
    #     rk8 = rkppmmin(rk8,press,tempk)
    # end
    rk7(tempk,press) = termolecular(1.8e-30,3.0,2.8e-11,0,tempk,press)

    # Reaction 8 is H₂O₂ → 2HO
    # this is a photolysis reaction, rate constant from original mechanism below
    # rk8(hv) = 0.0003 * hv / 60 # converting to seconds
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk8(hv) = 7.373314014460711E-006 * hv

    # Reaction 9 is H₂O₂ + OH → H₂O + HO₂
    # updated July 2024, JPL 19-5 page 73
    # "the recommendation is a temperature independent value of 1.8e-12 cm3 molec-1 s-1"
    # original was rkppmmin(3.3e-12 * exp(-200.0/tempk),press,tempk)
    rk9() = 1.8e-12

    # Reaction 10 is ALD2 + OH → MCO₃ + H₂O
    # updated July 2024, JPL 19-5 page 107 Table 1D
    # follows an Arrhenius equation with A = 4.63e-12, E/R = -350
    # note that GCARR_ac has the negative in it already so 
    # this is equivalent to GCARR_ac(4.63d-12, 350.0d0);
    rk10(tempk) = arrhenius(4.63e-12,-350,tempk)

    # Reaction 11 is MGLY + hv → MCO₃ + CO + HO₂ 
    # this is a photolysis reaction with rate constant sampled from the GEOS-CF
    # chemical state at the surface of Los Angeles, using the KPP Standalone interface
    rk11(hv) = 1.863860764162095E-004 * hv 


    # Reaction 12 is MCO₃ + NO₂ {+M} → PAN
    # updated July 2024, JPL 19-5 page 435 Table 2-1
    rk12(tempk,press) = termolecular(7.3e-29,4.1,9.5e-12,1.6,tempk,press) 

    # Reaction 13 is thermal decomposition of PAN → MCO₃ + NO₂ 
    # updated July 2024, JPL 19-5 page 521 Table 3-1
    # using the equilibrium constant divided by the same termolecular function as reaction 12
    # to get the backward rate constant
    rk13(tempk,press) = termolecular(7.3e-29,4.1,9.5e-12,1.6,tempk,press) / ( 9.0e-29*exp(14000/tempk) )




    # Evaluate the rate constant functions to get their values
    p = [rk1(hv), rk2(tempk), rk3(hv), rk4(hv), rk5(tempk),
        rk6(tempk), rk7(press,tempk), rk8(hv), rk9(), rk10(tempk), 
        rk11(hv), rk12(tempk,press), rk13(tempk,press)]


    rn = @reaction_network begin
        # Note, removing O₂ April 2025 for CRNNs. Treat instead as infinite source/sink
        rk1, NO₂ → NO + O₃
        rk2, O₃ + NO → NO₂
        rk3, HCHO → 2HO₂ + CO  
        rk4, HCHO → H₂ + CO
        rk5, HCHO + OH → HO₂ + CO + H₂O 
        rk6, HO₂ + NO → OH + NO₂
        rk7, OH + NO₂ → HNO₃
        rk8, H₂O₂ → 2OH
        rk9, H₂O₂ + OH → H₂O + HO₂
        rk10, ALD2 + OH → MCO₃ + H₂O
        rk11, MGLY → MCO₃ + CO + HO₂
        rk12, MCO₃ + NO₂ → PAN
        rk13, PAN → MCO₃ + NO₂
        # All species go to null for testing
        # rk2*1e12, NO₂ → ∅;
        # rk5*1e12, NO → ∅;
        # rk6*1e12, O₃ → ∅;
        # rk7*1e12, OH → ∅;
        # rk9*1e12, HO₂ → ∅;
        # rk2*1e12, HCHO → ∅;
        # rk6*1e12, H₂O₂ → ∅;
        # rk7*1e12, HNO₃ → ∅;
        # rk9*1e12, CO → ∅;
        # rk12*1e12, H₂ → ∅;
        # rk13*1e12, H₂O → ∅;
        # rk12*1e12, MCO₃ → ∅;
        # rk13*1e12, PAN → ∅;
        # rk13*1e12, ALD2 → ∅;
        # rk13*1e12, MGLY → ∅;
    end 

    nspecs = size(species(rn),1)

    specmap = Dict()
    for (k, v) in speciesmap(rn)
        specmap[replace(string(k), r"\(t\)$"=>"")] = v # remove '(t)' from names
    end
    natoms = 4 # C N H
    atoms = zeros(Float64, nspecs, natoms)
    atoms[specmap["O₃"],   :]   = [0 0 0 3]
    atoms[specmap["NO"],   :]   = [0 1 0 1]
    atoms[specmap["NO₂"],  :]   = [0 1 0 2]
    atoms[specmap["HCHO"], :]   = [1 0 2 1]
    atoms[specmap["HO₂"],  :]   = [0 0 1 2]
    atoms[specmap["H₂O₂"], :]   = [0 0 2 2]
    atoms[specmap["OH"],   :]   = [0 0 1 1]
    atoms[specmap["HNO₃"], :]   = [0 1 1 3]
    atoms[specmap["CO"],   :]   = [1 0 0 1]
    atoms[specmap["H₂"],   :]   = [0 0 2 0]
    atoms[specmap["H₂O"],  :]   = [0 0 2 1]
    atoms[specmap["ALD2"], :]   = [2 0 4 1]
    atoms[specmap["MGLY"], :]   = [3 0 4 2]
    atoms[specmap["MCO₃"], :]   = [2 0 3 3]
    atoms[specmap["PAN"],  :]   = [2 1 3 5]

    return rn, atoms, specmap, p
end

struct ChemSys
    rn
    atoms::Array{Float64, 2}
    specmap::Dict
    p::Array{Float64, 1}
end


function ppb_to_molec_cm3(ppb::Float64, tempk::Float64, press::Float64)
    molec_cm3 = ppb / 1e9 * press / .0821 / tempk * 6.022e23 / 1e3
    return molec_cm3
end
function molec_cm3_to_ppb(molec_cm3::Float64, tempk::Float64, press::Float64)
    ppb = molec_cm3 / 6.022e23 * 1e3 * .0821 * tempk / press * 1e9
    return ppb
end


ChemSys(press=1.0::Real, tempk=298.0::Real) = ChemSys(newmodel(press, tempk)...)
jpm = ChemSys()

function plot_trajectories(jpm, tempk, press, vars_to_plot=["O₃", "NO", "NO₂", "HCHO", "HO₂", "H₂O₂", "OH", "HNO₃", "CO", "H₂", "ALD2", "MGLY", "MCO₃", "PAN", "H₂O"])
    # Initialize c0 with zeros
    c0 = zeros(Float64, length(jpm.specmap))
    # Set specific initial conditions

    # Initialize concentrations from uniform distribution
    c0[jpm.specmap["O₃"], :]   .= 50*rand()         
    c0[jpm.specmap["NO"], :]   .= 20*rand()   
    c0[jpm.specmap["NO₂"], :]  .= 20*rand()
    c0[jpm.specmap["HCHO"], :] .= 50*rand()  
    c0[jpm.specmap["HO₂"], :]  .= 0.01*rand()
    c0[jpm.specmap["H₂O₂"], :] .= 10*rand()
    c0[jpm.specmap["OH"], :]   .= 0.001*rand() 
    c0[jpm.specmap["MCO₃"], :] .= 0.01*rand()
    c0[jpm.specmap["PAN"], :]  .= 25*rand()
    c0[jpm.specmap["ALD2"], :] .= 10*rand()
    c0[jpm.specmap["MGLY"], :] .= 50*rand()
    # Convert concentrations to molecules/cm^3
    c0 = ppb_to_molec_cm3.(c0, tempk, press)

    # Define the time span for the simulation
    tspan = (0, 3600)
    # Create the ODE problem
    prob = ODEProblem(jpm.rn, c0, tspan, jpm.p,combinatoric_ratelaws = false)

    # Solve the ODE problem
    sol = solve(prob, Rosenbrock23(), saveat=1.0, progress=true)

    # Get indices of variables to plot
    indices = [jpm.specmap[var] for var in vars_to_plot]
    # Extract the variables at every time step
    selected_vars = [state[indices] for state in sol.u]
    # Transpose the matrix to get each variable as a separate vector
    selected_vars_transposed = transpose(hcat(selected_vars...))
    # Convert to ppb
    selected_vars_transposed = molec_cm3_to_ppb.(selected_vars_transposed, tempk, press)

    # Plot each variable as its own subplot
    num_vars = length(vars_to_plot)
    fig_height = 200 * num_vars  # Adjust height per variable
    fig_width = 800  # Standard width
    p = plot(layout = (num_vars, 1), size = (fig_width, fig_height))
    # Add each variable to the plot
    for (i, var) in enumerate(vars_to_plot)
        plot!(p[i], sol.t, selected_vars_transposed[:,i], legend=false, left_margin=6Measures.mm)
        ylims!(p[i], (0, maximum(selected_vars_transposed[:,i])))
        xlabel!(p[i], "Time [seconds]")
        ylabel!(p[i], var * " [ppb]")
        title!(p[i], var)
    end
    display(p)
end

function plot_selected_vars_transposed(selected_vars_transposed, vars_to_plot, sol_t, tempk, press)
    # Plot each variable as its own subplot
    num_vars = length(vars_to_plot)
    fig_height = 200 * num_vars  # Adjust height per variable
    fig_width = 800  # Standard width
    p = plot(layout = (num_vars, 1), size = (fig_width, fig_height))
    # Add each variable to the plot
    for (i, var) in enumerate(vars_to_plot)
        plot!(p[i], sol_t, selected_vars_transposed[:,i], legend=false, left_margin=6Measures.mm)
        # ylims!(p[i], (0, maximum(selected_vars_transposed[:,i])))
        xlabel!(p[i], "Time [seconds]")
        ylabel!(p[i], var * " [ppb]")
        title!(p[i], var)
    end
    display(p)
end

# plot_trajectories(jpm, tempk, press)

# mass balance checker:
# does jpm.atoms * the tendencies equal zero?
# if so, mass is conserved
# tendencies are the difference in concentration at each time step
# M = jpm.atoms[indices,:]
# tendencies = diff(selected_vars_transposed, dims=1)
# mass_balance = tendencies*M



# function to run N simulations of 1 hour, with random initial conditions, sampling every 5 minutes
function run_simulations(num_tests::Int, duration::Float64, sampling_interval::Float64)
    # num_tests = 1
    # duration = 3600.0
    # sampling_interval = 300.0
    # print number of experiments
    
    indices = [jpm.specmap[var] for var in ["O₃", "NO", "NO₂", "HCHO", "HO₂", "H₂O₂", "OH", "HNO₃", "CO", "H₂", "ALD2", "MGLY", "MCO₃", "PAN", "H₂O"]]
    label = ["Experiment number", "Time (min)", "O3 [ppb]", "NO [ppb]", "NO2 [ppb]", "HCHO [ppb]", "HO2 [ppb]", "H2O2 [ppb]", "OH [ppb]", "HNO3 [ppb]", "CO [ppb]", "H2 [ppb]", "ALD2 [ppb]", "MGLY [ppb]", "MCO3 [ppb]", "PAN [ppb]", "H2O [ppb]"]
    num_time_points = Int(duration / sampling_interval)+1
    experiments_array = Array{Any}(undef, (num_tests * num_time_points + 1), length(label))
    experiments_array[1,:] = label

    row_index = 2
    for i in 1:num_tests
        println("Experiment number: ", i)
        # Initialize c0 with zeros
        c0 = zeros(Float64, length(jpm.specmap))
        # Set specific initial conditions
        # Initialize concentrations from uniform distribution
        c0[jpm.specmap["O₃"], :]   .= 50*rand()         
        c0[jpm.specmap["NO"], :]   .= 20*rand()   
        c0[jpm.specmap["NO₂"], :]  .= 20*rand()
        c0[jpm.specmap["HCHO"], :] .= 50*rand()  
        c0[jpm.specmap["HO₂"], :]  .= 0.01*rand()
        c0[jpm.specmap["H₂O₂"], :] .= 10*rand()
        c0[jpm.specmap["OH"], :]   .= 0.001*rand() 
        c0[jpm.specmap["MCO₃"], :] .= 0.01*rand()
        c0[jpm.specmap["PAN"], :]  .= 25*rand()
        c0[jpm.specmap["ALD2"], :] .= 10*rand()
        c0[jpm.specmap["MGLY"], :] .= 50*rand()
        # Convert concentrations to molecules/cm^3
        c0 = ppb_to_molec_cm3.(c0, tempk, press)
        # Define the time span for the simulation
        tspan = (0, duration)
        # Create the ODE problem
        # prob = ODEProblem(jpm.rn, c0, tspan, jpm.p)
        prob = ODEProblem(jpm.rn, c0, tspan, jpm.p, combinatoric_ratelaws = false)
        # Solve the ODE problem
        global sol = solve(prob, Rosenbrock23(), saveat=sampling_interval, progress=false)
        # Extract the variables at every time step
        selected_vars = [state[indices] for state in sol.u]
        # create array
        global selected_vars_transposed = transpose(hcat(selected_vars...))
        # Convert to ppb
        selected_vars_transposed = molec_cm3_to_ppb.(selected_vars_transposed, tempk, press)
        for (t, time_point) in enumerate(0:sampling_interval:duration)
            # print the time point
            # println("Time point: ", time_point)
            experiments_array[row_index, :] = [i, time_point, selected_vars_transposed[t, :]...]
            row_index += 1
        end
    end
    return experiments_array
end

# number_of_experiments = Integer(1e6)
number_of_experiments = Integer(1)
# number_of_experiments = 3
experiments_array = run_simulations(number_of_experiments, 3600.0, 300.0)

# mass balance checker
M = zeros(Int, 15, 3)
M[1,  :] = [0, 0, 0]  # O₃
M[2,  :] = [0, 1, 0]  # NO
M[3,  :] = [0, 1, 0]  # NO₂
M[4,  :] = [1, 0, 2]  # HCHO
M[5,  :] = [0, 0, 1]  # HO₂
M[6,  :] = [0, 0, 2]  # H₂O₂
M[7,  :] = [0, 0, 1]  # OH
M[8,  :] = [0, 1, 1]  # HNO₃
M[9,  :] = [1, 0, 0]  # CO
M[10, :] = [0, 0, 2]  # H₂
M[11, :] = [2, 0, 4]  # ALD2
M[12, :] = [3, 0, 4]  # MGLY
M[13, :] = [2, 0, 3]  # MCO₃
M[14, :] = [2, 1, 3]  # PAN
M[15, :] = [0, 0, 2]  # H₂O

# get the tendencies for all experiments
tendencies_diff = diff(experiments_array[2:end, 3:end], dims=1)
# skip every 12 rows to get the tendencies for each time point
tendencies = tendencies_diff[1:13:end, :]

maximum(tendencies*M, dims=1)

