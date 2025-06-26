% === Simscape Gas Pipeline Model Setup ===

% Create new Simulink model
gas_model = "gas_pipeline_model";
new_system(gas_model);
open_system(gas_model);

% Add Simscape and Simscape Gas libraries
load_system("simscape");
load_system("simscape_gas");

% Add Gas Properties block (Semiperfect gas: methane)
add_block("simscape_gas/Foundation/Gas/Gas Properties (G)", ...
    [gas_model "/Gas Properties"], "Position", [50 100 200 140]);

% Add Pipe (G) block
add_block("simscape_gas/Sources & Elements/Pipe (G)", ...
    [gas_model "/Pipe"], "Position", [300 100 450 140]);

% Add Controlled Mass Flow Rate Source (G)
add_block("simscape_gas/Sources & Elements/Controlled Mass Flow Rate Source (G)", ...
    [gas_model "/Mass Source"], "Position", [50 200 200 240]);

% Add Reservoir (G) at outlet
add_block("simscape_gas/Sources & Elements/Reservoir (G)", ...
    [gas_model "/Reservoir"], "Position", [600 100 750 140]);

% Connect blocks
add_line(gas_model, "Mass Source/1", "Pipe/1");
add_line(gas_model, "Pipe/2", "Reservoir/1");

% Connect Gas Properties block to Pipe and Mass Source
add_line(gas_model, "Gas Properties/1", "Pipe/3");
add_line(gas_model, "Gas Properties/1", "Mass Source/2");
add_line(gas_model, "Gas Properties/1", "Reservoir/2");

% Set pipe parameters (example: 2000 m steel pipe, diameter 1 m, wall thickness 0.02 m)
set_param([gas_model "/Pipe"], ...
    "Length", "2000", ...
    "Diameter", "1", ...
    "Wall thickness", "0.02", ...
    "Cross-sectional area", "pi*(1/2)^2", ...
    "Internal surface roughness", "1e-5", ...
    "Number of segments", "20");

% Set gas properties: methane
set_param([gas_model "/Gas Properties"], "Gas specification", "Semiperfect");
set_param([gas_model "/Gas Properties"], "Predefined gas", "Methane");

% Set boundary conditions
set_param([gas_model "/Mass Source"], "Mass flow rate", "@(t) 10", ...
    "Inlet pressure", "5000000", ...
    "Inlet temperature", "288.15");

set_param([gas_model "/Reservoir"], "Pressure", "4500000", ...
    "Temperature", "288.15");

% Save and display
save_system(gas_model);
open_system(gas_model);
