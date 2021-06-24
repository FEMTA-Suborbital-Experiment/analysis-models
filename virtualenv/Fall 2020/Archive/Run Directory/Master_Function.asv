clear
clc

% For the simulation time of 200 seconds, this function gathers the necessary outputs  
% to feed to the flight computer. The simulation outputs Pressure and Temperature
% in the Propellant Tank, Pressure Loss across the Orifice in the piping network,  
% and Pressure in the Collection Chamber for the duration of the experiment.

tankPressure;
propellant_tank_pressure=simData.physics.pressure.values(:,4);
propellant_tank_temnperature=simData.physics.temperature.values(:,4);
volume_leaving_propTank=simData.physics.volume.values(:,4);
volumetric_flow_rate=diff(volume_leaving_propTank)/0.03*-1;
PressLossThruOrifice;
pressure_loss_through_orifice=delP_inOrifice;
collection_chamber;
pressure_in_collection_chamber=CCPressure(0:length(CCPressure)-1)-delP_inOrifice;
