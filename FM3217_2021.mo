<<<<<<< Updated upstream
within ;
package FM3217_2021 "Course files of the 2021 tutorial"

  package Tutorial1 "Simple Pendulum"

    model SimplePendulum "First simple version"
      constant Real g(unit="m/s2")=9.81 "Gravitationl constant";
      parameter Real L(unit="m")=1 "Length of the pendulum";
      Real Theta(unit="rad",start=0.1) "Angle of the pendulum";
      Real ThetaDot "Helping variable for 2nd derivative";
    equation
       ThetaDot = der(Theta);
       der(ThetaDot) = - g/L * sin(Theta);
      annotation (experiment(StopTime=10));
    end SimplePendulum;

    model SimplePendulumSIunits "Model using SI units"
      constant Modelica.SIunits.Acceleration g=9.81 "Gravitationl constant";
      parameter Modelica.SIunits.Length L=1 "Length of the pendulum";
      Modelica.SIunits.Angle Theta(start=0.1) "Angle of the pendulum";
      Real ThetaDot  "Helping variable for 2nd derivative";
    equation
       ThetaDot = der(Theta);
       der(ThetaDot) = - g/L * sin(Theta);
      annotation (experiment(StopTime=10));
    end SimplePendulumSIunits;

    model SimplePendulumUsingImports "Model using import"
      import SI = Modelica.SIunits;

      constant SI.Acceleration g=9.81 "Gravitationl constant";
      parameter SI.Length L=1 "Length of the pendulum";
      SI.Angle Theta(start=0.1) "Angle of the pendulum";
      Real ThetaDot  "Helping variable for 2nd derivative";
    equation
       ThetaDot = der(Theta);
       der(ThetaDot) = - g/L * sin(Theta);
      annotation (experiment(StopTime=10));
    end SimplePendulumUsingImports;
  end Tutorial1;

  package Tutorial2 "Motor and Motor Drive"

    model Motor

      parameter Modelica.SIunits.Resistance Ra = 0.5 "Resistance of the armature";


      Modelica.Electrical.Analog.Basic.Resistor resistor(R=Ra)  annotation (Placement(transformation(extent={{-56,10},{-36,30}})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-80,-62},{-60,-42}})));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L=La)   annotation (Placement(transformation(extent={{-20,10},{0,30}})));
      Modelica.Electrical.Analog.Basic.EMF emf annotation (Placement(transformation(extent={{10,-10},{30,10}})));
      Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=-90,
            origin={-70,0})));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J=Ja)    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
      Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b1 "Flange of right shaft"
                                                annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      parameter Modelica.SIunits.Inductance La=0.05 "Inductance of the armature";
      parameter Modelica.SIunits.Inertia Ja=0.05 "Inertia of the motor";
    equation
      connect(signalVoltage.p, resistor.p) annotation (Line(points={{-70,10},{-70,20},{-56,20}}, color={0,0,255}));
      connect(resistor.n, inductor.p) annotation (Line(points={{-36,20},{-20,20}}, color={0,0,255}));
      connect(inductor.n, emf.p) annotation (Line(points={{0,20},{20,20},{20,10}}, color={0,0,255}));
      connect(emf.n, signalVoltage.n) annotation (Line(points={{20,-10},{20,-28},{-70,-28},{-70,-10}}, color={0,0,255}));
      connect(emf.n, ground.p) annotation (Line(points={{20,-10},{20,-28},{-70,-28},{-70,-42}}, color={0,0,255}));
      connect(emf.flange, inertia.flange_a) annotation (Line(points={{30,0},{40,0}}, color={0,0,0}));
      connect(signalVoltage.v, u) annotation (Line(points={{-82,0},{-120,0}}, color={0,0,127}));
      connect(inertia.flange_b, flange_b1) annotation (Line(points={{60,0},{100,0}}, color={0,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Bitmap(extent={{-98,-76},{102,128}}, fileName="modelica://FM3217_2021/Resources/Images/dc-motor.jpg")}),
                                                                     Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Motor;

    model MotorDrive
      Motor motor(Ja=2) annotation (Placement(transformation(extent={{-2,-10},{18,10}})));
      Modelica.Blocks.Sources.Step step(startTime=0.1) annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      Modelica.Blocks.Continuous.PID PID(Ti=0.1, Td=0) annotation (Placement(transformation(extent={{-38,-10},{-18,10}})));
      Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio=100) annotation (Placement(transformation(extent={{34,-10},{54,10}})));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J=5) annotation (Placement(transformation(extent={{70,-10},{90,10}})));
      Modelica.Mechanics.Rotational.Sensors.AngleSensor angleSensor annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={90,-30})));
    equation
      connect(step.y, feedback.u1) annotation (Line(points={{-79,0},{-68,0}}, color={0,0,127}));
      connect(feedback.y, PID.u) annotation (Line(points={{-51,0},{-40,0}}, color={0,0,127}));
      connect(idealGear.flange_a, motor.flange_b1) annotation (Line(points={{34,0},{18,0}}, color={0,0,0}));
      connect(angleSensor.flange, inertia.flange_b) annotation (Line(points={{90,-20},{90,0}}, color={0,0,0}));
      connect(angleSensor.phi, feedback.u2) annotation (Line(points={{90,-41},{90,-50},{-60,-50},{-60,-8}}, color={0,0,127}));
      connect(PID.y, motor.u) annotation (Line(points={{-17,0},{-4,0}}, color={0,0,127}));
      connect(idealGear.flange_b, inertia.flange_a) annotation (Line(points={{54,0},{70,0}}, color={0,0,0}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{120,100}})),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
        experiment(StopTime=100));
    end MotorDrive;
  end Tutorial2;

  package Tutorial3 "Simple machine with loads"
    package Components
      model Machine

         parameter Modelica.SIunits.Inductance La=0.05 "Inductance of the armature" annotation (Dialog(group="Electrical"));
        parameter Modelica.SIunits.Inertia Ja=0.05 "Inertia of the motor" annotation (Dialog(group="Mechanical"));

         parameter Modelica.SIunits.Resistance Ra=0.5 "Resistance of the armature" annotation (Dialog(group="Electrical"));

        Modelica.Electrical.Analog.Basic.Resistor resistor(R=Ra)  annotation (Placement(transformation(extent={{-56,10},{-36,30}})));
        Modelica.Electrical.Analog.Basic.Inductor inductor(L=La)   annotation (Placement(transformation(extent={{-20,10},{0,30}})));
        Modelica.Electrical.Analog.Basic.EMF emf annotation (Placement(transformation(extent={{10,-10},{30,10}})));
        Modelica.Mechanics.Rotational.Components.Inertia inertia(J=Ja)    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange "Flange of right shaft" annotation (Placement(transformation(extent={{90,-10},{110,10}})));
        Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive electrical pin" annotation (Placement(transformation(extent={{-110,70},{-90,90}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative electrical pin" annotation (Placement(transformation(extent={{-110,-90},{-90,-70}})));
      equation
        connect(resistor.n, inductor.p) annotation (Line(points={{-36,20},{-20,20}}, color={0,0,255}));
        connect(inductor.n, emf.p) annotation (Line(points={{0,20},{20,20},{20,10}}, color={0,0,255}));
        connect(emf.flange, inertia.flange_a) annotation (Line(points={{30,0},{40,0}}, color={0,0,0}));
        connect(inertia.flange_b, flange) annotation (Line(points={{60,0},{100,0}}, color={0,0,0}));
        connect(resistor.p, p) annotation (Line(points={{-56,20},{-80,20},{-80,80},{-100,80}}, color={0,0,255}));
        connect(emf.n, n) annotation (Line(points={{20,-10},{20,-80},{-100,-80}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Bitmap(extent={{-98,-76},{102,128}}, fileName="modelica://FM3217_2021/Resources/Images/dc-motor.jpg")}),
                                                                       Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Machine;

      model Turbine
        Modelica.Mechanics.Rotational.Sources.ConstantTorque constantTorque(tau_constant=Tt) annotation (Placement(transformation(extent={{76,-10},{56,10}})));
        Modelica.Mechanics.Rotational.Components.Inertia turbineWheel(J=Jt) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (Placement(transformation(rotation=0, extent={{-110,-10},{-90,10}})));
        parameter Modelica.SIunits.Inertia Jt=1 "Inertia of the turbine runner wheel";
        parameter Modelica.SIunits.Torque Tt=10 "Turbine torque";
      equation
        connect(constantTorque.flange, turbineWheel.flange_b) annotation (Line(points={{56,0},{10,0}}, color={0,0,0}));
        connect(flange_a, turbineWheel.flange_a) annotation (Line(points={{-100,0},{-10,0}}, color={0,0,0}));
        annotation (Icon(graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://FM3217_2021/Resources/Images/Turbine.png")}));
      end Turbine;

      model Rload
        Modelica.Electrical.Analog.Basic.Resistor resistorUpper(R=Rl/2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,18})));
        Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive electrical pin" annotation (Placement(transformation(extent={{-10,90},{10,110}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative electrical pin" annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
        parameter Modelica.SIunits.Resistance Rl=5 "Ohmic part of the load";
        Modelica.Electrical.Analog.Basic.Resistor resistorLower(R=Rl/2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-16})));
      equation
        connect(resistorUpper.p, p) annotation (Line(points={{0,28},{0,100}}, color={0,0,255}));
        connect(resistorUpper.n, resistorLower.p) annotation (Line(points={{0,8},{0,4},{1.77636e-15,4},{1.77636e-15,-6}}, color={0,0,255}));
        connect(resistorLower.n, n) annotation (Line(points={{-1.77636e-15,-26},{-1.77636e-15,-60},{0,-60},{0,-100}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Rload;

      model RLload
        extends Rload;
        Modelica.Electrical.Analog.Basic.Inductor inductor(L=Ll) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,0})));
        parameter Modelica.SIunits.Inductance Ll=5e-3 "Inductive part of the load";
      equation
        connect(inductor.p, p) annotation (Line(points={{40,10},{40,60},{0,60},{0,100}}, color={0,0,255}));
        connect(inductor.n, n) annotation (Line(points={{40,-10},{40,-60},{0,-60},{0,-100}}, color={0,0,255}));
      end RLload;

      model RLCload
        extends RLload;
        Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=Cl) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,0})));
        parameter Modelica.SIunits.Capacitance Cl=2 "Capacitive part of the load";
      equation
        connect(capacitor.p, p) annotation (Line(points={{-40,10},{-40,60},{0,60},{0,100}}, color={0,0,255}));
        connect(capacitor.n, n) annotation (Line(points={{-40,-10},{-40,-60},{0,-60},{0,-100}}, color={0,0,255}));
      end RLCload;
    end Components;

    package Tests
      model TestMachine
        extends Modelica.Icons.Example;


        Components.Machine machine annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=-90,
              origin={-30,0})));
        Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-40,-48},{-20,-28}})));
        Modelica.Blocks.Sources.Constant const(k=2)
                                               annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{-24,6},{-4,26}})));
      equation
        connect(signalVoltage.n, machine.n) annotation (Line(points={{-30,-10},{-30,-20},{0,-20},{0,-8}}, color={0,0,255}));
        connect(ground.p, signalVoltage.n) annotation (Line(points={{-30,-28},{-30,-10}}, color={0,0,255}));
        connect(const.y, signalVoltage.v) annotation (Line(points={{-49,0},{-42,0}}, color={0,0,127}));
        connect(signalVoltage.p, powerSensor.pc) annotation (Line(points={{-30,10},{-30,16},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nc, machine.p) annotation (Line(points={{-4,16},{0,16},{0,8}}, color={0,0,255}));
        connect(powerSensor.pv, powerSensor.pc) annotation (Line(points={{-14,26},{-24,26},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nv, signalVoltage.n) annotation (Line(points={{-14,6},{-14,-10},{-30,-10}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1.1));
      end TestMachine;

      model TestTurbine
        extends Modelica.Icons.Example;

        Components.Machine machine annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=-90,
              origin={-30,0})));
        Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-40,-48},{-20,-28}})));
        Modelica.Blocks.Sources.Constant const(k=2)
                                               annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{-24,6},{-4,26}})));
        Components.Turbine turbine annotation (Placement(transformation(rotation=0, extent={{56,-10},{76,10}})));
      equation
        connect(signalVoltage.n, machine.n) annotation (Line(points={{-30,-10},{-30,-20},{0,-20},{0,-8}}, color={0,0,255}));
        connect(ground.p, signalVoltage.n) annotation (Line(points={{-30,-28},{-30,-10}}, color={0,0,255}));
        connect(const.y, signalVoltage.v) annotation (Line(points={{-49,0},{-42,0}}, color={0,0,127}));
        connect(signalVoltage.p, powerSensor.pc) annotation (Line(points={{-30,10},{-30,16},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nc, machine.p) annotation (Line(points={{-4,16},{0,16},{0,8}}, color={0,0,255}));
        connect(powerSensor.pv, powerSensor.pc) annotation (Line(points={{-14,26},{-24,26},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nv, signalVoltage.n) annotation (Line(points={{-14,6},{-14,-10},{-30,-10}}, color={0,0,255}));
        connect(machine.flange, turbine.flange_a) annotation (Line(points={{20,0},{56,0}}, color={0,0,0}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1.1));
      end TestTurbine;

      model TestLoads
        extends Modelica.Icons.Example;

        Components.Machine machine annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        Components.Turbine turbine annotation (Placement(transformation(rotation=0, extent={{56,-10},{76,10}})));
        Components.Rload rload annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
        Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-30,-50},{-10,-30}})));
        Components.RLload rLload annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        Components.RLCload rLCload annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      equation
        connect(machine.flange, turbine.flange_a) annotation (Line(points={{20,0},{56,0}}, color={0,0,0}));
        connect(rload.p, machine.p) annotation (Line(points={{-20,10},{-20,20},{0,20},{0,8}}, color={0,0,255}));
        connect(rload.n, machine.n) annotation (Line(points={{-20,-10},{-20,-20},{0,-20},{0,-8}}, color={0,0,255}));
        connect(rload.n, ground.p) annotation (Line(points={{-20,-10},{-20,-30}}, color={0,0,255}));
        connect(rLload.p, machine.p) annotation (Line(points={{-50,10},{-50,20},{0,20},{0,8}}, color={0,0,255}));
        connect(rLload.n, ground.p) annotation (Line(points={{-50,-10},{-50,-20},{-20,-20},{-20,-30}}, color={0,0,255}));
        connect(rLCload.p, machine.p) annotation (Line(points={{-80,10},{-80,20},{0,20},{0,8}}, color={0,0,255}));
        connect(rLCload.n, ground.p) annotation (Line(points={{-80,-10},{-80,-20},{-20,-20},{-20,-30}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1.1));
      end TestLoads;

      model TestLoads20Nm
        extends TestLoads(turbine(Tt=20));
      end TestLoads20Nm;

      model TestLoads25Nm
        extends TestLoads(turbine(Tt=25), rLload(Rl=10, Ll=4e-3));
      end TestLoads25Nm;
    end Tests;
  end Tutorial3;

  package Tutorial4 "Electric kettle, coffee and milk"

    model ElectricKettle
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=230^2/2000, useHeatPort=true) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-20,0})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-90,-60},{-70,-40}})));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V=230*sqrt(2), freqHz=50) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,0})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{-50,10},{-30,30}})));
      Modelica.Blocks.Math.Mean mean(f=50) annotation (Placement(transformation(
            extent={{-6,-6},{6,6}},
            rotation=270,
            origin={-50,-6})));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor water(C=4.18*1000*1.7, T(
          start=283.15,
          fixed=true,
          displayUnit="degC")) annotation (Placement(transformation(extent={{0,12},{20,32}})));
      Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor annotation (Placement(transformation(extent={{22,-10},{42,10}})));
      Modelica.Electrical.Analog.Ideal.IdealOpeningSwitch switch annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
      Modelica.Blocks.Sources.Constant maxTemp(k=95) annotation (Placement(transformation(extent={{30,20},{50,40}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor kettleWall(G=5) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={10,-18})));
      Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature roomTemperature(T=21) annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
      Modelica.Blocks.Logical.OnOffController onOffController(bandwidth=3) annotation (Placement(transformation(extent={{60,-10},{80,10}})));
      Modelica.Blocks.Logical.Not invert annotation (Placement(transformation(extent={{-44,44},{-56,56}})));
    equation
      connect(sineVoltage.n, resistor.n) annotation (Line(points={{-90,-10},{-90,-20},{-20,-20},{-20,-10}}, color={0,0,255}));
      connect(ground.p, resistor.n) annotation (Line(points={{-80,-40},{-80,-20},{-20,-20},{-20,-10}}, color={0,0,255}));
      connect(powerSensor.nc, resistor.p) annotation (Line(points={{-30,20},{-20,20},{-20,10}}, color={0,0,255}));
      connect(powerSensor.pv, resistor.p) annotation (Line(points={{-40,30},{-20,30},{-20,10}}, color={0,0,255}));
      connect(powerSensor.nv, resistor.n) annotation (Line(points={{-40,10},{-40,-10},{-20,-10}}, color={0,0,255}));
      connect(powerSensor.power, mean.u) annotation (Line(points={{-50,9},{-50,1.2}}, color={0,0,127}));
      connect(water.port, resistor.heatPort) annotation (Line(points={{10,12},{10,0},{-10,0}}, color={191,0,0}));
      connect(temperatureSensor.port, water.port) annotation (Line(points={{22,0},{10,0},{10,12}}, color={191,0,0}));
      connect(sineVoltage.p, switch.p) annotation (Line(points={{-90,10},{-90,20},{-80,20}}, color={0,0,255}));
      connect(switch.n, powerSensor.pc) annotation (Line(points={{-60,20},{-50,20}}, color={0,0,255}));
      connect(kettleWall.port_a, water.port) annotation (Line(points={{10,-8},{10,12}}, color={191,0,0}));
      connect(roomTemperature.port, kettleWall.port_b) annotation (Line(points={{0,-40},{10,-40},{10,-28}}, color={191,0,0}));
      connect(temperatureSensor.T, onOffController.u) annotation (Line(points={{42,0},{50,0},{50,-6},{58,-6}}, color={0,0,127}));
      connect(maxTemp.y, onOffController.reference) annotation (Line(points={{51,30},{54,30},{54,6},{58,6}}, color={0,0,127}));
      connect(onOffController.y, invert.u) annotation (Line(points={{81,0},{90,0},{90,50},{-42.8,50}}, color={255,0,255}));
      connect(invert.y, switch.control) annotation (Line(points={{-56.6,50},{-70,50},{-70,32}}, color={255,0,255}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://FM3217_2021/Resources/Images/ElectricKettleImage.jpg")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=500, __Dymola_NumberOfIntervals=5000));
    end ElectricKettle;

    package CoffeeAndMilk
      extends Modelica.Icons.Package;

      model Coffee "Lumped thermal element storing heat"
        Modelica.SIunits.Temperature T(start = 353.15, fixed = true, displayUnit = "degC") "Temperature of element";
        Modelica.SIunits.TemperatureSlope der_T(start = 0) "Time derivative of temperature (= der(T))";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90), visible = true, iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        parameter Modelica.SIunits.SpecificHeatCapacity Ccoffee = 4180 "Specific heat capacity of the coffee.";
        parameter Modelica.SIunits.SpecificHeatCapacity Cmilk = 3770 "Specific heat capacity of the milk.";
        Modelica.SIunits.HeatCapacity C "Heat capacity of the coffee and milk mixture";
        parameter Modelica.SIunits.Volume Vcoffee(displayUnit = "ml") = 0.0002 "Volume of the coffee.";
        parameter Modelica.SIunits.Volume Vmilk(displayUnit = "ml") = 1e-05 "Volume of the added milk.";
        parameter Modelica.SIunits.Time AddTime = 300 "The time at which milk is added.";
        parameter Modelica.SIunits.Temperature MilkTemperature = 278.15 "Temperature of the added milk.";
      equation
        C = if time > AddTime then Ccoffee * Vcoffee * 1000 + Cmilk * Vmilk * 1000 else Ccoffee * Vcoffee * 1000;
        when time > AddTime then
          reinit(T, (Ccoffee * Vcoffee * 1000 * T + Cmilk * Vmilk * 1000 * MilkTemperature) / C);
        end when;
        T = port.T;
        der_T = der(T);
        C*der(T) = port.Q_flow;
        annotation (
          Icon( graphics={Text(
                origin={70,90},
                lineColor={64,64,64},
                extent={{-70,-30},{70,10}},
                textString="%C"), Bitmap(
                origin={0,0},
                extent={{-100,-100},{100,100}},
                fileName="modelica://FM3217_2020/Resources/Images/coffee.png")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},{0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,35}},
                lineColor={0,0,0},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-6,-1},{6,-12}},
                lineColor={255,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{11,13},{50,-25}},
                lineColor={0,0,0},
                textString="T"),
              Line(points={{0,-12},{0,-96}}, color={255,0,0})}),
          Documentation(info="<html>
<p>
This is a generic model for the heat capacity of a material.
No specific geometry is assumed beyond a total volume with
uniform temperature for the entire volume.
Furthermore, it is assumed that the heat capacity
is constant (independent of temperature).
</p>
<p>
The temperature T [Kelvin] of this component is a <b>state</b>.
A default of T = 25 degree Celsius (= SIunits.Conversions.from_degC(25))
is used as start value for initialization.
This usually means that at start of integration the temperature of this
component is 25 degrees Celsius. You may, of course, define a different
temperature as start value for initialization. Alternatively, it is possible
to set parameter <b>steadyStateStart</b> to <b>true</b>. In this case
the additional equation '<b>der</b>(T) = 0' is used during
initialization, i.e., the temperature T is computed in such a way that
the component starts in <b>steady state</b>. This is useful in cases,
where one would like to start simulation in a suitable operating
point without being forced to integrate for a long time to arrive
at this point.
</p>
<p>
Note, that parameter <b>steadyStateStart</b> is not available in
the parameter menu of the simulation window, because its value
is utilized during translation to generate quite different
equations depending on its setting. Therefore, the value of this
parameter can only be changed before translating the model.
</p>
<p>
This component may be used for complicated geometries where
the heat capacity C is determined my measurements. If the component
consists mainly of one type of material, the <b>mass m</b> of the
component may be measured or calculated and multiplied with the
<b>specific heat capacity cp</b> of the component material to
compute C:
</p>
<pre>
   C = cp*m.
   Typical values for cp at 20 degC in J/(kg.K):
      aluminium   896
      concrete    840
      copper      383
      iron        452
      silver      235
      steel       420 ... 500 (V2A)
      wood       2500
</pre>
</html>"));
      end Coffee;

      model CoffeeWithMilkPort "Lumped thermal element storing heat"
        Modelica.SIunits.Temperature T(start = 353.15, fixed = true, displayUnit = "degC") "Temperature of element";
        Modelica.SIunits.Enthalpy H annotation(Dialog(group = "Variables"));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90), visible = true, iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        parameter Modelica.SIunits.SpecificHeatCapacity Ccoffee = 4180 "Specific heat capacity of the coffee.";
        parameter Modelica.SIunits.SpecificHeatCapacity Cmilk = 3770 "Specific heat capacity of the milk.";
        Modelica.SIunits.HeatCapacity C "Heat capacity of the coffee and milk mixture";
        parameter Modelica.SIunits.Volume Vcoffee(displayUnit = "ml") = 0.0002 "Volume of the coffee.";
        parameter Modelica.SIunits.Temperature MilkTemperature = 278.15 "Temperature of the added milk.";
        Modelica.SIunits.Volume Vmilk(displayUnit = "ml", start = 0, fixed = true) "Volume of the added milk.";
        Modelica.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin={0,80},    extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin={0,80},    extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      equation
        C = Ccoffee * Vcoffee * 1000 + Cmilk * Vmilk * 1000;
        der(Vmilk) = u;
        T = port.T;
        T = H / C;
        der(H) = port.Q_flow + Cmilk*1000*u*MilkTemperature;
        annotation (
          Icon( graphics={Text(
                origin={70,90},
                lineColor={64,64,64},
                extent={{-70,-30},{70,10}},
                textString="%C"), Bitmap(
                origin={0,0},
                extent={{-100,-100},{100,100}},
                fileName="modelica://FM3217_2020/Resources/Images/coffee.png")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},{0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,35}},
                lineColor={0,0,0},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-6,-1},{6,-12}},
                lineColor={255,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{11,13},{50,-25}},
                lineColor={0,0,0},
                textString="T"),
              Line(points={{0,-12},{0,-96}}, color={255,0,0})}),
          Documentation(info="<html>
<p>
This is a generic model for the heat capacity of a material.
No specific geometry is assumed beyond a total volume with
uniform temperature for the entire volume.
Furthermore, it is assumed that the heat capacity
is constant (independent of temperature).
</p>
<p>
The temperature T [Kelvin] of this component is a <b>state</b>.
A default of T = 25 degree Celsius (= SIunits.Conversions.from_degC(25))
is used as start value for initialization.
This usually means that at start of integration the temperature of this
component is 25 degrees Celsius. You may, of course, define a different
temperature as start value for initialization. Alternatively, it is possible
to set parameter <b>steadyStateStart</b> to <b>true</b>. In this case
the additional equation '<b>der</b>(T) = 0' is used during
initialization, i.e., the temperature T is computed in such a way that
the component starts in <b>steady state</b>. This is useful in cases,
where one would like to start simulation in a suitable operating
point without being forced to integrate for a long time to arrive
at this point.
</p>
<p>
Note, that parameter <b>steadyStateStart</b> is not available in
the parameter menu of the simulation window, because its value
is utilized during translation to generate quite different
equations depending on its setting. Therefore, the value of this
parameter can only be changed before translating the model.
</p>
<p>
This component may be used for complicated geometries where
the heat capacity C is determined my measurements. If the component
consists mainly of one type of material, the <b>mass m</b> of the
component may be measured or calculated and multiplied with the
<b>specific heat capacity cp</b> of the component material to
compute C:
</p>
<pre>
   C = cp*m.
   Typical values for cp at 20 degC in J/(kg.K):
      aluminium   896
      concrete    840
      copper      383
      iron        452
      silver      235
      steel       420 ... 500 (V2A)
      wood       2500
</pre>
</html>"));
      end CoffeeWithMilkPort;

      model Milk
        Modelica.Blocks.Sources.Pulse pulse(
          nperiod=1,
          width=100,
          period=AddDuration,
          amplitude=AddedVolume/AddDuration,
          startTime=AddTime) annotation (Placement(visible=true, transformation(
              origin={0,0},
              extent={{-20,-20},{20,20}},
              rotation=0)));
        Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(
            visible=true,
            transformation(
              origin={60,0},
              extent={{-10,-10},{10,10}},
              rotation=0),
            iconTransformation(
              origin={0,-110},
              extent={{-10,-10},{10,10}},
              rotation=-90)));
        parameter Modelica.SIunits.Time AddTime=300 "The time at which milk is added.";
        parameter Modelica.SIunits.Time AddDuration=1 "Duration of the adding of milk.";
        parameter Modelica.SIunits.Volume AddedVolume(displayUnit="ml") = 1e-05 "Amount of milk added.";
      equation
        connect(pulse.y, y) annotation (Line(
            visible=true,
            points={{22,0},{60,0}},
            color={1,37,163}));
        annotation (Icon(graphics={Bitmap(extent={{-80,-100},{100,100}}, fileName="modelica://FM3217_2020/Resources/Images/milk.png")}));
      end Milk;

      package Scenarios
        extends Modelica.Icons.Package;

        model Approach1
          extends Modelica.Icons.Example;
          Coffee coffee(AddTime = AddTime) annotation(Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -360)));
          Modelica.Thermal.HeatTransfer.Components.ThermalConductor mug(G = 0.54) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambient(T = 293.15) annotation(Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          parameter Modelica.SIunits.Time AddTime=595   "The time at which milk is added.";
        equation
          connect(ambient.port, mug.port_b) annotation(Line(visible = true, points={{30,0},{10,0}},       color = {191, 0, 0}));
          connect(coffee.port, mug.port_a) annotation(Line(visible = true, points={{-30,0},{-10,0}},     color = {191, 0, 0}));
          annotation(experiment(StopTime=600));
        end Approach1;

        model Approach2
          extends Modelica.Icons.Example;
          CoffeeWithMilkPort coffee annotation(Placement(visible = true, transformation(origin = {-40, -0}, extent = {{-10, -10}, {10, 10}}, rotation = -360)));
          Modelica.Thermal.HeatTransfer.Components.ThermalConductor mug(G = 0.54) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambient(T = 293.15) annotation(Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          Milk milk(AddDuration = 1, AddTime=10)   annotation(Placement(visible = true, transformation(origin = {-40, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
          connect(mug.port_b, ambient.port) annotation (
            Line(points = {{10, 0}, {30, 0}, {30, 0}, {30, 0}}, color = {191, 0, 0}));
          connect(milk.y, coffee.u) annotation (
            Line(points={{-40,19},{-40,19},{-40,8},{-40,8}},          color = {0, 0, 127}));
          connect(coffee.port, mug.port_a) annotation (
            Line(points = {{-30, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {191, 0, 0}));
          connect(ambient.port, mug.port_b) annotation(Line(visible = true, points={{30,0},{10,0}},       color = {191, 0, 0}));
          connect(coffee.port, mug.port_a) annotation(Line(visible = true, points={{-30,0},{-10,0}},      color = {191, 0, 0}));
          connect(milk.y, coffee.u) annotation(Line(visible = true, points={{-40,19},{-40,8},{-40,8}},                          color = {1, 37, 163}));
        end Approach2;
        annotation(Diagram(coordinateSystem(extent = {{-148.5, -105}, {148.5, 105}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})));
      end Scenarios;
    end CoffeeAndMilk;

    package CupOfCoffee "Animation of the refilling process of a cup of coffee."
      model CupOfCoffee_1
        Real T(start=380);
      equation
        der(T) = -0.2*(T - 300);
      end CupOfCoffee_1;

      model CupOfCoffee_2 "Add some text"
        Real T(start=380) "Coffee Temperature";
      equation
        der(T) = -0.2*(T - 300) "Newton's Law of Cooling";
      end CupOfCoffee_2;

      model CupOfCoffee_3 "Use parameters"
        parameter Real T0=380 "Initial temp.";
        parameter Real Tamb=300 "Ambient temperature";
        parameter Real C=0.2;
        Real T(start=T0) "Coffee Temperature";
      equation
        der(T) = -C*(T - Tamb) "Newton's Law of Cooling";
      end CupOfCoffee_3;

      model CupOfCoffee_4 "More physical"
        import SI=Modelica.SIunits;
        parameter SI.Temperature T0=380 "Initial temp.";
        parameter SI.Temperature Tamb=300 "Ambient temperature";
        parameter SI.Density rho=1000 "Coffee density";
        parameter SI.SpecificHeatCapacity cv=4179 "Coffee specific heat";
        parameter SI.CoefficientOfHeatTransfer h=25 "Convection coefficient";
        parameter SI.Volume V=4e-4 "Volume of coffee";
        parameter SI.Area A=4e-3 "Area of coffee";
        SI.Temperature T(start=T0) "Coffee Temperature";
      equation
        rho*V*cv*der(T) = -h*A*(T - Tamb) "First law of thermodynamics";
      end CupOfCoffee_4;

      model CupOfCoffee_5
        ThermalCapacitance coffee(
          T0=380,
          rho=1000,
          V=4e-4,
          cv=4179);
        Convection cooling(
          h=25,
          A=4e-3,
          Tamb=300);
      equation
        connect(coffee.p, cooling.p);
      end CupOfCoffee_5;

      model CupOfCoffee_6
        ThermalCapacitance coffee(
          T0=380,
          rho=1000,
          V=4.08e-4,
          cv=4179);
        ThermalCapacitance cup(
          T0=300,
          rho=3700,
          V=8.45e-5,
          cv=880);
        Boundary cup2coffee(h=100, A=2.53e-2);
        Convection coffee_cooling(
          h=25,
          A=4e-3,
          Tamb=300);
        Convection cup_cooling(
          h=25,
          A=2.79e-2,
          Tamb=300);
      equation
        connect(coffee.p, cup2coffee.p1);
        connect(coffee.p, coffee_cooling.p);
        connect(cup.p, cup2coffee.p2);
        connect(cup.p, cup_cooling.p);
      end CupOfCoffee_6;

      model CupOfCoffee_7 "Additional physics"
        import SI=Modelica.SIunits;
        import Modelica.Constants.pi;
        parameter SI.CoefficientOfHeatTransfer h_air=25
          "Convection coefficient with air";
        parameter SI.CoefficientOfHeatTransfer h_fluid=100
          "Convection coefficient with fluid";
        parameter Customer customer;
        parameter Service service;
        parameter CoffeeProperties coffee;
        parameter CupProperties cup;
        SI.Temperature T "Coffee Temperature";
        SI.Temperature Tcup "Coffee cup temperature";
        SI.Mass M "Mass of coffee in the cup";
        SI.Area A "Surface area exposed to ambient";
        SI.Length H "Height of coffee in the cup";
        SI.Volume V "Volume of coffee in cup";
        Boolean drinking "true when drinking";
        Boolean empty "true when cup is empty";
        SI.MassFlowRate mdot_drink "drinking mass flow rate";
        discrete SI.MassFlowRate mdot_refill "refilling mass flow rate";
        SI.Energy U "coffee internal energy";
        SI.Area Acup_int "internal surface area of cup";
      initial equation
        H = 0.9*cup.H;
        T = service.Tcoffee;
        Tcup = service.Tamb;
      algorithm
        when sample(service.start_refill, service.dt_refill) then
          mdot_refill := service.dm_refill;
        end when;
        when H>=0.9*cup.H then
          mdot_refill := 0;
        end when;
      equation
        U = M*coffee.cv*T;
        der(U) = -h_air*A*(T - service.Tamb) - h_fluid*Acup_int*(T - Tcup) - mdot_drink*coffee.cp*
          T + mdot_refill*coffee.cp*service.Tcoffee
          "First law of thermodynamics for coffee";
        cup.M*cup.cv*der(Tcup) = -h_air*cup.A_ext*(Tcup - service.Tamb) - h_fluid*Acup_int*(
          Tcup - T) "First law of thermodynamics for cup";
        A = pi*cup.D^2/4 "Area of coffee exposed to air";
        Acup_int = pi*cup.D*(cup.H-H) "Area on inside of mug exposed to air";
        V = A*H "Volume of coffee in the cup";
        M = coffee.rho*V "Mass of coffee in the cup";
        der(M) = -mdot_drink + mdot_refill "Conservation of mass for coffee";
        empty = M <= 1e-9;
        drinking = mod(time, customer.dt_drink) <= customer.drink_duration and not empty;
        mdot_drink = if drinking then customer.sip_rate else 0;
      end CupOfCoffee_7;

      connector ThermalPort
        Modelica.SIunits.Temperature T;
        flow Modelica.SIunits.HeatFlowRate q;
      end ThermalPort;

      model ThermalCapacitance
        import SI=Modelica.SIunits;
        parameter SI.Temperature T0;
        parameter SI.Density rho;
        parameter SI.SpecificHeatCapacity cv;
        parameter SI.Volume V;
        ThermalPort p;
      initial equation
        p.T = T0;
      equation
        rho*V*cv*der(p.T) = p.q;
      end ThermalCapacitance;

      model Convection "Convection to the ambient"
        import SI=Modelica.SIunits;
        parameter SI.CoefficientOfHeatTransfer h;
        parameter SI.Temperature Tamb;
        parameter SI.Area A;
        ThermalPort p;
      equation
        p.q = h*A*(p.T - Tamb);
      end Convection;

      model Boundary "Convective boundary"
        import SI=Modelica.SIunits;
        parameter SI.CoefficientOfHeatTransfer h;
        parameter SI.Area A;
        ThermalPort p1;
        ThermalPort p2;
      equation
        p1.q + p2.q = 0 "No storage of energy";
        p1.q = h*A*(p1.T - p2.T);
      end Boundary;

      record Customer
        import SI=Modelica.SIunits;
        SI.Time dt_drink=60 "amount of time between sips";
        SI.MassFlowRate dm_drink=1.48e-2 "amount of mass consumed during each sip";
        SI.Time drink_duration=1 "duration for each drink";
        SI.MassFlowRate sip_rate=dm_drink/drink_duration "Rate of drinking coffee";
        SI.Time start_drink=1 "time when drinking starts";
      end Customer;

      record Service
        import SI=Modelica.SIunits;
        SI.Time start_refill=1750 "time when refilling starts";
        SI.Time dt_refill=1750 "amount of time between refills";
        SI.Mass dm_refill=0.3 "amount of mass added during each refill";
        SI.Time refill_duration=20 "duration for each refill";
        SI.Temperature Tcoffee=360 "Temperature of refill coffee";
        parameter SI.Temperature Tamb=300 "Ambient temperature";
      end Service;

      record CoffeeProperties
        import SI=Modelica.SIunits;
        parameter SI.Density rho=1000 "Coffee density";
        parameter SI.SpecificHeatCapacity cv=4179
          "Coffee specific heat (constant volume)";
        parameter SI.SpecificHeatCapacity cp=cv
          "Coffee specific heat (constant pressure)";
      end CoffeeProperties;

      record CupProperties
        import SI=Modelica.SIunits;
        import Modelica.Constants.pi;
        SI.Height H=0.127 "Height of coffee cup";
        SI.Diameter D=0.0508 "Diameter of cup base";
        SI.Density rho=3700 "density of the cup";
        SI.Mass M=rho*V "mass of the cup";
        SI.SpecificHeatCapacity cv=880 "Cup specific heat (constant volume)";
        SI.Length t=0.003175 "cup wall thickness";
        SI.Volume V=pi*t*H*(D + t)+pi*(D/2)^2*t "volume of the cup";
        SI.Area A_ext=pi*H*(D + 2*t) "external surface area of cup";
      end CupProperties;

      model CupOfCoffeeAnimation "animated model of cup of coffee"
        import SI=Modelica.SIunits;
        CupOfCoffee.CupOfCoffee_7 coffee_example;
        SI.Temperature coffeeTemp;
        SI.Temperature mugTemp;
        SI.Height h;
        function ComputeCoffeeColor
          input SI.Temperature T;
          output Real color[3];
        protected
          SI.Temperature Tmin=293;
          SI.Temperature Tmax=380;
          Real per=(min(max(Tmin, T), Tmax) - Tmin)/(Tmax - Tmin);
        algorithm
          color := {100 + 155*per,100 - 50*per,100 - 50*per};
        end ComputeCoffeeColor;

        function ComputeMugColor
          input SI.Temperature T;
          output Real color[3];
        protected
          SI.Temperature Tmin=293;
          SI.Temperature Tmax=380;
          Real per=(min(max(Tmin, T), Tmax) - Tmin)/(Tmax - Tmin);
        algorithm
          color := {100 + 155*per,100 - 50*per,100 - 50*per};
        end ComputeMugColor;
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape coffee(
          shapeType="cylinder",
          length=h,
          width=1.0,
          color=ComputeCoffeeColor(coffeeTemp),
          height=1.0,
          lengthDirection={0,1,0}) annotation (Placement(transformation(extent={{
                  -20,20},{0,40}}, rotation=0)));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape cup_bottom(
          shapeType="cylinder",
          length=-0.2,
          width=1.2,
          height=1.2,
          lengthDirection={0,1,0},
          color=ComputeMugColor(mugTemp));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape cup_wall(
          shapeType="pipe",
          length=1.2,
          width=1.2,
          height=1.2,
          lengthDirection={0,1,0},
          extra=0.8,
          color=ComputeMugColor(mugTemp));
      equation
        coffeeTemp = coffee_example.T;
        mugTemp = coffee_example.Tcup;
        h=coffee_example.H/coffee_example.cup.H/0.9;
        annotation (uses, experiment(StopTime=2000, __Dymola_NumberOfIntervals=5000));
      end CupOfCoffeeAnimation;

      model TestAll
        CupOfCoffee_1 cup1;
        CupOfCoffee_2 cup2;
        CupOfCoffee_3 cup3;
        CupOfCoffee_4 cup4;
        CupOfCoffee_5 cup5;
        CupOfCoffee_6 cup6;
        CupOfCoffee_7 cup7;
      end TestAll;
      annotation (
        conversion(noneFromVersion="", noneFromVersion="1"));
    end CupOfCoffee;
  end Tutorial4;

  package Tutorial5 "Introdcution to HPL"
    model ConnectingPipes
      inner HydroPower.System_HPL system_HPL(
        steadyState=true,
        Q_start(displayUnit="m3/s"),
        constantTemperature=true) annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.SinksAndSources.Fixed_pT source1(paraOption=false) annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
      HydroPower.SinksAndSources.Fixed_pT sink1(paraOption=false) annotation (Placement(transformation(extent={{60,60},{40,80}})));
      HydroPower.HydroSystems.Pipe pipe1(
        L=100,
        ZL=90,
        enable_dataVizPort_upper=false,
        enable_dataVizPort_lower=false) annotation (Placement(transformation(extent={{-8,60},{12,80}})));
    equation
      connect(source1.b, pipe1.a) annotation (Line(points={{-39,70},{-9,70}}, color={0,0,255}));
      connect(pipe1.b, sink1.b) annotation (Line(points={{13,70},{39,70}}, color={0,0,255}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=100,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end ConnectingPipes;

    model PipeWithValve
      extends ConnectingPipes;
      HydroPower.SinksAndSources.Fixed_pT source2(paraOption=false) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      HydroPower.SinksAndSources.Fixed_pT sink2(paraOption=false) annotation (Placement(transformation(extent={{60,20},{40,40}})));
      HydroPower.HydroSystems.PipeValve pipeValve2(
        m_dot_nom=115*1000,
        dp_nom=900000,
        ZL=90) annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      Modelica.Blocks.Sources.Ramp ramp(duration=10, startTime=10) annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
    equation
      connect(pipeValve2.a, source2.b) annotation (Line(points={{-11,30},{-39,30}}, color={0,0,255}));
      connect(pipeValve2.b, sink2.b) annotation (Line(points={{11,30},{39,30}}, color={0,0,255}));
      connect(ramp.y, pipeValve2.ValveCtrl) annotation (Line(points={{-79,50},{0,50},{0,41}}, color={0,0,127}));
      annotation (experiment(
          StopTime=100,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"), Documentation(info="<html>
<p>P = &rho; g Q H</p>
<p>P =100 MW</p>
<p>H = 90 m </p>
<p><br>Q = P/(&rho;*g*H) = 100e6 / (1e3 * 9.81 * 90)</p>
</html>"));
    end PipeWithValve;

    model SimpleWaterWay
      extends PipeWithValve(pipeValve2(L=110, ZL=100), ramp(
          height=-1,
          duration=1,
          offset=1));
      HydroPower.SinksAndSources.Fixed_pT source3(paraOption=false) annotation (Placement(transformation(extent={{-100,-20},{-80,0}})));
      HydroPower.SinksAndSources.Fixed_pT sink3(paraOption=false) annotation (Placement(transformation(extent={{100,-20},{80,0}})));
      HydroPower.HydroSystems.Pipe pipe3(
        horizontalIcon=true,
        L=10000,
        ZL=100,
        ZR=90) annotation (Placement(transformation(extent={{-50,-20},{-30,0}})));
      HydroPower.HydroSystems.HydroComponents.Containers.ClosedVolume closedVolume3(D=10) annotation (Placement(transformation(extent={{-10,-20},{10,0}})));
      HydroPower.HydroSystems.PipeValve pipeValve3(
        m_dot_nom=115*1000,
        dp_nom=900000,
        ZL=90) annotation (Placement(transformation(extent={{30,-20},{50,0}})));
    equation
      connect(sink3.b, pipeValve3.b) annotation (Line(points={{79,-10},{51,-10}}, color={0,0,255}));
      connect(pipeValve3.a, closedVolume3.b) annotation (Line(points={{29,-10},{10,-10}}, color={0,0,255}));
      connect(pipe3.b, closedVolume3.a) annotation (Line(points={{-29,-10},{-10,-10}}, color={0,0,255}));
      connect(pipe3.a, source3.b) annotation (Line(points={{-51,-10},{-79,-10}}, color={0,0,255}));
      connect(pipeValve3.ValveCtrl, pipeValve2.ValveCtrl) annotation (Line(points={{40,1},{40,10},{-70,10},{-70,50},{0,50},{0,41}}, color={0,0,127}));
      annotation (experiment(
          StopTime=100,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end SimpleWaterWay;

    model SimpleWaterWayWithSurgeTank
      extends SimpleWaterWay(ramp(
          height=-0.1,
          duration=1,
          offset=1,
          startTime=10), system_HPL(Q_start=119.4));
      HydroPower.SinksAndSources.Fixed_pT source4(paraOption=false) annotation (Placement(transformation(extent={{-100,-70},{-80,-50}})));
      HydroPower.HydroSystems.Pipe pipe4(
        horizontalIcon=true,
        L=10000,
        ZL=100,
        ZR=90) annotation (Placement(transformation(extent={{-50,-70},{-30,-50}})));
      HydroPower.SinksAndSources.Fixed_pT sink4(paraOption=false) annotation (Placement(transformation(extent={{100,-70},{80,-50}})));
      HydroPower.HydroSystems.PipeValve pipeValve4(
        m_dot_nom=115*1000,
        dp_nom=900000,
        ZL=90) annotation (Placement(transformation(extent={{30,-70},{50,-50}})));
      HydroPower.HydroSystems.SurgeTank surgeTank4(
        D=10,
        deltZ=100,
        Vol=100) annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
    equation
      connect(pipeValve4.b, sink4.b) annotation (Line(points={{51,-60},{79,-60}}, color={0,0,255}));
      connect(pipe4.a, source4.b) annotation (Line(points={{-51,-60},{-79,-60}}, color={0,0,255}));
      connect(pipe4.b, surgeTank4.a) annotation (Line(points={{-29,-60},{-11,-60}}, color={0,0,255}));
      connect(surgeTank4.b, pipeValve4.a) annotation (Line(points={{11,-60},{29,-60}}, color={0,0,255}));
      connect(pipeValve4.ValveCtrl, pipeValve2.ValveCtrl) annotation (Line(points={{40,-49},{40,-40},{-70,-40},{-70,50},{0,50},{0,41}}, color={0,0,127}));
      annotation (experiment(
          StopTime=1000,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end SimpleWaterWayWithSurgeTank;
  end Tutorial5;

  package Tutorial6 "Waterway with reservoir"
    model ReservoirBase
      inner HydroPower.System_HPL system_HPL(steadyState=true, constantTemperature=true) annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.HydroSystems.Reservoir headwater annotation (Placement(transformation(extent={{-70,40},{-50,60}})));
      HydroPower.HydroSystems.Reservoir tailwater annotation (Placement(transformation(extent={{50,20},{70,40}})));
      HydroPower.HydroSystems.Pipe conduit(horizontalIcon=true) annotation (Placement(transformation(extent={{-40,34},{-20,54}})));
    equation
      connect(headwater.a2_pipe, conduit.a) annotation (Line(points={{-49,44},{-41,44}}, color={0,0,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end ReservoirBase;

    model TwoReservoirs
      extends Modelica.Icons.Example;
      extends ReservoirBase(conduit(L=1000, ZL=100));
    equation
      connect(conduit.b, tailwater.a1_pipe) annotation (Line(points={{-19,44},{20,44},{20,24},{49,24}}, color={0,0,255}));
      annotation (experiment(
          StopTime=600,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end TwoReservoirs;

    model TwoReservoirsWithSource
      extends TwoReservoirs;
      HydroPower.SinksAndSources.Fixed_HT constantWaterHead(
        paraOption=false,
        H_const=headwater.H_start[1],
        Hmax=headwater.Hmax[1],
        depth=headwater.depth[1]) annotation (Placement(transformation(extent={{-100,46},{-80,66}})));
      HydroPower.SinksAndSources.Fixed_HT constantWaterTail(
        paraOption=false,
        H_const=tailwater.H_start[tailwater.n],
        Hmax=tailwater.Hmax[tailwater.n],
        depth=tailwater.depth[tailwater.n]) annotation (Placement(transformation(extent={{100,26},{80,46}})));
    equation
      connect(constantWaterHead.b, headwater.a1_open) annotation (Line(points={{-79,56},{-71,56}}, color={0,0,255}));
      connect(constantWaterTail.b, tailwater.a2_open) annotation (Line(points={{79,36},{71,36}}, color={0,0,255}));
    end TwoReservoirsWithSource;

    model WaterWayRes
      extends Modelica.Icons.Example;
      extends FM3217_2021.Tutorial6.ReservoirBase(conduit(
          L=1000,
          ZL=100,
          ZR=90));
      HydroPower.SinksAndSources.Fixed_HT constantWaterHead(
        paraOption=false,
        H_const=75,
        Hmax=headwater.Hmax[1],
        depth=headwater.depth[1]) annotation (Placement(transformation(extent={{-100,46},{-80,66}})));
      HydroPower.SinksAndSources.Fixed_HT constantWaterTail(
        paraOption=false,
        H_const=75,
        Hmax=tailwater.Hmax[tailwater.n],
        depth=tailwater.depth[tailwater.n]) annotation (Placement(transformation(extent={{100,26},{80,46}})));
      HydroPower.HydroSystems.PipeValve pipeValve(
        m_dot_nom=113.26*1000,
        dp_nom=900000,
        d_nom(displayUnit="kg/m3"),
        ZL=90) annotation (Placement(transformation(extent={{20,30},{40,50}})));
      HydroPower.HydroSystems.SurgeTank surgeTank(
        D=10,
        deltZ=50,
        Vol=100)                                                      annotation (Placement(transformation(extent={{-10,34},{10,54}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.9,                     duration=30,
        offset=0.1,
        startTime=100)                                           annotation (Placement(transformation(extent={{-10,70},{10,90}})));
    equation
      connect(constantWaterHead.b, headwater.a1_open) annotation (Line(points={{-79,56},{-71,56}}, color={0,0,255}));
      connect(tailwater.a2_open, constantWaterTail.b) annotation (Line(points={{71,36},{79,36}}, color={0,0,255}));
      connect(surgeTank.b, pipeValve.a) annotation (Line(points={{11,44},{16,44},{16,40},{19,40}}, color={0,0,255}));
      connect(pipeValve.b, tailwater.a1_pipe) annotation (Line(points={{41,40},{44,40},{44,24},{49,24}}, color={0,0,255}));
      connect(conduit.b, surgeTank.a) annotation (Line(points={{-19,44},{-11,44}}, color={0,0,255}));
      connect(ramp.y, pipeValve.ValveCtrl) annotation (Line(points={{11,80},{30,80},{30,51}}, color={0,0,127}));
    end WaterWayRes;

    model WaterWayResClosingValve
      extends FM3217_2021.Tutorial6.WaterWayRes(ramp(height=-0.9,
                                                               offset=1));
    end WaterWayResClosingValve;

    model SundsbarmWaterway
      extends Modelica.Icons.Example;
      inner HydroPower.System_HPL system_HPL(steadyState=true,
        Q_start=24,                                            constantTemperature=true) annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.HydroSystems.Reservoir headwater(
        Hmax=ones(headwater.n)*(564 + 48 + 5),
        depth=ones(headwater.n)*(48 + 5),
        H_start=ones(headwater.n)*(564 + 48))     annotation (Placement(transformation(extent={{-90,20},{-70,40}})));
      HydroPower.HydroSystems.Reservoir tailwater(
        Hmax=ones(tailwater.n)*(110 + 5 + 3),
        depth=ones(tailwater.n)*(5 + 3),
        H_start=ones(tailwater.n)*(110 + 5))      annotation (Placement(transformation(extent={{70,0},{90,20}})));
      HydroPower.HydroSystems.Pipe conduit(
        endD={5.8,5.8},
        ZL=564,
        ZR=541.5,
        horizontalIcon=true,
        L=6600)                                                         annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      HydroPower.SinksAndSources.Fixed_HT constantWaterHead(
        paraOption=false,
        H_const=564 + 48,
        Hmax=564 + 48 + 5,
        depth=48 + 5)
                  annotation (Placement(transformation(extent={{-70,50},{-90,70}})));
      HydroPower.SinksAndSources.Fixed_HT constantTailWater(
        paraOption=false,
        H_const=110 + 5,
        Hmax=110 + 5 + 3,
        depth=5 + 3)
                  annotation (Placement(transformation(extent={{70,28},{90,48}})));
      HydroPower.HydroSystems.SurgeTank surgeTank(
        D=3.6,
        deltZ=150,
        H2L=0.87,
        Vol=100) annotation (Placement(transformation(extent={{-30,14},{-10,34}})));
      HydroPower.HydroSystems.PipeValve pressureShaft(
        endD={3,3},
        m_dot_nom=24e3,
        dp_nom=4890000,
        L=724,
        ZL=541.5,
        ZR=112.5) annotation (Placement(transformation(extent={{0,10},{20,30}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=-0.9,
        duration=10,
        offset=1,
        startTime=100) annotation (Placement(transformation(extent={{-20,50},{0,70}})));
      HydroPower.HydroSystems.Pipe tailRace(
        endD={5.8,5.8},
        ZL=110.5,
        ZR=110,
        horizontalIcon=true,
        L=600)                                                          annotation (Placement(transformation(extent={{40,-6},{60,14}})));
      HydroPower.HydroSystems.HydroComponents.Containers.ClosedVolume turbineHouse(D=5.8, L=2) annotation (Placement(transformation(extent={{24,6},{36,18}})));
    equation
      connect(headwater.a2_pipe,conduit. a) annotation (Line(points={{-69,24},{-61,24}}, color={0,0,255}));
      connect(headwater.a1_open, constantWaterHead.b) annotation (Line(points={{-91,36},{-96,36},{-96,60},{-91,60}}, color={0,0,255}));
      connect(tailwater.a2_open, constantTailWater.b) annotation (Line(points={{91,16},{96,16},{96,38},{91,38}}, color={0,0,255}));
      connect(conduit.b, surgeTank.a) annotation (Line(points={{-39,24},{-31,24}}, color={0,0,255}));
      connect(surgeTank.b, pressureShaft.a) annotation (Line(points={{-9,24},{-4,24},{-4,20},{-1,20}}, color={0,0,255}));
      connect(pressureShaft.ValveCtrl, ramp.y) annotation (Line(points={{10,31},{10,60},{1,60}}, color={0,0,127}));
      connect(tailwater.a1_pipe, tailRace.b) annotation (Line(points={{69,4},{61,4}}, color={0,0,255}));
      connect(tailRace.a, turbineHouse.b) annotation (Line(points={{39,4},{38,4},{38,12},{36,12}}, color={0,0,255}));
      connect(pressureShaft.b, turbineHouse.a) annotation (Line(points={{21,20},{22,20},{22,12},{24,12}}, color={0,0,255}));
      annotation (
        experiment(
          StopTime=600,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end SundsbarmWaterway;
  end Tutorial6;

  package Tutorial7 "First model of Sundsbarm"
    model PlantConnectAndDisconnectToGrid
      extends HydroPower.Examples.PlantConnectAndDisconnectToGrid;
      annotation (experiment(
          StopTime=600,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end PlantConnectAndDisconnectToGrid;

    model Sundsbarm "Model of the Sundsbarm Power Station"
      extends Modelon.Icons.Experiment;

      HydroPower.MechanicalSystems.BasicTurbine turbine(
        np=12,
        H_nom=480,
        tableOnFile=true,
        LdraftTube=10,
        DavDraftTube=2,
        LscrollCase=5,
        DavScrollCasing=2,
        PUInFlowTables=true,
        QTableName="Qtab",
        Q_nom=24,
        H_start=564 + 48,
        H_start_draftTube=115,
        Ty=0.4,
        yvLim1=[-0.1,0.1],
        yvLim2=[-0.2,0.2],
        TurbineDataFile=Modelica.Utilities.Files.loadResource(HydroPower.TABLE_DIR
             + "TurbineDataFile.mat"),
        P_nom=103000000)
                        annotation (Placement(transformation(extent={{20,-50},{40,-30}},
                     rotation=0)));

      HydroPower.ElectricalSystems.PowerGrid powerGrid(
        startTime=1e6,
        distNoGen={-2,0,0},
        distTgen={150,1e6,1e6}) annotation (Placement(transformation(extent={{-90,40},{-70,60}}, rotation=0)));

      Modelica.Blocks.Sources.Ramp pwr_ref(
        duration=10,
        height=0,
        offset=45e6,
        startTime=1e6) annotation (Placement(transformation(extent={{-37,64},{-25,76}},
                      rotation=0)));
      HydroPower.ElectricalSystems.GeneratorAndMCB generator(
        np={12},
        f_start=0,
        J={83e3},
        timeMCB_close={150},
        timeMCB_open={200},
        P_nom={103000000})
                          annotation (Placement(transformation(extent={{-60,40},{-40,60}},
                       rotation=0)));
      HydroPower.ControllersAndSensors.TurbineGovernorAnalog turbineGovernor(
        K_noLoad=0.8,
        Ki_noLoad=0.2,
        Kd_noLoad=0.1,
        ep=1,
        tRamp=40,
        P_generator_nom=generator.P_nom[1],
        enableRamp=false,
        K_load=0.4,
        Kd_load=0.1,
        Ki_load=0.2)      annotation (Placement(transformation(extent={{-24,40},{-4,60}}, rotation=0)));
      HydroPower.Visualizers.RealValue turbinePower(precision=2, input_Value=turbine.summary.P_turbine*1e-6) annotation (Placement(transformation(extent={{64,10},{78,24}})));
      HydroPower.Visualizers.BooleanIndicator MCB(input_Value=turbineGovernor.summary.isMCB) annotation (Placement(transformation(extent={{64,73},{77,87}})));
      HydroPower.Visualizers.RealValue gridbalanceNum(precision=2, input_Value=generator.summary.P_grid_tot*1e-6) annotation (Placement(transformation(extent={{64,38},{78,52}})));
      inner HydroPower.System_HPL system_HPL(
        steadyState=true, constantTemperature=true)
                     annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.Visualizers.RealValue gridbalanceNum1(precision=2, input_Value=generator.summary.f[1]) annotation (Placement(transformation(extent={{64,53},{78,67}})));
      HydroPower.Visualizers.RealValue gridbalanceNum2(precision=2, input_Value=generator.summary.P_generator[1]*1e-6) annotation (Placement(transformation(extent={{64,25},{78,39}})));
      HydroPower.HydroSystems.Pipe      pressureShaft(
        L=724,
        ZR=112.5,
        endD={3,3},
        ZL=541.5,
        enable_dataVizPort_lower=false)
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
      HydroPower.HydroSystems.Pipe conduit(
        horizontalIcon=true,
        ZL=564,
        endD={5.8,5.8},
        L=6600,
        ZR=541.5)
        annotation (Placement(transformation(extent={{-70,-50},{-50,-30}})));
      HydroPower.HydroSystems.SurgeTank surgeTank(
        deltZ=150,
        H2L=0.87,
        Vol=100,
        D=3.6)
        annotation (Placement(transformation(extent={{-40,-50},{-20,-30}})));
      HydroPower.HydroSystems.Reservoir reservoir(
        depth=ones(reservoir.n)*(48 + 5),
        H_start=fill((564 + 48), reservoir.n),
        Hmax=fill((564 + 48 + 5), reservoir.n)) annotation (Placement(transformation(extent={{-95,-44},{-75,-24}})));
      HydroPower.HydroSystems.Reservoir river(
        depth=ones(river.n)*(5 + 3),
        H_start=fill((115), river.n),
        Hmax=fill((110 + 5 + 3), river.n)) annotation (Placement(transformation(extent={{75,-44},{95,-24}})));
      HydroPower.HydroSystems.Pipe downstream(
        horizontalIcon=true,
        L=600,
        ZL=110.5,
        endD={5.8,5.8},
        ZR=110) annotation (Placement(transformation(extent={{50,-50},{70,-30}})));
      HydroPower.SinksAndSources.Fixed_HT waterSource(
        paraOption=false,
        H_const=reservoir.H_start[reservoir.n],
        Hmax=reservoir.Hmax[reservoir.n],
        depth=reservoir.depth[1])
                  annotation (Placement(transformation(extent={{-70,-18},{-90,2}})));
      HydroPower.SinksAndSources.Fixed_HT waterSink(
        paraOption=false,
        H_const=river.H_start[river.n],
        Hmax=river.Hmax[river.n],
        depth=river.depth[river.n])
        annotation (Placement(transformation(extent={{70,-18},{90,2}})));
    equation

      connect(powerGrid.f_grid, generator.f_grid) annotation (Line(
          points={{-69,43},{-61,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(powerGrid.P_grid_balance, generator.P_grid_balance) annotation (Line(
          points={{-69,57},{-61,57}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pwr_ref.y, turbineGovernor.P_reference) annotation (Line(
          points={{-24.4,70},{-20,70},{-20,61}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.f_out[1], turbineGovernor.f) annotation (Line(
          points={{-39,43},{-25,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.onMCB, powerGrid.MCB) annotation (Line(
          points={{-50,61},{-50,70},{-80,70},{-80,61}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.onMCB[1], turbineGovernor.isMCB) annotation (Line(
          points={{-50,61},{-50,80},{-14,80},{-14,61}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.P_out[1], turbineGovernor.P_generator) annotation (Line(
          points={{-39,57},{-25,57}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(generator.f_out, powerGrid.f) annotation (Line(
          points={{-39,43},{-33,43},{-33,6},{-94,6},{-94,43},{-91,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(powerGrid.J_grid, generator.J_grid)
        annotation (Line(points={{-69,50},{-61,50}},          color={0,0,127}));
      connect(turbineGovernor.y, turbine.yGV)
        annotation (Line(points={{-3,50},{36,50},{36,-29}}, color={0,0,127}));
      connect(generator.f_out[1], turbine.f_generator) annotation (Line(points={{-39,43},{-33,43},{-33,6},{24,6},{24,-29}},
                                                         color={0,0,127}));
      connect(turbine.TurbineData[1], generator.P_turbine[1]) annotation (Line(
            points={{30,-29.6667},{30,34},{-56,34},{-56,39}},            color={0,0,
              127}));
      connect(surgeTank.a,conduit. b) annotation (Line(
          points={{-41,-40},{-49,-40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(surgeTank.b,pressureShaft. a) annotation (Line(
          points={{-19,-40},{-11,-40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(reservoir.a2_pipe,conduit. a) annotation (Line(
          points={{-74,-40},{-71,-40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(downstream.b,river. a1_pipe) annotation (Line(points={{71,-40},{74,-40}},
                                  color={0,0,255}));
      connect(waterSource.b,reservoir. a1_open) annotation (Line(points={{-91,-8},{-100,-8},{-100,-28},{-96,-28}},
                                      color={0,0,255}));
      connect(waterSink.b,river. a2_open) annotation (Line(points={{91,-8},{100,-8},{100,-28},{96,-28}},
                                       color={0,0,255}));
      connect(turbine.a, pressureShaft.b) annotation (Line(points={{19,-40},{11,-40}}, color={0,0,255}));
      connect(turbine.b, downstream.a) annotation (Line(points={{41,-40},{49,-40}}, color={0,0,255}));
      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{50,91},{90,7}},
              lineColor={215,215,215},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              radius=2),           Text(
              extent={{54,13},{90,9}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="turbine power [MW]"),Text(
              extent={{57,73},{85,65}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="Main Circuit Breaker"),
                                Text(
              extent={{39,41},{103,37}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="grid power balance [MW]"),
                                Text(
              extent={{39,56},{103,52}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="generator frequency [Hz]"),
                                Text(
              extent={{38,27},{102,23}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="generator power [MW]"),
            Bitmap(extent={{-11,-99},{77,-56}}, fileName="modelica://FM3217_2021/Resources/Images/SundsbarmTurbine.png"),
            Bitmap(extent={{-90,-15},{20,58}}, fileName="modelica://FM3217_2021/Resources/Images/SundsbarmGenerator.png")}),
        experiment(
          StopTime=600,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"),
        Documentation(info="<html>
<p>This example illustrates a hydro power plant acting under no load and when connected to the power grid. </p>
<h4>Model experiment description</h4>
<p>When there is no load present the governor will only have the frequency error as input signal which will have the effect that the frequency of the hydro plant generator is controlled to equal the nominal frequency. This behaviour can be seen during the first 150s of simulation. </p>
<p>When 150s has passed and the frequency of the generator is synchronized to the grid frequency, the MCB is closed, the power reference is set to 45MW and new PID parameters are applied. </p>
<p>At time=350 load rejection takes place and the MCB opens once again. </p>
<h4>Simulation setup</h4>
<p>Simulate for 600s using solver Radau with a tolerance set to 1e-6.</p>
<h4>Output</h4>
<p>The most interesting variables are:</p>
<ul>
<li>generator frequency - generator.summary.f[1]</li>
<li>generated power - generator.summary.P_generator[1]</li>
</ul>
</html>",     revisions="<html>
<hr><p><font color=\"#E72614\"><b>Copyright &copy; 2004-2017, MODELON AB</b></font> <font color=\"#AFAFAF\">The use of this software component is regulated by the licensing conditions for the HydroPower Library. This copyright notice must, unaltered, accompany all components that are derived from, copied from, or by other means have their origin from the HydroPower Library. </font>
</html>"));
    end Sundsbarm;
  end Tutorial7;

  package Tutorial8 "Power flow calculations"
    model EquivalentPowerFlow
      Modelica.Mechanics.Rotational.Components.Inertia grid(J=100, w(displayUnit="Hz", start=314.15926535898)) annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
      Modelica.Mechanics.Rotational.Sources.ConstantTorque generator(tau_constant=10) annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      Modelica.Mechanics.Rotational.Sources.ConstantTorque load(tau_constant=-15) annotation (Placement(transformation(extent={{90,-10},{70,10}})));
    equation
      connect(generator.flange, grid.flange_a) annotation (Line(points={{-70,0},{-40,0}}, color={0,0,0}));
      connect(load.flange, grid.flange_b) annotation (Line(points={{70,0},{40,0}}, color={0,0,0}));
    end EquivalentPowerFlow;

    model FrequencyCalc
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J=2*200e6/(2*Modelica.Constants.pi*50)^2, w(
          displayUnit="Hz",
          fixed=true,
          start=314.15926535898)) annotation (Placement(transformation(extent={{-62,70},{-42,90}})));
      Modelica.Mechanics.Rotational.Sources.ConstantTorque load(tau_constant=-10e6/(2*Modelica.Constants.pi*50)) annotation (Placement(transformation(extent={{80,70},{60,90}})));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{10,70},{30,90}})));
      Modelica.Blocks.Continuous.Integrator energy annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={12,50})));
    equation
      connect(inertia.flange_b, powerSensor.flange_a) annotation (Line(points={{-42,80},{10,80}}, color={0,0,0}));
      connect(powerSensor.flange_b, load.flange) annotation (Line(points={{30,80},{60,80}}, color={0,0,0}));
      connect(energy.u, powerSensor.power) annotation (Line(points={{12,62},{12,69}}, color={0,0,127}));
      annotation (Diagram( graphics={Text(
              extent={{-74,66},{-22,56}},
              lineColor={238,46,47},
              textString="= 200 MJ @50Hz"), Text(
              extent={{6,74},{128,22}},
              lineColor={238,46,47},
              textString="= 10 MW 
load for 1 sec
= 10 MJ")}),
        Documentation(info="<html>
<h4>Example</h4>
<p>K<sub>1</sub> = 200 MJ</p>
<p>P<sub>L</sub> = 10 MW @ 1sec</p>
<h5>Question:</h5>
<p>If the Inertia was rotating at 50 Hz in the beginning what is the speed/frequency after 1 sec?</p>
<h5>Solution:</h5>
<p><br>K<sub>2</sub> = K<sub>1</sub> - P<sub>L</sub> &middot; 1 s = 200 MJ - 10 MW &middot; 1 s = 200 MJ - 10 MJ = 190 MJ</p>
<p><br>K = 1/2 &middot; J &middot; w<sup>2</sup> = 2 &middot; J &middot; &pi;<sup>2</sup> &middot; f<sup>2</sup></p>
<p><br>K<sub>1</sub>/K<sub>2</sub> = 200 MJ / 190 MJ = f<sub>1<sup>2</sup>/f<sub>2<sup>2</sup> </p>
<blockquote>
<p><br>&rArr; f<sub>2</sub> = &radic;(50<sup>2</sup> &middot; 190/200) = <b>48.73 Hz</b></p>
</blockquote>
<h5>Tip:</h5>
<p>J = 2 &middot; K /w<sup>2</sup></p>
<p>T = P/w</p>
</html>"));
    end FrequencyCalc;

    model FrequencyCalcCorrect
      extends FrequencyCalc;
      Modelica.Mechanics.Rotational.Components.Inertia inertia1(J=2*200e6/(2*Modelica.Constants.pi*50)^2, w(
          displayUnit="Hz",
          start=314.15926535898,
          fixed=true)) annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor1 annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
      Modelica.Mechanics.Rotational.Sources.Torque torque annotation (Placement(transformation(extent={{30,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant const(k=-10e6) annotation (Placement(transformation(extent={{92,0},{72,20}})));
      Modelica.Blocks.Math.Division division annotation (Placement(transformation(extent={{56,-6},{44,6}})));
      Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-50,-20})));
      Modelica.Blocks.Continuous.Integrator energy1 annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-18,-30})));
    equation
      connect(inertia1.flange_b,powerSensor1. flange_a) annotation (Line(points={{-60,0},{-20,0}},     color={0,0,0}));
      connect(powerSensor1.flange_b,torque. flange) annotation (Line(points={{0,0},{10,0}},      color={0,0,0}));
      connect(const.y,division. u1) annotation (Line(points={{71,10},{66,10},{66,3.6},{57.2,3.6}},       color={0,0,127}));
      connect(torque.tau,division. y) annotation (Line(points={{32,0},{43.4,0}},     color={0,0,127}));
      connect(inertia1.flange_b,speedSensor. flange) annotation (Line(points={{-60,0},{-50,0},{-50,-10}},     color={0,0,0}));
      connect(speedSensor.w,division. u2) annotation (Line(points={{-50,-31},{-50,-60},{66,-60},{66,-3.6},{57.2,-3.6}},   color={0,0,127}));
      connect(powerSensor1.power, energy1.u) annotation (Line(points={{-18,-11},{-18,-18}}, color={0,0,127}));
    end FrequencyCalcCorrect;

    model ElectricalPowerFlowCalc
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      Modelica.Electrical.Analog.Basic.Inductor Xs annotation (Placement(transformation(extent={{-20,60},{0,80}})));
      Modelica.Electrical.Analog.Sources.SineVoltage Ea(
        V=sqrt(2)*230,
        phase=0,
        freqHz=50) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-60,50})));
      Modelica.Electrical.Analog.Sources.SineVoltage Vterminal(V=sqrt(2)*230, freqHz=50) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={60,50})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Modelica.Electrical.MultiPhase.Basic.Inductor inductor(L={10,10,10}) annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
      Modelica.Electrical.MultiPhase.Sources.SineVoltage EaM(
        V=0.9*sqrt(2)*{230,230,230},
        phase=-Modelica.Electrical.MultiPhase.Functions.symmetricOrientation(3) + 0*{30,30,30}/180*Modelica.Constants.pi,
        freqHz={50,50,50}) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,-40})));
      Modelica.Electrical.MultiPhase.Sources.SineVoltage VterminalM(V=sqrt(2)*{230,230,230}, freqHz={50,50,50}) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={60,-40})));
      Modelica.Electrical.MultiPhase.Basic.Star star annotation (Placement(transformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-60,-75})));
      Modelica.Electrical.Analog.Basic.Ground ground1
                                                     annotation (Placement(transformation(extent={{-70,-106},{-50,-86}})));
      Modelica.Electrical.MultiPhase.Sensors.AronSensor P annotation (Placement(transformation(extent={{10,-30},{30,-10}})));
      Modelica.Electrical.MultiPhase.Sensors.ReactivePowerSensor Q annotation (Placement(transformation(extent={{34,-30},{54,-10}})));
    equation
      connect(Xs.p, Ea.p) annotation (Line(points={{-20,70},{-60,70},{-60,60}}, color={0,0,255}));
      connect(Ea.n, ground.p) annotation (Line(points={{-60,40},{-60,36},{-60,36},{-60,20}},
                                                                           color={0,0,255}));
      connect(Vterminal.n, ground.p) annotation (Line(points={{60,40},{60,20},{-60,20}}, color={0,0,255}));
      connect(Xs.n, powerSensor.pc) annotation (Line(points={{0,70},{20,70}}, color={0,0,255}));
      connect(powerSensor.nc, Vterminal.p) annotation (Line(points={{40,70},{60,70},{60,60}}, color={0,0,255}));
      connect(powerSensor.nv, ground.p) annotation (Line(points={{30,60},{30,20},{-60,20}}, color={0,0,255}));
      connect(powerSensor.pv, Vterminal.p) annotation (Line(points={{30,80},{30,90},{60,90},{60,60}}, color={0,0,255}));
      connect(EaM.plug_p, inductor.plug_p) annotation (Line(points={{-60,-30},{-60,-20},{-20,-20}}, color={0,0,255}));
      connect(star.plug_p, EaM.plug_n) annotation (Line(points={{-60,-65},{-60,-50}}, color={0,0,255}));
      connect(VterminalM.plug_n, EaM.plug_n) annotation (Line(points={{60,-50},{60,-60},{-60,-60},{-60,-50}}, color={0,0,255}));
      connect(star.pin_n, ground1.p) annotation (Line(points={{-60,-85},{-60,-86}}, color={0,0,255}));
      connect(inductor.plug_n, P.plug_p) annotation (Line(points={{0,-20},{10,-20}}, color={0,0,255}));
      connect(P.plug_n, Q.plug_p) annotation (Line(points={{30,-20},{34,-20}}, color={0,0,255}));
      connect(Q.plug_n, VterminalM.plug_p) annotation (Line(points={{54,-20},{60,-20},{60,-30}}, color={0,0,255}));
      annotation (experiment(StopTime=0.1));
    end ElectricalPowerFlowCalc;
  end Tutorial8;

  package Tutorial9

    model SampleCalculations
      extends Modelica.Icons.Information;
      annotation (preferredView="info",Documentation(info="<html>
<h4>Sample calculations for production balance in the Power Grid</h4>
<h5>Default values</h5>
<p>P_grid = 1000 MW</p>
<p>P_1 = 0.45 * P_grid = 450 MW</p>
<p>P_1units = 20</p>
<p>P_1perUnit = P_1 / P_1units = 450 MW / 20 = 22.5 MW</p>
<p><u><b>&Delta;P</b></u> = distNoGen * P_1perUnit = -2 * 22.5 =<b> <u>45 MW</b></u></p>
<h5>Manual Droop calculations for DroopSimulation</h5>
<p>R = - (&Delta;F/F_r) / (&Delta;P/P_r)</p>
<p>&Delta;P = - P_r * (&Delta;F/F_r) / R</p>
<p>&Delta;P = - 45 MW * (-0.4043/50) / 0.1</p>
<p><u>&Delta;P = 3.639 MW (theory)</u></p>
<p><u>&Delta;P = 3.167 MW (simulation)</u></p>
</html>"));
    end SampleCalculations;

    model ResitiveLoad
      extends HydroPower.Examples.PlantConnectAndDisconnectToGrid(
        powerGrid(
          loadDiv={1,0,0},
          enableDroop=false,
          distTgen={150,1000,1e6},
          distNoGen={-1,-50,0}),
        pwr_ref(offset=45e6),
        turbineGovernor(enableDroop=false),
        generator(timeMCB_open={1e6}));

      annotation (experiment(
          StopTime=2000,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"),                                                    preferredView="info", Documentation(info="<html>
<p>The model demonstrates the behaviour of an electrical generator in combination with an electrical grid where we only have pure resistive loads present. </p>
<h4>Task</h4>
<p>Apply the following changes:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor: </li>
<li><ul>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid: </li>
<li><ul>
<li>Only <b>resistive</b> loads: <span style=\"font-family: monospace;\">loadDiv = {1,0,0}</span></li>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li>Remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>Remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s</span> from production</li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>How does the frequency behave?</li>
<li>What is the power difference &Delta;P of the generator when comparing the power produced at <span style=\"font-family: monospace;\">t=1000s</span> and the power produced by the generator at <span style=\"font-family: monospace;\">t=2000s</span>?</li>
<li>Change the power set point of the turbine governor to 45 MW and look at the frequency behaviour.</li>
</ol>
</html>"));
    end ResitiveLoad;

    model DroopSimulations
      extends ResitiveLoad(
        pwr_ref(offset=22.5e6),
        turbineGovernor(enableDroop=true, ep=0.05),
        powerGrid(P_grid=100000000, distNoGen={-10,-50,0}));
      annotation (experiment(
          StopTime=10000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"),                                                   preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Enable the droop: <span style=\"font-family: monospace;\">enableDroop = true</span> </li>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid: </li>
<li><ul>
<li>Only resistive loads: loadDiv = {1,0,0}</li>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li>Remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>Remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s </span>from production</li>
</ul></li>
</ul>
<p><br>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine governor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>What is the power difference &Delta;P of the generator when comparing the power produced at <span style=\"font-family: monospace;\">t=1000s</span> and the power produced by the generator at <span style=\"font-family: monospace;\">t=2000s</span>?</li>
<li>What is &Delta;f in this case and does &Delta;P fit with the theoretical obtained one from the calculation with the droop factor <span style=\"font-family: monospace;\">ep</span>?</li>
<li>Now re-run the simulation but set the droop for the turbine governor to <span style=\"font-family: monospace;\">ep = 0.05</span>. What difference does this make?</li>
<li>What does the <span style=\"font-family: monospace;\">DeadBand</span> setting in the turbine governor mean and what difference does it make to simulate with <span style=\"font-family: monospace;\">DeadBand=0</span>?</li>
<li>Try the simulation for a smaller power grid, e.g., <span style=\"font-family: monospace;\">P_grid = 100 MW</span>. Pay attention that you also need to adjust the number of units that are changed in order to maintain the 22.5 MW change and make the second change only 22.5 MW. How does the result differ now?</li>
</ol>
</html>"),
        experiment(StopTime=2000, Algorithm="Radau"));
    end DroopSimulations;
  end Tutorial9;

  package Tutorial10
    model LoadChanges
      extends HydroPower.Examples.PlantConnectAndDisconnectToGrid;
      annotation (preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid:</li>
<li><ul>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li><b>resistive, f and f*f dependent loads </b>(default <span style=\"font-family: monospace;\">loadDiv={0.5,0.25,0.25}</span>)</li>
<li>remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s </span>from production</li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine govenor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>Compare the behaviour of the frequency with the one obtained from model a load disttribution of <span style=\"font-family: monospace;\">loadDiv={1,0,0}</span>?</li>
<li>In case you do not see much difference run another simulation with <span style=\"font-family: monospace;\">loadDiv={0,0,1}</span> and compare again the frequency with the previous simulations.</li>
</ol>
</html>"),
        experiment(
          StopTime=2000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end LoadChanges;

    model ProductionDroop
      extends LoadChanges;
      annotation (preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid:</li>
<li><ul>
<li>enable the droop for the production: <span style=\"font-family: monospace;\">enableDroop = true</span></li>
<li><b>set the droop of the production to <span style=\"font-family: monospace;\">ep = {0.1,0.08,0.04}</span></b></li>
<li>default resistive, f and f*f dependent loads: <span style=\"font-family: monospace;\">loadDiv={0.5,0.25,0.25}</span></li>
<li>remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s </span>from production</li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine govenor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>Compare the behaviour of the frequency with the one with the droop disabled?</li>
<li>What does the power balance (<code>P_grid_tot</code>) look like compared to no droop within the power grid?</li>
</ol>
</html>"),
        experiment(
          StopTime=2000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end ProductionDroop;

    model RandomLoad
      extends LoadChanges;
                                                                   annotation (preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid:</li>
<li><ul>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li>default resistive, f and f*f dependent loads: <span style=\"font-family: monospace;\">loadDiv={0.5,0.25,0.25}</span></li>
<li>remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li><b>Add a random load disturbance with a time interval of <span style=\"font-family: monospace;\">h = 10s</span> after 1000 seconds of simulation.</b></li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine governor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f</span></li>

</ul>

<h4>Questions</h4>
<ol>
<li>Look at the power balance <code>P_grid_tot</code> and try to explain the behaviour!</li>
<li>What happens when you change the interval <code>h</code> for example to <code>h = 1s</code>?</li>
<li>What happens when you enable droop in the Power Grid again? E.g., <code>ep = {0.1,0.08,0.04}</code>?</li>
</ol>
</html>"),
        experiment(
          StopTime=2000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end RandomLoad;
  end Tutorial10;
  
  package Tutorial11
    record Data "Data record for all simulation parameter"
      extends Modelica.Icons.Record;
  
      // General
      parameter Modelica.SIunits.Power Pn=45e6 "Nominal turbine power"
        annotation(Dialog(group="General"));
      parameter Integer np(min=2)=12
                                    "Number of poles of generator"
        annotation(Dialog(group="General"));
  
      // Mechanical System
      parameter Modelica.SIunits.Inertia J=183e3 "Turbine and generator inertia"
        annotation(Dialog(group="Mechanical System"));
      parameter Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm rpm_n=3000/np*2 "Nominal turbine speed"
        annotation(Dialog(group="Mechanical System", enable=false));
      parameter Modelica.SIunits.AngularVelocity wn = rpm_n / 60 * 2 * Modelica.Constants.pi "Nominal angular velocity"
        annotation(Dialog(group="Mechanical System", enable=false));
        parameter Modelica.SIunits.Time Ta=J*wn^2/Pn "Mechanical time constant"
        annotation(Dialog(group="Mechanical System", enable=false));
  
         // WaterWay
      parameter Modelica.SIunits.VolumeFlowRate Q_n=92 "Nominal flow rate"
        annotation (Dialog(group="WaterWay"));
      parameter Modelica.SIunits.Length H_n=50 "Nominal water head"
        annotation (Dialog(group="WaterWay"));
        parameter Modelica.SIunits.Length L[2]={200,100} "Total length of the waterway"
        annotation (Dialog(group="WaterWay"));
      parameter Modelica.SIunits.Length d[2]= {5.5,5.5} "Average pipe diameter"
        annotation (Dialog(group="WaterWay"));
  
        parameter Modelica.SIunits.Area A[2]=Modelica.Constants.pi/4*d .^ 2 "Average pipe area"
        annotation (Dialog(group="WaterWay", enable=false));
  
      parameter Modelica.SIunits.Time Tw=Q_n/(Modelica.Constants.g_n*H_n)*(L[1]/A[1]+L[2]/A[2])
        "Time constant of the waterway"
        annotation (Dialog(group="WaterWay", enable=false));
  
      annotation (preferredView="info",Documentation(info="<html>
  <p>
  The transfer function of the hydro power system can be represented by:
  </p>
  <img src=\"modelica://FM3217_2021/Resources/Images/TFsystem.png\">
  <p>For the determining the time constants of the transfer-function based model we need to determine the following two time constants:</p>
  <ul>
  <li><img src=\"modelica://FM3217_2021/Resources/Images/equations/equation-Ta.png\" alt=\"T_a = J * w_n^2 / P_n\"/></li>
  <li><img src=\"modelica://FM3217_2021/Resources/Images/equations/equation-Tw.png\" alt=\"T_w = Q_n/(g*H_n) * (L/A)\"/></li>
  </ul>
  <p>We do not worry about the controller parameters since we are going 
  to use a standard PID model from the standard libary.</p>
  <p>See also the list of parameters below.</p>
  </html>"));
    end Data;
  
    model MechanicalTF
  
        Data data                        annotation(Placement(visible = true, transformation(origin={-50,90},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.TransferFunction MechanicalSystem(a = {data.Ta, 0}, b = {1})  annotation(Placement(visible = true, transformation(origin={-10,-20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Add add annotation(Placement(visible = true, transformation(origin = {30, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant one(k = 1)  annotation(Placement(visible = true, transformation(origin={-10,-60},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Gain w(k = data.wn, y(unit="rad/s"))
                                              annotation(Placement(visible = true, transformation(origin = {70, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Trapezoid powerChange(amplitude = 0.2, falling = 50, offset = -0.1, period = 200, rising = 50, startTime = 25, width = 50)  annotation(Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Gain P(k = data.Pn)  annotation(Placement(visible = true, transformation(origin={-30,46},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Division division annotation(Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.Inertia inertia(J = data.J, w(fixed = true, start = data.wn))  annotation(Placement(visible = true, transformation(origin = {58, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Torque torque annotation(Placement(visible = true, transformation(origin = {30, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(visible = true, transformation(origin = {80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
  connect(powerChange.y, MechanicalSystem.u) annotation(
        Line(points = {{-79, 0}, {-60.5, 0}, {-60.5, -30}, {-24, -30}}, color = {0, 0, 127}));
        connect(P.u, powerChange.y) annotation(Line(points={{-42,46},{-60,46},{-60,0},{-79,0}},                    color = {0, 0, 127}));
        connect(division.u1, P.y) annotation(Line(points={{-12,46},{-19,46}},                            color = {0, 0, 127}));
        connect(speedSensor.w, division.u2) annotation(Line(points={{80,19},
                {79.75,19},{79.75,10},{-20,10},{-20,34},{-12,34}},                                                                          color = {0, 0, 127}));
        connect(inertia.flange_b, speedSensor.flange) annotation(Line(points = {{68, 40}, {80, 40}, {80, 40}, {80, 40}}));
    connect(torque.flange, inertia.flange_a) annotation(Line(points = {{40, 40}, {48, 40}}));
    connect(division.y, torque.tau) annotation(Line(points={{11,40},{18,
                40}},                                                              color = {0, 0, 127}));
    connect(add.y, w.u) annotation(Line(points = {{41, -40}, {58, -40}}, color = {0, 0, 127}));
    connect(one.y, add.u2) annotation(Line(points={{1,-60},{11,-60},{11,-46},{18,-46}},          color = {0, 0, 127}));
    connect(MechanicalSystem.y, add.u1) annotation(Line(points={{1,-20},{9.5,-20},{9.5,-34},{18,-34}},          color = {0, 0, 127}));
      annotation (experiment(
          StopTime=600,
          Tolerance=1e-05));
    end MechanicalTF;
  
    model WaterWayTF
      HydroPower.HydroSystems.PipeValve pipeValve(
          endD=data.d,
          L(displayUnit="m") = data.L[1] + data.L[2],
          ZL=data.H_n,
          m_dot_nom=data.Q_n*1000,
          dp_nom=550000)
        annotation (Placement(transformation(extent={{-2,30},{18,50}})));
      HydroPower.SinksAndSources.Fixed_pT Source_pT(paraOption=false)
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
      HydroPower.SinksAndSources.Fixed_pT Sink_pT(paraOption=false)
        annotation (Placement(transformation(extent={{70,30},{50,50}})));
      Modelica.Blocks.Sources.Ramp Opening(
        offset=1,
        startTime=50,
        height=-0.9,
        duration=10)
        annotation (Placement(transformation(extent={{-100,-70},{-80,
                  -50}})));
      inner HydroPower.System_HPL system_HPL(
          Q_start=data.Q_n,                  constantTemperature=true,
        pipeRoughness=0,
        steadyState=true)
        annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      Modelica.Blocks.Continuous.TransferFunction waterWay(
          b={-data.Tw,1},
          a={data.Tw/2,1},
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          y_start=1)
        annotation (Placement(transformation(extent={{-40,-70},{-20,-50}})));
      Modelica.Blocks.Sources.RealExpression Q_valve(y=pipeValve.summary.Q_valve)
        annotation (Placement(transformation(extent={{-58,-12},{0,12}})));
      Modelica.Blocks.Sources.RealExpression dp_valve(y=pipeValve.summary.dp_valve)
        annotation (Placement(transformation(extent={{-58,-32},{0,-10}})));
      Modelica.Blocks.Math.Product Ph
        annotation (Placement(transformation(extent={{20,-20},{40,0}})));
        Modelica.Blocks.Math.Gain Ph_TF(k=data.Pn) annotation (
            Placement(transformation(extent={{20,-70},{40,-50}})));
        Data data(Pn(displayUnit="W"))
                       annotation (Placement(transformation(extent={{
                  -60,80},{-40,100}})));
        Modelica.Blocks.Math.Feedback error annotation (Placement(transformation(extent={{60,-20},{80,0}})));
    equation
      connect(Source_pT.b, pipeValve.a) annotation (Line(
          points={{-29,40},{-3,40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(Opening.y,pipeValve. ValveCtrl) annotation (Line(
          points={{-79,-60},{-60,-60},{-60,60},{8,60},{8,51}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(waterWay.u,Opening. y) annotation (Line(
          points={{-42,-60},{-79,-60}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pipeValve.b, Sink_pT.b)
        annotation (Line(points={{19,40},{49,40}},         color={0,0,255}));
        connect(waterWay.y, Ph_TF.u) annotation (Line(points={{-19,-60},{18,-60}},
                           color={0,0,127}));
        connect(Q_valve.y, Ph.u1) annotation (Line(points={{2.9,0},{10,0},{10,-4},{18,-4}}, color={0,0,127}));
        connect(dp_valve.y, Ph.u2) annotation (Line(points={{2.9,-21},{9.45,-21},{9.45,-16},{18,-16}}, color={0,0,127}));
        connect(Ph.y, error.u1) annotation (Line(points={{41,-10},{62,-10}}, color={0,0,127}));
        connect(Ph_TF.y, error.u2) annotation (Line(points={{41,-60},{70,-60},{70,-18}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=300, __Dymola_Algorithm="Radau"));
    end WaterWayTF;
  
    model TFbasedModel
  
      Modelica.Blocks.Continuous.TransferFunction mechanicalSystem(a={data.Ta,
              0})
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
      Modelica.Blocks.Continuous.TransferFunction waterWay(b={-data.Tw,
              1}, a={data.Tw/2,1})
        annotation (Placement(transformation(extent={{-20,40},{0,60}})));
      Modelica.Blocks.Math.Feedback feedback
        annotation (Placement(transformation(extent={{10,40},{30,60}})));
      HydroPower.ControllersAndSensors.TurbineGovernorAnalog turbineGovernorTF(
          enableDroop=false,
        Kd_noLoad=0.05,
        Ki_noLoad=0.025,
        K_noLoad=0.2,
        tRamp=40,
          DeadBand=0,
          P_generator_nom=data.Pn)
                    annotation (Placement(transformation(extent={{-60,40},{-40,
                60}},
              rotation=0)));
      Modelica.Blocks.Sources.Constant zero(k=0)
        annotation (Placement(transformation(extent={{-92,54},{-80,66}})));
      Modelica.Blocks.Sources.BooleanConstant MCBon(k=false)
        annotation (Placement(transformation(extent={{-82,72},{-70,84}})));
      Modelica.Blocks.Sources.RealExpression turbineLosses(y=2.26e6)
        annotation (Placement(transformation(extent={{-60,16},{-30,36}})));
      Modelica.Blocks.Sources.RealExpression gridLoad(y=0e6)
        annotation (Placement(transformation(extent={{-60,4},{-30,24}})));
        Data data(
          Pn(displayUnit="W") = 45e6,
          Q_n=92,
          H_n=50,
          d={5.5,5.5},
          L={200,100},
          J=183000,
          rpm_n=500) annotation (Placement(transformation(extent={{-60,
                  80},{-40,100}})));
        Modelica.Blocks.Math.Add add(k1=1/data.Pn, k2=1/data.Pn)
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
    equation
  
      connect(waterWay.y, feedback.u1) annotation (Line(
          points={{1,50},{12,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(mechanicalSystem.u, feedback.y) annotation (Line(
          points={{38,50},{29,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(turbineGovernorTF.y, waterWay.u) annotation (Line(
          points={{-39,50},{-22,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(mechanicalSystem.y, turbineGovernorTF.f) annotation (Line(
          points={{61,50},{72,50},{72,70},{-66,70},{-66,43},{-61,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(zero.y, turbineGovernorTF.P_generator) annotation (Line(
          points={{-79.4,60},{-70,60},{-70,57},{-61,57}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(MCBon.y, turbineGovernorTF.isMCB) annotation (Line(
          points={{-69.4,78},{-50,78},{-50,61}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(turbineGovernorTF.P_reference, zero.y) annotation (Line(
          points={{-56,61},{-56,64},{-70,64},{-70,60},{-79.4,60}},
          color={0,0,127},
          smooth=Smooth.None));
        connect(turbineLosses.y, add.u1) annotation (Line(points={{
                -28.5,26},{-28.5,26},{-12,26}}, color={0,0,127}));
        connect(gridLoad.y, add.u2) annotation (Line(points={{-28.5,14},
                {-28.5,14},{-12,14}}, color={0,0,127}));
        connect(add.y, feedback.u2) annotation (Line(points={{11,20},{
                20,20},{20,42}}, color={0,0,127}));
  
      annotation (experiment(
          StopTime=600,
          Tolerance=1e-05));
    end TFbasedModel;
  
    model CleanExample
        "Hydro plant model connecting to grid at t=150s and disconnect at t=350s"
      extends Modelon.Icons.Experiment;
  
        HydroPower.HydroSystems.Reservoir reservoir1(
          L=500,
          depthIntake={0,15},
          steadyState=false,
          H_start=fill((100), (4)),
          Hmax=fill((105), (4)),
          n=4) annotation (Placement(transformation(extent={{-60,-74},{-40,-54}}, rotation=0)));
        HydroPower.HydroSystems.Reservoir reservoir2(
          L=500,
          nSegmentIntake={1,3},
          depth={70,70,70,70},
          depthIntake={0,1},
          steadyState=false,
          H_start=fill((60), (4)),
          Hmax=fill((70), (4)),
          n=4) annotation (Placement(transformation(extent={{70,-90},{90,-70}}, rotation=0)));
      HydroPower.HydroSystems.Pipe headrace(
        L=200,
        ZL=90,
        ZR=40,
        endL={5,5},
        endD={5.5,5.5},
        Q_start=0.1,
        n=5,
        p_start=540000,
        enable_dataVizPort_lower=false) annotation (Placement(
            transformation(extent={{-30,-80},{-10,-60}}, rotation=0)));
      HydroPower.MechanicalSystems.BasicTurbine turbine(
        np=32,
        H_nom=50,
        tableOnFile=true,
        LdraftTube=10,
        DavDraftTube=2,
        LscrollCase=5,
        DavScrollCasing=2,
        PUInFlowTables=true,
        QTableName="Qtab",
        Q_nom=92,
        H_start=100,
        H_start_draftTube=40,
        Ty=0.4,
        yvLim1=[-0.1, 0.1],
        yvLim2=[-0.2, 0.2],
        TurbineDataFile=Modelica.Utilities.Files.loadResource(HydroPower.TABLE_DIR
             + "TurbineDataFile.mat"),
        P_nom=45000000) annotation (Placement(transformation(extent={{-2,-80},
                  {18,-60}},
                     rotation=0)));
  
      HydroPower.HydroSystems.Pipe tailrace(
        n=4,
        endL={5,5},
        ZL=40,
        ZR=0,
        L=100,
        endD={5.5,5.5},
        Q_start=0.1,
        p_start=400000,
        enable_dataVizPort_lower=false) annotation (Placement(
            transformation(extent={{30,-80},{50,-60}}, rotation=0)));
      HydroPower.ElectricalSystems.PowerGrid powerGrid(
        startTime=1e6,
        unitsJ={122000,5.5e6,8000},
        NoLoadUnits={200,400,1000},
        distNoGen={-2,0,0},
        distTgen={150,1e6,1e6}) annotation (Placement(transformation(
              extent={{-80,-40},{-60,-20}}, rotation=0)));
  
      Modelica.Blocks.Sources.Ramp pwr_ref(
        duration=10,
        height=0,
        offset=45e6,
        startTime=1e6) annotation (Placement(transformation(extent={{-25,-16},
                {-13,-4}},
                      rotation=0)));
      HydroPower.ElectricalSystems.GeneratorAndMCB generator(
        Kdmp={0.05},
        f_start=0,
        timeMCB_open={200},
          J={183000},
          np={12},
          P_nom={45000000},
          timeMCB_close={1e6})
                          annotation (Placement(transformation(extent={{-50,-40},
                {-30,-20}},
                       rotation=0)));
      HydroPower.ControllersAndSensors.TurbineGovernorAnalog turbineGovernor(
        ep=1,
        DeadBand=0.001,
        Ki_load=0.1,
        Kd_load=0.5,
        Kd_noLoad=0.05,
        Ki_noLoad=0.025,
        K_noLoad=0.2,
        K_load=0.4,
        tRamp=40,
        P_generator_nom=generator.P_nom[1],
        enableRamp=false) annotation (Placement(transformation(extent=
               {{-10,-40},{10,-20}}, rotation=0)));
      inner HydroPower.System_HPL system_HPL(
        steadyState=true,
        pipeRoughness=0.1,
        T_start=293) annotation (Placement(transformation(extent={{
                -100,-80},{-80,-60}})));
    equation
  
      connect(powerGrid.f_grid, generator.f_grid) annotation (Line(
          points={{-59,-37},{-51,-37}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(powerGrid.P_grid_balance, generator.P_grid_balance) annotation (Line(
          points={{-59,-23},{-51,-23}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pwr_ref.y, turbineGovernor.P_reference) annotation (Line(
          points={{-12.4,-10},{-6,-10},{-6,-19}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.f_out[1], turbineGovernor.f) annotation (Line(
          points={{-29,-37},{-11,-37}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.onMCB, powerGrid.MCB) annotation (Line(
          points={{-40,-19},{-40,0},{-70,0},{-70,-19}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.onMCB[1], turbineGovernor.isMCB) annotation (Line(
          points={{-40,-19},{-40,0},{0,0},{0,-19}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.P_out[1], turbineGovernor.P_generator) annotation (Line(
          points={{-29,-23},{-11,-23}},
          color={0,0,127},
          smooth=Smooth.None));
  
      connect(generator.f_out, powerGrid.f) annotation (Line(
          points={{-29,-37},{-20,-37},{-20,-52},{-90,-52},{-90,-37},{
              -81,-37}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(reservoir1.a2_pipe, headrace.a) annotation (Line(
          points={{-39,-70},{-31,-70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(turbine.a, headrace.b) annotation (Line(
          points={{-3,-70},{-9,-70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(turbine.b, tailrace.a) annotation (Line(
          points={{19,-70},{29,-70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(tailrace.b, reservoir2.a1_pipe) annotation (Line(
          points={{51,-70},{60,-70},{60,-86},{69,-86}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(powerGrid.J_grid, generator.J_grid)
        annotation (Line(points={{-59,-30},{-59,-30},{-51,-30}},
                                                              color={0,0,127}));
      connect(turbineGovernor.y, turbine.yGV)
        annotation (Line(points={{11,-30},{12,-30},{16,-30},{16,-60},{
                16,-59},{14,-59}},                      color={0,0,127}));
      connect(generator.f_out[1], turbine.f_generator) annotation (Line(points={{-29,-37},
                {-29,-37},{-20,-37},{-20,-38},{-20,-36},{-20,-52},{2,
                -52},{2,-59}},                           color={0,0,127}));
      connect(generator.P_turbine[1], turbine.TurbineData[1])
        annotation (Line(points={{-48,-41},{-48,-48},{8,-48},{8,-59.6667}},
                          color={0,0,127}));
      annotation (
         experiment(
          StopTime=600,
          Tolerance=1e-005,
          __Dymola_Algorithm="Radau"));
    end CleanExample;
  
    model CombinedModel
        extends TFbasedModel;
        extends CleanExample;
        annotation (experiment(
            StopTime=600,
            Tolerance=1e-05,
            __Dymola_Algorithm="Radau"));
    end CombinedModel;
  
      package OpenModelica
        extends Modelica.Icons.Package;
        record Data
          extends Modelica.Icons.Record;
          extends OpenHPL.Data;
          parameter Modelica.SIunits.Power Pn = 45e6 "Nominal Power";
          parameter Modelica.SIunits.Frequency fn = 50 "Nominal frequency";
          parameter Integer p = 12 "Number of poles of the generator";
          parameter Modelica.SIunits.AngularVelocity wn = 2*Modelica.Constants.pi * fn /(p/2) "Nominal speed";
          parameter Modelica.SIunits.Inertia J = 183e3 "Inertia of gen and turbine";
          parameter Modelica.SIunits.VolumeFlowRate Qn = 92 "Nominal flow rate";
          parameter Modelica.SIunits.Height Hn = 50 "Nominal head";
          parameter Modelica.SIunits.Length L[2] = {200,100} "Lenght of the water way";
          parameter Modelica.SIunits.Length d[2] = {5.5,5.5} "Pipe diameter";
          parameter Modelica.SIunits.Area A[2] = Modelica.Constants.pi/4*d.^ 2 "Area of pipe segments";
          parameter Modelica.SIunits.Time Ta = J * wn^2/Pn "Mechanical time constant";
          parameter Modelica.SIunits.Time Tw = Qn/ (Modelica.Constants.g_n *Hn) * (L[1]/A[1] + L[2]/A[2]) "Waterway time constant";
  annotation(
          Documentation(info = "<html>
  <p>
  The transfer function of the hydro power system can be represented by:
  </p>
  <img src=\"modelica://FM3217_2021/Resources/Images/TFsystem.png\">
  <p>For the determining the time constants of the transfer-function based model we need to determine the following two time constants:</p>
  <ul>
  <li><img src=\"modelica://FM3217_2021/Resources/Images/equations/equation-Ta.png\" alt=\"T_a = J * w_n^2 / P_n\"/></li>
  <li><img src=\"modelica://FM3217_2021/Resources/Images/equations/equation-Tw.png\" alt=\"T_w = Q_n/(g*H_n) * (L/A)\"/></li>
  </ul>
  <p>We do not worry about the controller parameters since we are going 
  to use a standard PID model from the standard libary.</p>
</html>"));
        end Data;
  
        model MechanicalTF
          Modelica.Blocks.Continuous.TransferFunction mechanical(a={data.Ta,0}) annotation (
            Placement(visible = true, transformation(origin = {0, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Sources.Trapezoid powerChange(
            amplitude=0.2,
            falling=50,
            offset=-0.1,
            period=200,
            rising=50,
            startTime=25,
            width=50) annotation (
            Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Sources.Constant one(k=1) annotation (
            Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Add add annotation (
            Placement(visible = true, transformation(origin = {40, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Gain gain(k=data.wn) annotation (
            Placement(visible = true, transformation(origin = {70, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Division division annotation (
            Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Gain P(k=data.Pn) annotation (
            Placement(visible = true, transformation(origin = {-30, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (
            Placement(visible = true, transformation(origin = {80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
          Modelica.Mechanics.Rotational.Sources.Torque torque annotation (
            Placement(visible = true, transformation(origin = {30, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Mechanics.Rotational.Components.Inertia inertia(J=data.J, w(fixed=true, start=data.wn)) annotation (
            Placement(visible = true, transformation(origin = {58, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FM3217_2021.Tutorial11.OpenModelica.Data data annotation(
          Placement(visible = true, transformation(origin = {-50, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
          connect(powerChange.y, mechanical.u) annotation (
            Line(points = {{-78, 0}, {-60, 0}, {-60, -52}, {-12, -52}, {-12, -52}}, color = {0, 0, 127}));
          connect(add.u1, mechanical.y) annotation (
            Line(points = {{28, -62}, {20, -62}, {20, -52}, {12, -52}}, color = {0, 0, 127}));
          connect(one.y, add.u2) annotation (
            Line(points = {{12, -80}, {20, -80}, {20, -74}, {28, -74}}, color = {0, 0, 127}));
          connect(gain.u, add.y) annotation (
            Line(points = {{58, -68}, {52, -68}, {52, -68}, {52, -68}}, color = {0, 0, 127}));
          connect(division.y, torque.tau) annotation (
            Line(points = {{11, 40}, {18, 40}}, color = {0, 0, 127}));
          connect(speedSensor.w, division.u2) annotation (
            Line(points = {{80, 19}, {79.75, 19}, {79.75, 10}, {-20, 10}, {-20, 34}, {-12, 34}}, color = {0, 0, 127}));
          connect(torque.flange, inertia.flange_a) annotation (
            Line(points = {{40, 40}, {48, 40}}));
          connect(division.u1, P.y) annotation (
            Line(points = {{-12, 46}, {-19, 46}}, color = {0, 0, 127}));
          connect(inertia.flange_b, speedSensor.flange) annotation (
            Line(points = {{68, 40}, {80, 40}, {80, 40}, {80, 40}}));
          connect(
            P.u, powerChange.y) annotation (
            Line(points = {{-42, 46}, {-60, 46}, {-60, 0}, {-78, 0}, {-78, 0}}, color = {0, 0, 127}));
          annotation (
            Icon(coordinateSystem(grid = {2, 0})),
          experiment(StartTime = 0, StopTime = 600, Tolerance = 1e-6, Interval = 1.2));
        end MechanicalTF;
  
        model WaterWayTF
          Modelica.Blocks.Math.Gain Ph_TF(k = data.Pn) annotation (
            Placement(visible = true, transformation(extent = {{20, -70}, {40, -50}}, rotation = 0)));
          Modelica.Blocks.Continuous.TransferFunction waterWay(a = {data.Tw / 2, 1}, b = {-data.Tw, 1}, initType = Modelica.Blocks.Types.Init.InitialOutput, y_start = 1) annotation (
            Placement(visible = true, transformation(extent = {{-40, -70}, {-20, -50}}, rotation = 0)));
          Modelica.Blocks.Sources.Ramp Opening(duration = 10, height = -0.9, offset = 1, startTime = 50) annotation (
            Placement(visible = true, transformation(extent = {{-90, -70}, {-70, -50}}, rotation = 0)));
          OpenHPL.Waterway.Reservoir reservoir(H_r=10) annotation (
            Placement(visible = true, transformation(origin = {-70, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          OpenHPL.Waterway.Reservoir reservoir1(H_r=10) annotation (
            Placement(visible = true, transformation(origin = {70, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          OpenHPL.Waterway.Pipe pipe(
            D_i=data.d[1],
            H=data.Hn - 10,
            L=data.L[1]) annotation (
            Placement(visible = true, transformation(origin = {-30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          OpenHPL.Waterway.Pipe pipe1(
            D_i=data.d[2],
            H=10,
            L=data.L[2]) annotation (
            Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          OpenHPL.ElectroMech.Turbines.Turbine2 turbine2(
            ConstEfficiency=true,
            H_n=data.Hn,
            Jt=data.J,
            ValveCapacity=false,
            Vdot_n=data.Qn, enableP_out = true,
            eta_h=1,
            p=data.p) annotation (
            Placement(visible = true, transformation(origin = {0, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          inner FM3217_2021.Tutorial11.OpenModelica.Data data(
            Steady=true,
            V_0=data.Qn,
            gamma_air=1) annotation (
            Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
  connect(waterWay.u, Opening.y) annotation(
          Line(points = {{-42, -60}, {-69, -60}}, color = {0, 0, 127}));
          connect(waterWay.y, Ph_TF.u) annotation (
            Line(points = {{-19, -60}, {18, -60}}, color = {0, 0, 127}));
  connect(pipe.i, reservoir.o) annotation(
          Line(points = {{-40, 30}, {-49, 30}, {-49, 40}, {-60, 40}}, color = {28, 108, 200}));
  connect(pipe1.o, reservoir1.o) annotation(
          Line(points = {{40, 30}, {50, 30}, {50, 20}, {60, 20}}, color = {28, 108, 200}));
  connect(pipe.o, turbine2.i) annotation(
          Line(points = {{-20, 30}, {-10, 30}}, color = {28, 108, 200}));
  connect(turbine2.o, pipe1.i) annotation(
          Line(points = {{10, 30}, {20, 30}}, color = {28, 108, 200}));
  connect(turbine2.u_t, Opening.y) annotation(
          Line(points = {{-8, 42}, {-8, 70}, {-90, 70}, {-90, -20}, {-60, -20}, {-60, -60}, {-69, -60}}, color = {0, 0, 127}));
          annotation (
            Icon(coordinateSystem(grid = {2, 0})),
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
        end WaterWayTF;
  
        model TFbasedModel
  
          Modelica.Blocks.Continuous.TransferFunction mechanicalSystem(a={data.Ta,
                  0})
            annotation (Placement(visible = true, transformation(extent = {{50, -50}, {70, -30}}, rotation = 0)));
          Modelica.Blocks.Continuous.TransferFunction waterWay(b={-data.Tw,
                  1}, a={data.Tw/2,1})
            annotation (Placement(visible = true, transformation(extent = {{-10, -50}, {10, -30}}, rotation = 0)));
          Modelica.Blocks.Math.Feedback feedback
            annotation (Placement(visible = true, transformation(extent = {{20, -50}, {40, -30}}, rotation = 0)));
          Modelica.Blocks.Sources.RealExpression turbineLosses(y=2.26e6)
            annotation (Placement(visible = true, transformation(extent = {{-50, -74}, {-20, -54}}, rotation = 0)));
            Modelica.Blocks.Math.Add add(k1=1/data.Pn, k2=1/data.Pn)
              annotation (Placement(visible = true, transformation(extent = {{0, -80}, {20, -60}}, rotation = 0)));
        Modelica.Blocks.Continuous.PI piTF(T=10) annotation (Placement(visible=true, transformation(
                  origin={-40,-40},
                  extent={{-10,-10},{10,10}},
                  rotation=0)));
            Modelica.Blocks.Math.Gain gain(k=-1) annotation (Placement(visible = true, transformation(extent = {{-72, -44}, {-64, -36}}, rotation = 0)));
        Modelica.Blocks.Sources.Pulse gridLoad(amplitude = 10e6, period = 400, startTime = 200)  annotation (
              Placement(visible = true, transformation(origin = {-30, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          inner FM3217_2021.Tutorial11.OpenModelica.Data data(Steady=true, V_0=data.Qn) annotation (
            Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
            connect(waterWay.y, feedback.u1) annotation (
              Line(points = {{11, -40}, {22, -40}}, color = {0, 0, 127}));
            connect(mechanicalSystem.u, feedback.y) annotation (
              Line(points = {{48, -40}, {39, -40}}, color = {0, 0, 127}));
            connect(turbineLosses.y, add.u1) annotation (
              Line(points = {{-18.5, -64}, {-18.5, -64}, {-2, -64}}, color = {0, 0, 127}));
            connect(add.y, feedback.u2) annotation (
              Line(points = {{21, -70}, {30, -70}, {30, -48}}, color = {0, 0, 127}));
            connect(piTF.y, waterWay.u) annotation (
              Line(points = {{-29, -40}, {-12, -40}, {-12, -40}, {-12, -40}}, color = {0, 0, 127}));
            connect(mechanicalSystem.y, gain.u) annotation (
              Line(points = {{71, -40}, {80, -40}, {80, -20}, {-80, -20}, {-80, -40}, {-72.8, -40}}, color = {0, 0, 127}));
            connect(gain.y, piTF.u) annotation (
              Line(points = {{-63.6, -40}, {-52, -40}}, color = {0, 0, 127}));
        connect(gridLoad.y, add.u2) annotation (
              Line(points={{-19,-86},{-14,-86},{-14,-76},{-2,-76}},          color = {0, 0, 127}));
            annotation (experiment(
              StopTime=600,
              Tolerance=1e-05, StartTime = 0, Interval = 1.2));
        end TFbasedModel;
  
        model CombinedModel
            extends TFbasedModel(data(Steady = false), mechanicalSystem.initType = Modelica.Blocks.Types.Init.SteadyState, mechanicalSystem.y_start = 0);
        OpenHPL.Waterway.Pipe pipe1(D_i = data.d[2], H = 10, L = data.L[2]) annotation (
              Placement(visible = true, transformation(origin={30,20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OpenHPL.ElectroMech.Turbines.Turbine2 turbine2(
            ConstEfficiency=true,
            H_n=data.Hn, Jt = data.J / 2, Ploss = 1.13e6,
            ValveCapacity=false,
            Vdot_n=92, eta_h = 1, p = data.p) annotation (Placement(visible=true, transformation(
                origin={0,20},
                extent={{-10,-10},{10,10}},
                rotation=0)));
        OpenHPL.Waterway.Reservoir reservoir1(H_r = 10) annotation (
              Placement(visible = true, transformation(origin={62,20},    extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        OpenHPL.Waterway.Reservoir reservoir(H_r = 10) annotation (
              Placement(visible = true, transformation(origin={-70,20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OpenHPL.Waterway.Pipe pipe(D_i = data.d[1], H = data.Hn - 10, L = data.L[1], vertical = false) annotation (
              Placement(visible = true, transformation(origin={-30,20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OpenHPL.ElectroMech.Generators.SimpleGen2 simpleGen2(Jg = data.J / 2, Ploss = 1.13e6,
            enable_f=true, enable_w = false)                                                                                          annotation (
              Placement(visible = true, transformation(origin={0,50},    extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.PI piReal(T = 10, k = 1) annotation (Placement(visible=true, transformation(
                  origin={-30,68},
                  extent={{-10,-10},{10,10}},
                  rotation=0)));
        Modelica.Blocks.Math.Feedback feedback2 annotation (
              Placement(visible = true, transformation(extent={{-74,58},{-54,78}},      rotation = 0)));
        Modelica.Blocks.Sources.Constant sysF(k = 50)  annotation (
              Placement(visible = true, transformation(origin={-92,68},    extent = {{-6, -6}, {6, 6}}, rotation = 0)));
        Modelica.Blocks.Math.Gain df_pu(k = 1 / 50) annotation (
              Placement(visible = true, transformation(extent={{-36, 86},{-28, 94}},      rotation = 0)));
        equation
          connect(reservoir1.o, pipe1.o) annotation (
            Line(points = {{52, 20}, {40, 20}}, color = {28, 108, 200}));
          connect(reservoir.o, pipe.i) annotation (
            Line(points = {{-60, 20}, {-40, 20}}, color = {28, 108, 200}));
          connect(pipe.o, turbine2.i) annotation (
            Line(points = {{-20, 20}, {-10, 20}}, color = {28, 108, 200}));
          connect(turbine2.o, pipe1.i) annotation (
            Line(points = {{10, 20}, {20, 20}}, color = {28, 108, 200}));
          connect(sysF.y, feedback2.u1) annotation (
            Line(points = {{-85, 68}, {-72, 68}}, color = {0, 0, 127}));
          connect(simpleGen2.f, feedback2.u2) annotation (
            Line(points = {{-11, 46}, {-64, 46}, {-64, 60}}, color = {0, 0, 127}));
          connect(feedback2.y, piReal.u) annotation (
            Line(points = {{-55, 68}, {-42, 68}}, color = {0, 0, 127}));
          connect(piReal.y, turbine2.u_t) annotation (
            Line(points = {{-19, 68}, {-8, 68}, {-8, 32}}, color = {0, 0, 127}));
          connect(df_pu.u, feedback2.y) annotation (
            Line(points = {{-37, 90}, {-46, 90}, {-46, 68}, {-55, 68}}, color = {0, 0, 127}));
          connect(turbine2.flange, simpleGen2.flange) annotation (
            Line(points = {{0, 29.8}, {0, 40.2}}, color = {0, 0, 0}));
          connect(
            simpleGen2.Pload, gridLoad.y) annotation (
            Line(points = {{0, 62}, {0, 62}, {0, 80}, {90, 80}, {90, -86}, {-18, -86}, {-18, -86}}, color = {0, 0, 127}));
          annotation (experiment(
                StopTime = 600, StartTime = 0, Tolerance = 1e-06, Interval = 1.2));
        end CombinedModel;
        annotation (
          Icon(coordinateSystem(grid = {2, 0})));
      end OpenModelica;
  
  end Tutorial11;
  annotation (uses(Modelica(version="3.2.3"), HydroPower(version="2.13"),
      Modelon(version="3.7")));
end FM3217_2021;
=======
within ;
package FM3217_2021 "Course files of the 2021 tutorial"

  package Tutorial1 "Simple Pendulum"

    model SimplePendulum "First simple version"
      constant Real g(unit="m/s2")=9.81 "Gravitationl constant";
      parameter Real L(unit="m")=1 "Length of the pendulum";
      Real Theta(unit="rad",start=0.1) "Angle of the pendulum";
      Real ThetaDot "Helping variable for 2nd derivative";
    equation
       ThetaDot = der(Theta);
       der(ThetaDot) = - g/L * sin(Theta);
      annotation (experiment(StopTime=10));
    end SimplePendulum;

    model SimplePendulumSIunits "Model using SI units"
      constant Modelica.SIunits.Acceleration g=9.81 "Gravitationl constant";
      parameter Modelica.SIunits.Length L=1 "Length of the pendulum";
      Modelica.SIunits.Angle Theta(start=0.1) "Angle of the pendulum";
      Real ThetaDot  "Helping variable for 2nd derivative";
    equation
       ThetaDot = der(Theta);
       der(ThetaDot) = - g/L * sin(Theta);
      annotation (experiment(StopTime=10));
    end SimplePendulumSIunits;

    model SimplePendulumUsingImports "Model using import"
      import SI = Modelica.SIunits;

      constant SI.Acceleration g=9.81 "Gravitationl constant";
      parameter SI.Length L=1 "Length of the pendulum";
      SI.Angle Theta(start=0.1) "Angle of the pendulum";
      Real ThetaDot  "Helping variable for 2nd derivative";
    equation
       ThetaDot = der(Theta);
       der(ThetaDot) = - g/L * sin(Theta);
      annotation (experiment(StopTime=10));
    end SimplePendulumUsingImports;
  end Tutorial1;

  package Tutorial2 "Motor and Motor Drive"

    model Motor

      parameter Modelica.SIunits.Resistance Ra = 0.5 "Resistance of the armature";


      Modelica.Electrical.Analog.Basic.Resistor resistor(R=Ra)  annotation (Placement(transformation(extent={{-56,10},{-36,30}})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-80,-62},{-60,-42}})));
      Modelica.Electrical.Analog.Basic.Inductor inductor(L=La)   annotation (Placement(transformation(extent={{-20,10},{0,30}})));
      Modelica.Electrical.Analog.Basic.EMF emf annotation (Placement(transformation(extent={{10,-10},{30,10}})));
      Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=-90,
            origin={-70,0})));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J=Ja)    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
      Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
      Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_b1 "Flange of right shaft"
                                                annotation (Placement(transformation(extent={{90,-10},{110,10}})));
      parameter Modelica.SIunits.Inductance La=0.05 "Inductance of the armature";
      parameter Modelica.SIunits.Inertia Ja=0.05 "Inertia of the motor";
    equation
      connect(signalVoltage.p, resistor.p) annotation (Line(points={{-70,10},{-70,20},{-56,20}}, color={0,0,255}));
      connect(resistor.n, inductor.p) annotation (Line(points={{-36,20},{-20,20}}, color={0,0,255}));
      connect(inductor.n, emf.p) annotation (Line(points={{0,20},{20,20},{20,10}}, color={0,0,255}));
      connect(emf.n, signalVoltage.n) annotation (Line(points={{20,-10},{20,-28},{-70,-28},{-70,-10}}, color={0,0,255}));
      connect(emf.n, ground.p) annotation (Line(points={{20,-10},{20,-28},{-70,-28},{-70,-42}}, color={0,0,255}));
      connect(emf.flange, inertia.flange_a) annotation (Line(points={{30,0},{40,0}}, color={0,0,0}));
      connect(signalVoltage.v, u) annotation (Line(points={{-82,0},{-120,0}}, color={0,0,127}));
      connect(inertia.flange_b, flange_b1) annotation (Line(points={{60,0},{100,0}}, color={0,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Bitmap(extent={{-98,-76},{102,128}}, fileName="modelica://FM3217_2021/Resources/Images/dc-motor.jpg")}),
                                                                     Diagram(coordinateSystem(preserveAspectRatio=false)));
    end Motor;

    model MotorDrive
      Motor motor(Ja=2) annotation (Placement(transformation(extent={{-2,-10},{18,10}})));
      Modelica.Blocks.Sources.Step step(startTime=0.1) annotation (Placement(transformation(extent={{-100,-10},{-80,10}})));
      Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
      Modelica.Blocks.Continuous.PID PID(Ti=0.1, Td=0) annotation (Placement(transformation(extent={{-38,-10},{-18,10}})));
      Modelica.Mechanics.Rotational.Components.IdealGear idealGear(ratio=100) annotation (Placement(transformation(extent={{34,-10},{54,10}})));
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J=5) annotation (Placement(transformation(extent={{70,-10},{90,10}})));
      Modelica.Mechanics.Rotational.Sensors.AngleSensor angleSensor annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={90,-30})));
    equation
      connect(step.y, feedback.u1) annotation (Line(points={{-79,0},{-68,0}}, color={0,0,127}));
      connect(feedback.y, PID.u) annotation (Line(points={{-51,0},{-40,0}}, color={0,0,127}));
      connect(idealGear.flange_a, motor.flange_b1) annotation (Line(points={{34,0},{18,0}}, color={0,0,0}));
      connect(angleSensor.flange, inertia.flange_b) annotation (Line(points={{90,-20},{90,0}}, color={0,0,0}));
      connect(angleSensor.phi, feedback.u2) annotation (Line(points={{90,-41},{90,-50},{-60,-50},{-60,-8}}, color={0,0,127}));
      connect(PID.y, motor.u) annotation (Line(points={{-17,0},{-4,0}}, color={0,0,127}));
      connect(idealGear.flange_b, inertia.flange_a) annotation (Line(points={{54,0},{70,0}}, color={0,0,0}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{120,100}})),
        Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})),
        experiment(StopTime=100));
    end MotorDrive;
  end Tutorial2;

  package Tutorial3 "Simple machine with loads"
    package Components
      model Machine

         parameter Modelica.SIunits.Inductance La=0.05 "Inductance of the armature" annotation (Dialog(group="Electrical"));
        parameter Modelica.SIunits.Inertia Ja=0.05 "Inertia of the motor" annotation (Dialog(group="Mechanical"));

         parameter Modelica.SIunits.Resistance Ra=0.5 "Resistance of the armature" annotation (Dialog(group="Electrical"));

        Modelica.Electrical.Analog.Basic.Resistor resistor(R=Ra)  annotation (Placement(transformation(extent={{-56,10},{-36,30}})));
        Modelica.Electrical.Analog.Basic.Inductor inductor(L=La)   annotation (Placement(transformation(extent={{-20,10},{0,30}})));
        Modelica.Electrical.Analog.Basic.EMF emf annotation (Placement(transformation(extent={{10,-10},{30,10}})));
        Modelica.Mechanics.Rotational.Components.Inertia inertia(J=Ja)    annotation (Placement(transformation(extent={{40,-10},{60,10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_b flange "Flange of right shaft" annotation (Placement(transformation(extent={{90,-10},{110,10}})));
        Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive electrical pin" annotation (Placement(transformation(extent={{-110,70},{-90,90}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative electrical pin" annotation (Placement(transformation(extent={{-110,-90},{-90,-70}})));
      equation
        connect(resistor.n, inductor.p) annotation (Line(points={{-36,20},{-20,20}}, color={0,0,255}));
        connect(inductor.n, emf.p) annotation (Line(points={{0,20},{20,20},{20,10}}, color={0,0,255}));
        connect(emf.flange, inertia.flange_a) annotation (Line(points={{30,0},{40,0}}, color={0,0,0}));
        connect(inertia.flange_b, flange) annotation (Line(points={{60,0},{100,0}}, color={0,0,0}));
        connect(resistor.p, p) annotation (Line(points={{-56,20},{-80,20},{-80,80},{-100,80}}, color={0,0,255}));
        connect(emf.n, n) annotation (Line(points={{20,-10},{20,-80},{-100,-80}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Bitmap(extent={{-98,-76},{102,128}}, fileName="modelica://FM3217_2021/Resources/Images/dc-motor.jpg")}),
                                                                       Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Machine;

      model Turbine
        Modelica.Mechanics.Rotational.Sources.ConstantTorque constantTorque(tau_constant=Tt) annotation (Placement(transformation(extent={{76,-10},{56,10}})));
        Modelica.Mechanics.Rotational.Components.Inertia turbineWheel(J=Jt) annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (Placement(transformation(rotation=0, extent={{-110,-10},{-90,10}})));
        parameter Modelica.SIunits.Inertia Jt=1 "Inertia of the turbine runner wheel";
        parameter Modelica.SIunits.Torque Tt=10 "Turbine torque";
      equation
        connect(constantTorque.flange, turbineWheel.flange_b) annotation (Line(points={{56,0},{10,0}}, color={0,0,0}));
        connect(flange_a, turbineWheel.flange_a) annotation (Line(points={{-100,0},{-10,0}}, color={0,0,0}));
        annotation (Icon(graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://FM3217_2021/Resources/Images/Turbine.png")}));
      end Turbine;

      model Rload
        Modelica.Electrical.Analog.Basic.Resistor resistorUpper(R=Rl/2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,18})));
        Modelica.Electrical.Analog.Interfaces.PositivePin p "Positive electrical pin" annotation (Placement(transformation(extent={{-10,90},{10,110}})));
        Modelica.Electrical.Analog.Interfaces.NegativePin n "Negative electrical pin" annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
        parameter Modelica.SIunits.Resistance Rl=5 "Ohmic part of the load";
        Modelica.Electrical.Analog.Basic.Resistor resistorLower(R=Rl/2) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={0,-16})));
      equation
        connect(resistorUpper.p, p) annotation (Line(points={{0,28},{0,100}}, color={0,0,255}));
        connect(resistorUpper.n, resistorLower.p) annotation (Line(points={{0,8},{0,4},{1.77636e-15,4},{1.77636e-15,-6}}, color={0,0,255}));
        connect(resistorLower.n, n) annotation (Line(points={{-1.77636e-15,-26},{-1.77636e-15,-60},{0,-60},{0,-100}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
      end Rload;

      model RLload
        extends Rload;
        Modelica.Electrical.Analog.Basic.Inductor inductor(L=Ll) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={40,0})));
        parameter Modelica.SIunits.Inductance Ll=5e-3 "Inductive part of the load";
      equation
        connect(inductor.p, p) annotation (Line(points={{40,10},{40,60},{0,60},{0,100}}, color={0,0,255}));
        connect(inductor.n, n) annotation (Line(points={{40,-10},{40,-60},{0,-60},{0,-100}}, color={0,0,255}));
      end RLload;

      model RLCload
        extends RLload;
        Modelica.Electrical.Analog.Basic.Capacitor capacitor(C=Cl) annotation (Placement(transformation(
              extent={{-10,-10},{10,10}},
              rotation=270,
              origin={-40,0})));
        parameter Modelica.SIunits.Capacitance Cl=2 "Capacitive part of the load";
      equation
        connect(capacitor.p, p) annotation (Line(points={{-40,10},{-40,60},{0,60},{0,100}}, color={0,0,255}));
        connect(capacitor.n, n) annotation (Line(points={{-40,-10},{-40,-60},{0,-60},{0,-100}}, color={0,0,255}));
      end RLCload;
    end Components;

    package Tests
      model TestMachine
        extends Modelica.Icons.Example;


        Components.Machine machine annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=-90,
              origin={-30,0})));
        Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-40,-48},{-20,-28}})));
        Modelica.Blocks.Sources.Constant const(k=2)
                                               annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{-24,6},{-4,26}})));
      equation
        connect(signalVoltage.n, machine.n) annotation (Line(points={{-30,-10},{-30,-20},{0,-20},{0,-8}}, color={0,0,255}));
        connect(ground.p, signalVoltage.n) annotation (Line(points={{-30,-28},{-30,-10}}, color={0,0,255}));
        connect(const.y, signalVoltage.v) annotation (Line(points={{-49,0},{-42,0}}, color={0,0,127}));
        connect(signalVoltage.p, powerSensor.pc) annotation (Line(points={{-30,10},{-30,16},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nc, machine.p) annotation (Line(points={{-4,16},{0,16},{0,8}}, color={0,0,255}));
        connect(powerSensor.pv, powerSensor.pc) annotation (Line(points={{-14,26},{-24,26},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nv, signalVoltage.n) annotation (Line(points={{-14,6},{-14,-10},{-30,-10}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1.1));
      end TestMachine;

      model TestTurbine
        extends Modelica.Icons.Example;

        Components.Machine machine annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (Placement(transformation(
              extent={{-10,10},{10,-10}},
              rotation=-90,
              origin={-30,0})));
        Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-40,-48},{-20,-28}})));
        Modelica.Blocks.Sources.Constant const(k=2)
                                               annotation (Placement(transformation(extent={{-70,-10},{-50,10}})));
        Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{-24,6},{-4,26}})));
        Components.Turbine turbine annotation (Placement(transformation(rotation=0, extent={{56,-10},{76,10}})));
      equation
        connect(signalVoltage.n, machine.n) annotation (Line(points={{-30,-10},{-30,-20},{0,-20},{0,-8}}, color={0,0,255}));
        connect(ground.p, signalVoltage.n) annotation (Line(points={{-30,-28},{-30,-10}}, color={0,0,255}));
        connect(const.y, signalVoltage.v) annotation (Line(points={{-49,0},{-42,0}}, color={0,0,127}));
        connect(signalVoltage.p, powerSensor.pc) annotation (Line(points={{-30,10},{-30,16},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nc, machine.p) annotation (Line(points={{-4,16},{0,16},{0,8}}, color={0,0,255}));
        connect(powerSensor.pv, powerSensor.pc) annotation (Line(points={{-14,26},{-24,26},{-24,16}}, color={0,0,255}));
        connect(powerSensor.nv, signalVoltage.n) annotation (Line(points={{-14,6},{-14,-10},{-30,-10}}, color={0,0,255}));
        connect(machine.flange, turbine.flange_a) annotation (Line(points={{20,0},{56,0}}, color={0,0,0}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1.1));
      end TestTurbine;

      model TestLoads
        extends Modelica.Icons.Example;

        Components.Machine machine annotation (Placement(transformation(extent={{0,-10},{20,10}})));
        Components.Turbine turbine annotation (Placement(transformation(rotation=0, extent={{56,-10},{76,10}})));
        Components.Rload rload annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
        Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-30,-50},{-10,-30}})));
        Components.RLload rLload annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
        Components.RLCload rLCload annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      equation
        connect(machine.flange, turbine.flange_a) annotation (Line(points={{20,0},{56,0}}, color={0,0,0}));
        connect(rload.p, machine.p) annotation (Line(points={{-20,10},{-20,20},{0,20},{0,8}}, color={0,0,255}));
        connect(rload.n, machine.n) annotation (Line(points={{-20,-10},{-20,-20},{0,-20},{0,-8}}, color={0,0,255}));
        connect(rload.n, ground.p) annotation (Line(points={{-20,-10},{-20,-30}}, color={0,0,255}));
        connect(rLload.p, machine.p) annotation (Line(points={{-50,10},{-50,20},{0,20},{0,8}}, color={0,0,255}));
        connect(rLload.n, ground.p) annotation (Line(points={{-50,-10},{-50,-20},{-20,-20},{-20,-30}}, color={0,0,255}));
        connect(rLCload.p, machine.p) annotation (Line(points={{-80,10},{-80,20},{0,20},{0,8}}, color={0,0,255}));
        connect(rLCload.n, ground.p) annotation (Line(points={{-80,-10},{-80,-20},{-20,-20},{-20,-30}}, color={0,0,255}));
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
          experiment(StopTime=1.1));
      end TestLoads;

      model TestLoads20Nm
        extends TestLoads(turbine(Tt=20));
      end TestLoads20Nm;

      model TestLoads25Nm
        extends TestLoads(turbine(Tt=25), rLload(Rl=10, Ll=4e-3));
      end TestLoads25Nm;
    end Tests;
  end Tutorial3;

  package Tutorial4 "Electric kettle, coffee and milk"

    model ElectricKettle
      Modelica.Electrical.Analog.Basic.Resistor resistor(R=230^2/2000, useHeatPort=true) annotation (Placement(transformation(
            extent={{-10,10},{10,-10}},
            rotation=270,
            origin={-20,0})));
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-90,-60},{-70,-40}})));
      Modelica.Electrical.Analog.Sources.SineVoltage sineVoltage(V=230*sqrt(2), freqHz=50) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-90,0})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{-50,10},{-30,30}})));
      Modelica.Blocks.Math.Mean mean(f=50) annotation (Placement(transformation(
            extent={{-6,-6},{6,6}},
            rotation=270,
            origin={-50,-6})));
      Modelica.Thermal.HeatTransfer.Components.HeatCapacitor water(C=4.18*1000*1.7, T(
          start=283.15,
          fixed=true,
          displayUnit="degC")) annotation (Placement(transformation(extent={{0,12},{20,32}})));
      Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor temperatureSensor annotation (Placement(transformation(extent={{22,-10},{42,10}})));
      Modelica.Electrical.Analog.Ideal.IdealOpeningSwitch switch annotation (Placement(transformation(extent={{-80,10},{-60,30}})));
      Modelica.Blocks.Sources.Constant maxTemp(k=95) annotation (Placement(transformation(extent={{30,20},{50,40}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor kettleWall(G=5) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={10,-18})));
      Modelica.Thermal.HeatTransfer.Celsius.FixedTemperature roomTemperature(T=21) annotation (Placement(transformation(extent={{-20,-50},{0,-30}})));
      Modelica.Blocks.Logical.OnOffController onOffController(bandwidth=3) annotation (Placement(transformation(extent={{60,-10},{80,10}})));
      Modelica.Blocks.Logical.Not invert annotation (Placement(transformation(extent={{-44,44},{-56,56}})));
    equation
      connect(sineVoltage.n, resistor.n) annotation (Line(points={{-90,-10},{-90,-20},{-20,-20},{-20,-10}}, color={0,0,255}));
      connect(ground.p, resistor.n) annotation (Line(points={{-80,-40},{-80,-20},{-20,-20},{-20,-10}}, color={0,0,255}));
      connect(powerSensor.nc, resistor.p) annotation (Line(points={{-30,20},{-20,20},{-20,10}}, color={0,0,255}));
      connect(powerSensor.pv, resistor.p) annotation (Line(points={{-40,30},{-20,30},{-20,10}}, color={0,0,255}));
      connect(powerSensor.nv, resistor.n) annotation (Line(points={{-40,10},{-40,-10},{-20,-10}}, color={0,0,255}));
      connect(powerSensor.power, mean.u) annotation (Line(points={{-50,9},{-50,1.2}}, color={0,0,127}));
      connect(water.port, resistor.heatPort) annotation (Line(points={{10,12},{10,0},{-10,0}}, color={191,0,0}));
      connect(temperatureSensor.port, water.port) annotation (Line(points={{22,0},{10,0},{10,12}}, color={191,0,0}));
      connect(sineVoltage.p, switch.p) annotation (Line(points={{-90,10},{-90,20},{-80,20}}, color={0,0,255}));
      connect(switch.n, powerSensor.pc) annotation (Line(points={{-60,20},{-50,20}}, color={0,0,255}));
      connect(kettleWall.port_a, water.port) annotation (Line(points={{10,-8},{10,12}}, color={191,0,0}));
      connect(roomTemperature.port, kettleWall.port_b) annotation (Line(points={{0,-40},{10,-40},{10,-28}}, color={191,0,0}));
      connect(temperatureSensor.T, onOffController.u) annotation (Line(points={{42,0},{50,0},{50,-6},{58,-6}}, color={0,0,127}));
      connect(maxTemp.y, onOffController.reference) annotation (Line(points={{51,30},{54,30},{54,6},{58,6}}, color={0,0,127}));
      connect(onOffController.y, invert.u) annotation (Line(points={{81,0},{90,0},{90,50},{-42.8,50}}, color={255,0,255}));
      connect(invert.y, switch.control) annotation (Line(points={{-56.6,50},{-70,50},{-70,32}}, color={255,0,255}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false), graphics={Bitmap(extent={{-100,-100},{100,100}}, fileName="modelica://FM3217_2021/Resources/Images/ElectricKettleImage.jpg")}),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=500, __Dymola_NumberOfIntervals=5000));
    end ElectricKettle;

    package CoffeeAndMilk
      extends Modelica.Icons.Package;

      model Coffee "Lumped thermal element storing heat"
        Modelica.SIunits.Temperature T(start = 353.15, fixed = true, displayUnit = "degC") "Temperature of element";
        Modelica.SIunits.TemperatureSlope der_T(start = 0) "Time derivative of temperature (= der(T))";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90), visible = true, iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        parameter Modelica.SIunits.SpecificHeatCapacity Ccoffee = 4180 "Specific heat capacity of the coffee.";
        parameter Modelica.SIunits.SpecificHeatCapacity Cmilk = 3770 "Specific heat capacity of the milk.";
        Modelica.SIunits.HeatCapacity C "Heat capacity of the coffee and milk mixture";
        parameter Modelica.SIunits.Volume Vcoffee(displayUnit = "ml") = 0.0002 "Volume of the coffee.";
        parameter Modelica.SIunits.Volume Vmilk(displayUnit = "ml") = 1e-05 "Volume of the added milk.";
        parameter Modelica.SIunits.Time AddTime = 300 "The time at which milk is added.";
        parameter Modelica.SIunits.Temperature MilkTemperature = 278.15 "Temperature of the added milk.";
      equation
        C = if time > AddTime then Ccoffee * Vcoffee * 1000 + Cmilk * Vmilk * 1000 else Ccoffee * Vcoffee * 1000;
        when time > AddTime then
          reinit(T, (Ccoffee * Vcoffee * 1000 * T + Cmilk * Vmilk * 1000 * MilkTemperature) / C);
        end when;
        T = port.T;
        der_T = der(T);
        C*der(T) = port.Q_flow;
        annotation (
          Icon( graphics={Text(
                origin={70,90},
                lineColor={64,64,64},
                extent={{-70,-30},{70,10}},
                textString="%C"), Bitmap(
                origin={0,0},
                extent={{-100,-100},{100,100}},
                fileName="modelica://FM3217_2020/Resources/Images/coffee.png")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},{0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,35}},
                lineColor={0,0,0},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-6,-1},{6,-12}},
                lineColor={255,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{11,13},{50,-25}},
                lineColor={0,0,0},
                textString="T"),
              Line(points={{0,-12},{0,-96}}, color={255,0,0})}),
          Documentation(info="<html>
<p>
This is a generic model for the heat capacity of a material.
No specific geometry is assumed beyond a total volume with
uniform temperature for the entire volume.
Furthermore, it is assumed that the heat capacity
is constant (independent of temperature).
</p>
<p>
The temperature T [Kelvin] of this component is a <b>state</b>.
A default of T = 25 degree Celsius (= SIunits.Conversions.from_degC(25))
is used as start value for initialization.
This usually means that at start of integration the temperature of this
component is 25 degrees Celsius. You may, of course, define a different
temperature as start value for initialization. Alternatively, it is possible
to set parameter <b>steadyStateStart</b> to <b>true</b>. In this case
the additional equation '<b>der</b>(T) = 0' is used during
initialization, i.e., the temperature T is computed in such a way that
the component starts in <b>steady state</b>. This is useful in cases,
where one would like to start simulation in a suitable operating
point without being forced to integrate for a long time to arrive
at this point.
</p>
<p>
Note, that parameter <b>steadyStateStart</b> is not available in
the parameter menu of the simulation window, because its value
is utilized during translation to generate quite different
equations depending on its setting. Therefore, the value of this
parameter can only be changed before translating the model.
</p>
<p>
This component may be used for complicated geometries where
the heat capacity C is determined my measurements. If the component
consists mainly of one type of material, the <b>mass m</b> of the
component may be measured or calculated and multiplied with the
<b>specific heat capacity cp</b> of the component material to
compute C:
</p>
<pre>
   C = cp*m.
   Typical values for cp at 20 degC in J/(kg.K):
      aluminium   896
      concrete    840
      copper      383
      iron        452
      silver      235
      steel       420 ... 500 (V2A)
      wood       2500
</pre>
</html>"));
      end Coffee;

      model CoffeeWithMilkPort "Lumped thermal element storing heat"
        Modelica.SIunits.Temperature T(start = 353.15, fixed = true, displayUnit = "degC") "Temperature of element";
        Modelica.SIunits.Enthalpy H annotation(Dialog(group = "Variables"));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation(Placement(transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 90), visible = true, iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        parameter Modelica.SIunits.SpecificHeatCapacity Ccoffee = 4180 "Specific heat capacity of the coffee.";
        parameter Modelica.SIunits.SpecificHeatCapacity Cmilk = 3770 "Specific heat capacity of the milk.";
        Modelica.SIunits.HeatCapacity C "Heat capacity of the coffee and milk mixture";
        parameter Modelica.SIunits.Volume Vcoffee(displayUnit = "ml") = 0.0002 "Volume of the coffee.";
        parameter Modelica.SIunits.Temperature MilkTemperature = 278.15 "Temperature of the added milk.";
        Modelica.SIunits.Volume Vmilk(displayUnit = "ml", start = 0, fixed = true) "Volume of the added milk.";
        Modelica.Blocks.Interfaces.RealInput u annotation(Placement(visible = true, transformation(origin={0,80},    extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin={0,80},    extent = {{-20, -20}, {20, 20}}, rotation = -90)));
      equation
        C = Ccoffee * Vcoffee * 1000 + Cmilk * Vmilk * 1000;
        der(Vmilk) = u;
        T = port.T;
        T = H / C;
        der(H) = port.Q_flow + Cmilk*1000*u*MilkTemperature;
        annotation (
          Icon( graphics={Text(
                origin={70,90},
                lineColor={64,64,64},
                extent={{-70,-30},{70,10}},
                textString="%C"), Bitmap(
                origin={0,0},
                extent={{-100,-100},{100,100}},
                fileName="modelica://FM3217_2020/Resources/Images/coffee.png")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},{0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,35}},
                lineColor={0,0,0},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-6,-1},{6,-12}},
                lineColor={255,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{11,13},{50,-25}},
                lineColor={0,0,0},
                textString="T"),
              Line(points={{0,-12},{0,-96}}, color={255,0,0})}),
          Documentation(info="<html>
<p>
This is a generic model for the heat capacity of a material.
No specific geometry is assumed beyond a total volume with
uniform temperature for the entire volume.
Furthermore, it is assumed that the heat capacity
is constant (independent of temperature).
</p>
<p>
The temperature T [Kelvin] of this component is a <b>state</b>.
A default of T = 25 degree Celsius (= SIunits.Conversions.from_degC(25))
is used as start value for initialization.
This usually means that at start of integration the temperature of this
component is 25 degrees Celsius. You may, of course, define a different
temperature as start value for initialization. Alternatively, it is possible
to set parameter <b>steadyStateStart</b> to <b>true</b>. In this case
the additional equation '<b>der</b>(T) = 0' is used during
initialization, i.e., the temperature T is computed in such a way that
the component starts in <b>steady state</b>. This is useful in cases,
where one would like to start simulation in a suitable operating
point without being forced to integrate for a long time to arrive
at this point.
</p>
<p>
Note, that parameter <b>steadyStateStart</b> is not available in
the parameter menu of the simulation window, because its value
is utilized during translation to generate quite different
equations depending on its setting. Therefore, the value of this
parameter can only be changed before translating the model.
</p>
<p>
This component may be used for complicated geometries where
the heat capacity C is determined my measurements. If the component
consists mainly of one type of material, the <b>mass m</b> of the
component may be measured or calculated and multiplied with the
<b>specific heat capacity cp</b> of the component material to
compute C:
</p>
<pre>
   C = cp*m.
   Typical values for cp at 20 degC in J/(kg.K):
      aluminium   896
      concrete    840
      copper      383
      iron        452
      silver      235
      steel       420 ... 500 (V2A)
      wood       2500
</pre>
</html>"));
      end CoffeeWithMilkPort;

      model Milk
        Modelica.Blocks.Sources.Pulse pulse(
          nperiod=1,
          width=100,
          period=AddDuration,
          amplitude=AddedVolume/AddDuration,
          startTime=AddTime) annotation (Placement(visible=true, transformation(
              origin={0,0},
              extent={{-20,-20},{20,20}},
              rotation=0)));
        Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(
            visible=true,
            transformation(
              origin={60,0},
              extent={{-10,-10},{10,10}},
              rotation=0),
            iconTransformation(
              origin={0,-110},
              extent={{-10,-10},{10,10}},
              rotation=-90)));
        parameter Modelica.SIunits.Time AddTime=300 "The time at which milk is added.";
        parameter Modelica.SIunits.Time AddDuration=1 "Duration of the adding of milk.";
        parameter Modelica.SIunits.Volume AddedVolume(displayUnit="ml") = 1e-05 "Amount of milk added.";
      equation
        connect(pulse.y, y) annotation (Line(
            visible=true,
            points={{22,0},{60,0}},
            color={1,37,163}));
        annotation (Icon(graphics={Bitmap(extent={{-80,-100},{100,100}}, fileName="modelica://FM3217_2020/Resources/Images/milk.png")}));
      end Milk;

      package Scenarios
        extends Modelica.Icons.Package;

        model Approach1
          extends Modelica.Icons.Example;
          Coffee coffee(AddTime = AddTime) annotation(Placement(visible = true, transformation(origin = {-40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -360)));
          Modelica.Thermal.HeatTransfer.Components.ThermalConductor mug(G = 0.54) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambient(T = 293.15) annotation(Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          parameter Modelica.SIunits.Time AddTime=595   "The time at which milk is added.";
        equation
          connect(ambient.port, mug.port_b) annotation(Line(visible = true, points={{30,0},{10,0}},       color = {191, 0, 0}));
          connect(coffee.port, mug.port_a) annotation(Line(visible = true, points={{-30,0},{-10,0}},     color = {191, 0, 0}));
          annotation(experiment(StopTime=600));
        end Approach1;

        model Approach2
          extends Modelica.Icons.Example;
          CoffeeWithMilkPort coffee annotation(Placement(visible = true, transformation(origin = {-40, -0}, extent = {{-10, -10}, {10, 10}}, rotation = -360)));
          Modelica.Thermal.HeatTransfer.Components.ThermalConductor mug(G = 0.54) annotation(Placement(visible = true, transformation(origin = {0, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ambient(T = 293.15) annotation(Placement(visible = true, transformation(origin = {40, 0}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          Milk milk(AddDuration = 1, AddTime=10)   annotation(Placement(visible = true, transformation(origin = {-40, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
          connect(mug.port_b, ambient.port) annotation (
            Line(points = {{10, 0}, {30, 0}, {30, 0}, {30, 0}}, color = {191, 0, 0}));
          connect(milk.y, coffee.u) annotation (
            Line(points={{-40,19},{-40,19},{-40,8},{-40,8}},          color = {0, 0, 127}));
          connect(coffee.port, mug.port_a) annotation (
            Line(points = {{-30, 0}, {-10, 0}, {-10, 0}, {-10, 0}}, color = {191, 0, 0}));
          connect(ambient.port, mug.port_b) annotation(Line(visible = true, points={{30,0},{10,0}},       color = {191, 0, 0}));
          connect(coffee.port, mug.port_a) annotation(Line(visible = true, points={{-30,0},{-10,0}},      color = {191, 0, 0}));
          connect(milk.y, coffee.u) annotation(Line(visible = true, points={{-40,19},{-40,8},{-40,8}},                          color = {1, 37, 163}));
        end Approach2;
        annotation(Diagram(coordinateSystem(extent = {{-148.5, -105}, {148.5, 105}}, preserveAspectRatio = true, initialScale = 0.1, grid = {5, 5})));
      end Scenarios;
    end CoffeeAndMilk;

    package CupOfCoffee "Animation of the refilling process of a cup of coffee."
      model CupOfCoffee_1
        Real T(start=380);
      equation
        der(T) = -0.2*(T - 300);
      end CupOfCoffee_1;

      model CupOfCoffee_2 "Add some text"
        Real T(start=380) "Coffee Temperature";
      equation
        der(T) = -0.2*(T - 300) "Newton's Law of Cooling";
      end CupOfCoffee_2;

      model CupOfCoffee_3 "Use parameters"
        parameter Real T0=380 "Initial temp.";
        parameter Real Tamb=300 "Ambient temperature";
        parameter Real C=0.2;
        Real T(start=T0) "Coffee Temperature";
      equation
        der(T) = -C*(T - Tamb) "Newton's Law of Cooling";
      end CupOfCoffee_3;

      model CupOfCoffee_4 "More physical"
        import SI=Modelica.SIunits;
        parameter SI.Temperature T0=380 "Initial temp.";
        parameter SI.Temperature Tamb=300 "Ambient temperature";
        parameter SI.Density rho=1000 "Coffee density";
        parameter SI.SpecificHeatCapacity cv=4179 "Coffee specific heat";
        parameter SI.CoefficientOfHeatTransfer h=25 "Convection coefficient";
        parameter SI.Volume V=4e-4 "Volume of coffee";
        parameter SI.Area A=4e-3 "Area of coffee";
        SI.Temperature T(start=T0) "Coffee Temperature";
      equation
        rho*V*cv*der(T) = -h*A*(T - Tamb) "First law of thermodynamics";
      end CupOfCoffee_4;

      model CupOfCoffee_5
        ThermalCapacitance coffee(
          T0=380,
          rho=1000,
          V=4e-4,
          cv=4179);
        Convection cooling(
          h=25,
          A=4e-3,
          Tamb=300);
      equation
        connect(coffee.p, cooling.p);
      end CupOfCoffee_5;

      model CupOfCoffee_6
        ThermalCapacitance coffee(
          T0=380,
          rho=1000,
          V=4.08e-4,
          cv=4179);
        ThermalCapacitance cup(
          T0=300,
          rho=3700,
          V=8.45e-5,
          cv=880);
        Boundary cup2coffee(h=100, A=2.53e-2);
        Convection coffee_cooling(
          h=25,
          A=4e-3,
          Tamb=300);
        Convection cup_cooling(
          h=25,
          A=2.79e-2,
          Tamb=300);
      equation
        connect(coffee.p, cup2coffee.p1);
        connect(coffee.p, coffee_cooling.p);
        connect(cup.p, cup2coffee.p2);
        connect(cup.p, cup_cooling.p);
      end CupOfCoffee_6;

      model CupOfCoffee_7 "Additional physics"
        import SI=Modelica.SIunits;
        import Modelica.Constants.pi;
        parameter SI.CoefficientOfHeatTransfer h_air=25
          "Convection coefficient with air";
        parameter SI.CoefficientOfHeatTransfer h_fluid=100
          "Convection coefficient with fluid";
        parameter Customer customer;
        parameter Service service;
        parameter CoffeeProperties coffee;
        parameter CupProperties cup;
        SI.Temperature T "Coffee Temperature";
        SI.Temperature Tcup "Coffee cup temperature";
        SI.Mass M "Mass of coffee in the cup";
        SI.Area A "Surface area exposed to ambient";
        SI.Length H "Height of coffee in the cup";
        SI.Volume V "Volume of coffee in cup";
        Boolean drinking "true when drinking";
        Boolean empty "true when cup is empty";
        SI.MassFlowRate mdot_drink "drinking mass flow rate";
        discrete SI.MassFlowRate mdot_refill "refilling mass flow rate";
        SI.Energy U "coffee internal energy";
        SI.Area Acup_int "internal surface area of cup";
      initial equation
        H = 0.9*cup.H;
        T = service.Tcoffee;
        Tcup = service.Tamb;
      algorithm
        when sample(service.start_refill, service.dt_refill) then
          mdot_refill := service.dm_refill;
        end when;
        when H>=0.9*cup.H then
          mdot_refill := 0;
        end when;
      equation
        U = M*coffee.cv*T;
        der(U) = -h_air*A*(T - service.Tamb) - h_fluid*Acup_int*(T - Tcup) - mdot_drink*coffee.cp*
          T + mdot_refill*coffee.cp*service.Tcoffee
          "First law of thermodynamics for coffee";
        cup.M*cup.cv*der(Tcup) = -h_air*cup.A_ext*(Tcup - service.Tamb) - h_fluid*Acup_int*(
          Tcup - T) "First law of thermodynamics for cup";
        A = pi*cup.D^2/4 "Area of coffee exposed to air";
        Acup_int = pi*cup.D*(cup.H-H) "Area on inside of mug exposed to air";
        V = A*H "Volume of coffee in the cup";
        M = coffee.rho*V "Mass of coffee in the cup";
        der(M) = -mdot_drink + mdot_refill "Conservation of mass for coffee";
        empty = M <= 1e-9;
        drinking = mod(time, customer.dt_drink) <= customer.drink_duration and not empty;
        mdot_drink = if drinking then customer.sip_rate else 0;
      end CupOfCoffee_7;

      connector ThermalPort
        Modelica.SIunits.Temperature T;
        flow Modelica.SIunits.HeatFlowRate q;
      end ThermalPort;

      model ThermalCapacitance
        import SI=Modelica.SIunits;
        parameter SI.Temperature T0;
        parameter SI.Density rho;
        parameter SI.SpecificHeatCapacity cv;
        parameter SI.Volume V;
        ThermalPort p;
      initial equation
        p.T = T0;
      equation
        rho*V*cv*der(p.T) = p.q;
      end ThermalCapacitance;

      model Convection "Convection to the ambient"
        import SI=Modelica.SIunits;
        parameter SI.CoefficientOfHeatTransfer h;
        parameter SI.Temperature Tamb;
        parameter SI.Area A;
        ThermalPort p;
      equation
        p.q = h*A*(p.T - Tamb);
      end Convection;

      model Boundary "Convective boundary"
        import SI=Modelica.SIunits;
        parameter SI.CoefficientOfHeatTransfer h;
        parameter SI.Area A;
        ThermalPort p1;
        ThermalPort p2;
      equation
        p1.q + p2.q = 0 "No storage of energy";
        p1.q = h*A*(p1.T - p2.T);
      end Boundary;

      record Customer
        import SI=Modelica.SIunits;
        SI.Time dt_drink=60 "amount of time between sips";
        SI.MassFlowRate dm_drink=1.48e-2 "amount of mass consumed during each sip";
        SI.Time drink_duration=1 "duration for each drink";
        SI.MassFlowRate sip_rate=dm_drink/drink_duration "Rate of drinking coffee";
        SI.Time start_drink=1 "time when drinking starts";
      end Customer;

      record Service
        import SI=Modelica.SIunits;
        SI.Time start_refill=1750 "time when refilling starts";
        SI.Time dt_refill=1750 "amount of time between refills";
        SI.Mass dm_refill=0.3 "amount of mass added during each refill";
        SI.Time refill_duration=20 "duration for each refill";
        SI.Temperature Tcoffee=360 "Temperature of refill coffee";
        parameter SI.Temperature Tamb=300 "Ambient temperature";
      end Service;

      record CoffeeProperties
        import SI=Modelica.SIunits;
        parameter SI.Density rho=1000 "Coffee density";
        parameter SI.SpecificHeatCapacity cv=4179
          "Coffee specific heat (constant volume)";
        parameter SI.SpecificHeatCapacity cp=cv
          "Coffee specific heat (constant pressure)";
      end CoffeeProperties;

      record CupProperties
        import SI=Modelica.SIunits;
        import Modelica.Constants.pi;
        SI.Height H=0.127 "Height of coffee cup";
        SI.Diameter D=0.0508 "Diameter of cup base";
        SI.Density rho=3700 "density of the cup";
        SI.Mass M=rho*V "mass of the cup";
        SI.SpecificHeatCapacity cv=880 "Cup specific heat (constant volume)";
        SI.Length t=0.003175 "cup wall thickness";
        SI.Volume V=pi*t*H*(D + t)+pi*(D/2)^2*t "volume of the cup";
        SI.Area A_ext=pi*H*(D + 2*t) "external surface area of cup";
      end CupProperties;

      model CupOfCoffeeAnimation "animated model of cup of coffee"
        import SI=Modelica.SIunits;
        CupOfCoffee.CupOfCoffee_7 coffee_example;
        SI.Temperature coffeeTemp;
        SI.Temperature mugTemp;
        SI.Height h;
        function ComputeCoffeeColor
          input SI.Temperature T;
          output Real color[3];
        protected
          SI.Temperature Tmin=293;
          SI.Temperature Tmax=380;
          Real per=(min(max(Tmin, T), Tmax) - Tmin)/(Tmax - Tmin);
        algorithm
          color := {100 + 155*per,100 - 50*per,100 - 50*per};
        end ComputeCoffeeColor;

        function ComputeMugColor
          input SI.Temperature T;
          output Real color[3];
        protected
          SI.Temperature Tmin=293;
          SI.Temperature Tmax=380;
          Real per=(min(max(Tmin, T), Tmax) - Tmin)/(Tmax - Tmin);
        algorithm
          color := {100 + 155*per,100 - 50*per,100 - 50*per};
        end ComputeMugColor;
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape coffee(
          shapeType="cylinder",
          length=h,
          width=1.0,
          color=ComputeCoffeeColor(coffeeTemp),
          height=1.0,
          lengthDirection={0,1,0}) annotation (Placement(transformation(extent={{
                  -20,20},{0,40}}, rotation=0)));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape cup_bottom(
          shapeType="cylinder",
          length=-0.2,
          width=1.2,
          height=1.2,
          lengthDirection={0,1,0},
          color=ComputeMugColor(mugTemp));
        Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape cup_wall(
          shapeType="pipe",
          length=1.2,
          width=1.2,
          height=1.2,
          lengthDirection={0,1,0},
          extra=0.8,
          color=ComputeMugColor(mugTemp));
      equation
        coffeeTemp = coffee_example.T;
        mugTemp = coffee_example.Tcup;
        h=coffee_example.H/coffee_example.cup.H/0.9;
        annotation (uses, experiment(StopTime=2000, __Dymola_NumberOfIntervals=5000));
      end CupOfCoffeeAnimation;

      model TestAll
        CupOfCoffee_1 cup1;
        CupOfCoffee_2 cup2;
        CupOfCoffee_3 cup3;
        CupOfCoffee_4 cup4;
        CupOfCoffee_5 cup5;
        CupOfCoffee_6 cup6;
        CupOfCoffee_7 cup7;
      end TestAll;
      annotation (
        conversion(noneFromVersion="", noneFromVersion="1"));
    end CupOfCoffee;
  end Tutorial4;

  package Tutorial5 "Introdcution to HPL"
    model ConnectingPipes
      inner HydroPower.System_HPL system_HPL(
        steadyState=true,
        Q_start(displayUnit="m3/s"),
        constantTemperature=true) annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.SinksAndSources.Fixed_pT source1(paraOption=false) annotation (Placement(transformation(extent={{-60,60},{-40,80}})));
      HydroPower.SinksAndSources.Fixed_pT sink1(paraOption=false) annotation (Placement(transformation(extent={{60,60},{40,80}})));
      HydroPower.HydroSystems.Pipe pipe1(
        L=100,
        ZL=90,
        enable_dataVizPort_upper=false,
        enable_dataVizPort_lower=false) annotation (Placement(transformation(extent={{-8,60},{12,80}})));
    equation
      connect(source1.b, pipe1.a) annotation (Line(points={{-39,70},{-9,70}}, color={0,0,255}));
      connect(pipe1.b, sink1.b) annotation (Line(points={{13,70},{39,70}}, color={0,0,255}));
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false)),
        Diagram(coordinateSystem(preserveAspectRatio=false)),
        experiment(
          StopTime=100,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end ConnectingPipes;

    model PipeWithValve
      extends ConnectingPipes;
      HydroPower.SinksAndSources.Fixed_pT source2(paraOption=false) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
      HydroPower.SinksAndSources.Fixed_pT sink2(paraOption=false) annotation (Placement(transformation(extent={{60,20},{40,40}})));
      HydroPower.HydroSystems.PipeValve pipeValve2(
        m_dot_nom=115*1000,
        dp_nom=900000,
        ZL=90) annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      Modelica.Blocks.Sources.Ramp ramp(duration=10, startTime=10) annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
    equation
      connect(pipeValve2.a, source2.b) annotation (Line(points={{-11,30},{-39,30}}, color={0,0,255}));
      connect(pipeValve2.b, sink2.b) annotation (Line(points={{11,30},{39,30}}, color={0,0,255}));
      connect(ramp.y, pipeValve2.ValveCtrl) annotation (Line(points={{-79,50},{0,50},{0,41}}, color={0,0,127}));
      annotation (experiment(
          StopTime=100,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"), Documentation(info="<html>
<p>P = &rho; g Q H</p>
<p>P =100 MW</p>
<p>H = 90 m </p>
<p><br>Q = P/(&rho;*g*H) = 100e6 / (1e3 * 9.81 * 90)</p>
</html>"));
    end PipeWithValve;

    model SimpleWaterWay
      extends PipeWithValve(pipeValve2(L=110, ZL=100), ramp(
          height=-1,
          duration=1,
          offset=1));
      HydroPower.SinksAndSources.Fixed_pT source3(paraOption=false) annotation (Placement(transformation(extent={{-100,-20},{-80,0}})));
      HydroPower.SinksAndSources.Fixed_pT sink3(paraOption=false) annotation (Placement(transformation(extent={{100,-20},{80,0}})));
      HydroPower.HydroSystems.Pipe pipe3(
        horizontalIcon=true,
        L=10000,
        ZL=100,
        ZR=90) annotation (Placement(transformation(extent={{-50,-20},{-30,0}})));
      HydroPower.HydroSystems.HydroComponents.Containers.ClosedVolume closedVolume3(D=10) annotation (Placement(transformation(extent={{-10,-20},{10,0}})));
      HydroPower.HydroSystems.PipeValve pipeValve3(
        m_dot_nom=115*1000,
        dp_nom=900000,
        ZL=90) annotation (Placement(transformation(extent={{30,-20},{50,0}})));
    equation
      connect(sink3.b, pipeValve3.b) annotation (Line(points={{79,-10},{51,-10}}, color={0,0,255}));
      connect(pipeValve3.a, closedVolume3.b) annotation (Line(points={{29,-10},{10,-10}}, color={0,0,255}));
      connect(pipe3.b, closedVolume3.a) annotation (Line(points={{-29,-10},{-10,-10}}, color={0,0,255}));
      connect(pipe3.a, source3.b) annotation (Line(points={{-51,-10},{-79,-10}}, color={0,0,255}));
      connect(pipeValve3.ValveCtrl, pipeValve2.ValveCtrl) annotation (Line(points={{40,1},{40,10},{-70,10},{-70,50},{0,50},{0,41}}, color={0,0,127}));
      annotation (experiment(
          StopTime=100,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end SimpleWaterWay;

    model SimpleWaterWayWithSurgeTank
      extends SimpleWaterWay(ramp(
          height=-0.1,
          duration=1,
          offset=1,
          startTime=10), system_HPL(Q_start=119.4));
      HydroPower.SinksAndSources.Fixed_pT source4(paraOption=false) annotation (Placement(transformation(extent={{-100,-70},{-80,-50}})));
      HydroPower.HydroSystems.Pipe pipe4(
        horizontalIcon=true,
        L=10000,
        ZL=100,
        ZR=90) annotation (Placement(transformation(extent={{-50,-70},{-30,-50}})));
      HydroPower.SinksAndSources.Fixed_pT sink4(paraOption=false) annotation (Placement(transformation(extent={{100,-70},{80,-50}})));
      HydroPower.HydroSystems.PipeValve pipeValve4(
        m_dot_nom=115*1000,
        dp_nom=900000,
        ZL=90) annotation (Placement(transformation(extent={{30,-70},{50,-50}})));
      HydroPower.HydroSystems.SurgeTank surgeTank4(
        D=10,
        deltZ=100,
        Vol=100) annotation (Placement(transformation(extent={{-10,-70},{10,-50}})));
    equation
      connect(pipeValve4.b, sink4.b) annotation (Line(points={{51,-60},{79,-60}}, color={0,0,255}));
      connect(pipe4.a, source4.b) annotation (Line(points={{-51,-60},{-79,-60}}, color={0,0,255}));
      connect(pipe4.b, surgeTank4.a) annotation (Line(points={{-29,-60},{-11,-60}}, color={0,0,255}));
      connect(surgeTank4.b, pipeValve4.a) annotation (Line(points={{11,-60},{29,-60}}, color={0,0,255}));
      connect(pipeValve4.ValveCtrl, pipeValve2.ValveCtrl) annotation (Line(points={{40,-49},{40,-40},{-70,-40},{-70,50},{0,50},{0,41}}, color={0,0,127}));
      annotation (experiment(
          StopTime=1000,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end SimpleWaterWayWithSurgeTank;
  end Tutorial5;

  package Tutorial6 "Waterway with reservoir"
    model ReservoirBase
      inner HydroPower.System_HPL system_HPL(steadyState=true, constantTemperature=true) annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.HydroSystems.Reservoir headwater annotation (Placement(transformation(extent={{-70,40},{-50,60}})));
      HydroPower.HydroSystems.Reservoir tailwater annotation (Placement(transformation(extent={{50,20},{70,40}})));
      HydroPower.HydroSystems.Pipe conduit(horizontalIcon=true) annotation (Placement(transformation(extent={{-40,34},{-20,54}})));
    equation
      connect(headwater.a2_pipe, conduit.a) annotation (Line(points={{-49,44},{-41,44}}, color={0,0,255}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)));
    end ReservoirBase;

    model TwoReservoirs
      extends Modelica.Icons.Example;
      extends ReservoirBase(conduit(L=1000, ZL=100));
    equation
      connect(conduit.b, tailwater.a1_pipe) annotation (Line(points={{-19,44},{20,44},{20,24},{49,24}}, color={0,0,255}));
      annotation (experiment(
          StopTime=600,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end TwoReservoirs;

    model TwoReservoirsWithSource
      extends TwoReservoirs;
      HydroPower.SinksAndSources.Fixed_HT constantWaterHead(
        paraOption=false,
        H_const=headwater.H_start[1],
        Hmax=headwater.Hmax[1],
        depth=headwater.depth[1]) annotation (Placement(transformation(extent={{-100,46},{-80,66}})));
      HydroPower.SinksAndSources.Fixed_HT constantWaterTail(
        paraOption=false,
        H_const=tailwater.H_start[tailwater.n],
        Hmax=tailwater.Hmax[tailwater.n],
        depth=tailwater.depth[tailwater.n]) annotation (Placement(transformation(extent={{100,26},{80,46}})));
    equation
      connect(constantWaterHead.b, headwater.a1_open) annotation (Line(points={{-79,56},{-71,56}}, color={0,0,255}));
      connect(constantWaterTail.b, tailwater.a2_open) annotation (Line(points={{79,36},{71,36}}, color={0,0,255}));
    end TwoReservoirsWithSource;

    model WaterWayRes
      extends Modelica.Icons.Example;
      extends FM3217_2021.Tutorial6.ReservoirBase(conduit(
          L=1000,
          ZL=100,
          ZR=90));
      HydroPower.SinksAndSources.Fixed_HT constantWaterHead(
        paraOption=false,
        H_const=75,
        Hmax=headwater.Hmax[1],
        depth=headwater.depth[1]) annotation (Placement(transformation(extent={{-100,46},{-80,66}})));
      HydroPower.SinksAndSources.Fixed_HT constantWaterTail(
        paraOption=false,
        H_const=75,
        Hmax=tailwater.Hmax[tailwater.n],
        depth=tailwater.depth[tailwater.n]) annotation (Placement(transformation(extent={{100,26},{80,46}})));
      HydroPower.HydroSystems.PipeValve pipeValve(
        m_dot_nom=113.26*1000,
        dp_nom=900000,
        d_nom(displayUnit="kg/m3"),
        ZL=90) annotation (Placement(transformation(extent={{20,30},{40,50}})));
      HydroPower.HydroSystems.SurgeTank surgeTank(
        D=10,
        deltZ=50,
        Vol=100)                                                      annotation (Placement(transformation(extent={{-10,34},{10,54}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=0.9,                     duration=30,
        offset=0.1,
        startTime=100)                                           annotation (Placement(transformation(extent={{-10,70},{10,90}})));
    equation
      connect(constantWaterHead.b, headwater.a1_open) annotation (Line(points={{-79,56},{-71,56}}, color={0,0,255}));
      connect(tailwater.a2_open, constantWaterTail.b) annotation (Line(points={{71,36},{79,36}}, color={0,0,255}));
      connect(surgeTank.b, pipeValve.a) annotation (Line(points={{11,44},{16,44},{16,40},{19,40}}, color={0,0,255}));
      connect(pipeValve.b, tailwater.a1_pipe) annotation (Line(points={{41,40},{44,40},{44,24},{49,24}}, color={0,0,255}));
      connect(conduit.b, surgeTank.a) annotation (Line(points={{-19,44},{-11,44}}, color={0,0,255}));
      connect(ramp.y, pipeValve.ValveCtrl) annotation (Line(points={{11,80},{30,80},{30,51}}, color={0,0,127}));
    end WaterWayRes;

    model WaterWayResClosingValve
      extends FM3217_2021.Tutorial6.WaterWayRes(ramp(height=-0.9,
                                                               offset=1));
    end WaterWayResClosingValve;

    model SundsbarmWaterway
      extends Modelica.Icons.Example;
      inner HydroPower.System_HPL system_HPL(steadyState=true,
        Q_start=24,                                            constantTemperature=true) annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.HydroSystems.Reservoir headwater(
        Hmax=ones(headwater.n)*(564 + 48 + 5),
        depth=ones(headwater.n)*(48 + 5),
        H_start=ones(headwater.n)*(564 + 48))     annotation (Placement(transformation(extent={{-90,20},{-70,40}})));
      HydroPower.HydroSystems.Reservoir tailwater(
        Hmax=ones(tailwater.n)*(110 + 5 + 3),
        depth=ones(tailwater.n)*(5 + 3),
        H_start=ones(tailwater.n)*(110 + 5))      annotation (Placement(transformation(extent={{70,0},{90,20}})));
      HydroPower.HydroSystems.Pipe conduit(
        endD={5.8,5.8},
        ZL=564,
        ZR=541.5,
        horizontalIcon=true,
        L=6600)                                                         annotation (Placement(transformation(extent={{-60,14},{-40,34}})));
      HydroPower.SinksAndSources.Fixed_HT constantWaterHead(
        paraOption=false,
        H_const=564 + 48,
        Hmax=564 + 48 + 5,
        depth=48 + 5)
                  annotation (Placement(transformation(extent={{-70,50},{-90,70}})));
      HydroPower.SinksAndSources.Fixed_HT constantTailWater(
        paraOption=false,
        H_const=110 + 5,
        Hmax=110 + 5 + 3,
        depth=5 + 3)
                  annotation (Placement(transformation(extent={{70,28},{90,48}})));
      HydroPower.HydroSystems.SurgeTank surgeTank(
        D=3.6,
        deltZ=150,
        H2L=0.87,
        Vol=100) annotation (Placement(transformation(extent={{-30,14},{-10,34}})));
      HydroPower.HydroSystems.PipeValve pressureShaft(
        endD={3,3},
        m_dot_nom=24e3,
        dp_nom=4890000,
        L=724,
        ZL=541.5,
        ZR=112.5) annotation (Placement(transformation(extent={{0,10},{20,30}})));
      Modelica.Blocks.Sources.Ramp ramp(
        height=-0.9,
        duration=10,
        offset=1,
        startTime=100) annotation (Placement(transformation(extent={{-20,50},{0,70}})));
      HydroPower.HydroSystems.Pipe tailRace(
        endD={5.8,5.8},
        ZL=110.5,
        ZR=110,
        horizontalIcon=true,
        L=600)                                                          annotation (Placement(transformation(extent={{40,-6},{60,14}})));
      HydroPower.HydroSystems.HydroComponents.Containers.ClosedVolume turbineHouse(D=5.8, L=2) annotation (Placement(transformation(extent={{24,6},{36,18}})));
    equation
      connect(headwater.a2_pipe,conduit. a) annotation (Line(points={{-69,24},{-61,24}}, color={0,0,255}));
      connect(headwater.a1_open, constantWaterHead.b) annotation (Line(points={{-91,36},{-96,36},{-96,60},{-91,60}}, color={0,0,255}));
      connect(tailwater.a2_open, constantTailWater.b) annotation (Line(points={{91,16},{96,16},{96,38},{91,38}}, color={0,0,255}));
      connect(conduit.b, surgeTank.a) annotation (Line(points={{-39,24},{-31,24}}, color={0,0,255}));
      connect(surgeTank.b, pressureShaft.a) annotation (Line(points={{-9,24},{-4,24},{-4,20},{-1,20}}, color={0,0,255}));
      connect(pressureShaft.ValveCtrl, ramp.y) annotation (Line(points={{10,31},{10,60},{1,60}}, color={0,0,127}));
      connect(tailwater.a1_pipe, tailRace.b) annotation (Line(points={{69,4},{61,4}}, color={0,0,255}));
      connect(tailRace.a, turbineHouse.b) annotation (Line(points={{39,4},{38,4},{38,12},{36,12}}, color={0,0,255}));
      connect(pressureShaft.b, turbineHouse.a) annotation (Line(points={{21,20},{22,20},{22,12},{24,12}}, color={0,0,255}));
      annotation (
        experiment(
          StopTime=600,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end SundsbarmWaterway;
  end Tutorial6;

  package Tutorial7 "First model of Sundsbarm"
    model PlantConnectAndDisconnectToGrid
      extends HydroPower.Examples.PlantConnectAndDisconnectToGrid;
      annotation (experiment(
          StopTime=600,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end PlantConnectAndDisconnectToGrid;

    model Sundsbarm "Model of the Sundsbarm Power Station"
      extends Modelon.Icons.Experiment;

      HydroPower.MechanicalSystems.BasicTurbine turbine(
        np=12,
        H_nom=480,
        tableOnFile=true,
        LdraftTube=10,
        DavDraftTube=2,
        LscrollCase=5,
        DavScrollCasing=2,
        PUInFlowTables=true,
        QTableName="Qtab",
        Q_nom=24,
        H_start=564 + 48,
        H_start_draftTube=115,
        Ty=0.4,
        yvLim1=[-0.1,0.1],
        yvLim2=[-0.2,0.2],
        TurbineDataFile=Modelica.Utilities.Files.loadResource(HydroPower.TABLE_DIR
             + "TurbineDataFile.mat"),
        P_nom=103000000)
                        annotation (Placement(transformation(extent={{20,-50},{40,-30}},
                     rotation=0)));

      HydroPower.ElectricalSystems.PowerGrid powerGrid(
        startTime=1e6,
        distNoGen={-2,0,0},
        distTgen={150,1e6,1e6}) annotation (Placement(transformation(extent={{-90,40},{-70,60}}, rotation=0)));

      Modelica.Blocks.Sources.Ramp pwr_ref(
        duration=10,
        height=0,
        offset=45e6,
        startTime=1e6) annotation (Placement(transformation(extent={{-37,64},{-25,76}},
                      rotation=0)));
      HydroPower.ElectricalSystems.GeneratorAndMCB generator(
        np={12},
        f_start=0,
        J={83e3},
        timeMCB_close={150},
        timeMCB_open={200},
        P_nom={103000000})
                          annotation (Placement(transformation(extent={{-60,40},{-40,60}},
                       rotation=0)));
      HydroPower.ControllersAndSensors.TurbineGovernorAnalog turbineGovernor(
        K_noLoad=0.8,
        Ki_noLoad=0.2,
        Kd_noLoad=0.1,
        ep=1,
        tRamp=40,
        P_generator_nom=generator.P_nom[1],
        enableRamp=false,
        K_load=0.4,
        Kd_load=0.1,
        Ki_load=0.2)      annotation (Placement(transformation(extent={{-24,40},{-4,60}}, rotation=0)));
      HydroPower.Visualizers.RealValue turbinePower(precision=2, input_Value=turbine.summary.P_turbine*1e-6) annotation (Placement(transformation(extent={{64,10},{78,24}})));
      HydroPower.Visualizers.BooleanIndicator MCB(input_Value=turbineGovernor.summary.isMCB) annotation (Placement(transformation(extent={{64,73},{77,87}})));
      HydroPower.Visualizers.RealValue gridbalanceNum(precision=2, input_Value=generator.summary.P_grid_tot*1e-6) annotation (Placement(transformation(extent={{64,38},{78,52}})));
      inner HydroPower.System_HPL system_HPL(
        steadyState=true, constantTemperature=true)
                     annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      HydroPower.Visualizers.RealValue gridbalanceNum1(precision=2, input_Value=generator.summary.f[1]) annotation (Placement(transformation(extent={{64,53},{78,67}})));
      HydroPower.Visualizers.RealValue gridbalanceNum2(precision=2, input_Value=generator.summary.P_generator[1]*1e-6) annotation (Placement(transformation(extent={{64,25},{78,39}})));
      HydroPower.HydroSystems.Pipe      pressureShaft(
        L=724,
        ZR=112.5,
        endD={3,3},
        ZL=541.5,
        enable_dataVizPort_lower=false)
        annotation (Placement(transformation(extent={{-10,-50},{10,-30}})));
      HydroPower.HydroSystems.Pipe conduit(
        horizontalIcon=true,
        ZL=564,
        endD={5.8,5.8},
        L=6600,
        ZR=541.5)
        annotation (Placement(transformation(extent={{-70,-50},{-50,-30}})));
      HydroPower.HydroSystems.SurgeTank surgeTank(
        deltZ=150,
        H2L=0.87,
        Vol=100,
        D=3.6)
        annotation (Placement(transformation(extent={{-40,-50},{-20,-30}})));
      HydroPower.HydroSystems.Reservoir reservoir(
        depth=ones(reservoir.n)*(48 + 5),
        H_start=fill((564 + 48), reservoir.n),
        Hmax=fill((564 + 48 + 5), reservoir.n)) annotation (Placement(transformation(extent={{-95,-44},{-75,-24}})));
      HydroPower.HydroSystems.Reservoir river(
        depth=ones(river.n)*(5 + 3),
        H_start=fill((115), river.n),
        Hmax=fill((110 + 5 + 3), river.n)) annotation (Placement(transformation(extent={{75,-44},{95,-24}})));
      HydroPower.HydroSystems.Pipe downstream(
        horizontalIcon=true,
        L=600,
        ZL=110.5,
        endD={5.8,5.8},
        ZR=110) annotation (Placement(transformation(extent={{50,-50},{70,-30}})));
      HydroPower.SinksAndSources.Fixed_HT waterSource(
        paraOption=false,
        H_const=reservoir.H_start[reservoir.n],
        Hmax=reservoir.Hmax[reservoir.n],
        depth=reservoir.depth[1])
                  annotation (Placement(transformation(extent={{-70,-18},{-90,2}})));
      HydroPower.SinksAndSources.Fixed_HT waterSink(
        paraOption=false,
        H_const=river.H_start[river.n],
        Hmax=river.Hmax[river.n],
        depth=river.depth[river.n])
        annotation (Placement(transformation(extent={{70,-18},{90,2}})));
    equation

      connect(powerGrid.f_grid, generator.f_grid) annotation (Line(
          points={{-69,43},{-61,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(powerGrid.P_grid_balance, generator.P_grid_balance) annotation (Line(
          points={{-69,57},{-61,57}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pwr_ref.y, turbineGovernor.P_reference) annotation (Line(
          points={{-24.4,70},{-20,70},{-20,61}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.f_out[1], turbineGovernor.f) annotation (Line(
          points={{-39,43},{-25,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.onMCB, powerGrid.MCB) annotation (Line(
          points={{-50,61},{-50,70},{-80,70},{-80,61}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.onMCB[1], turbineGovernor.isMCB) annotation (Line(
          points={{-50,61},{-50,80},{-14,80},{-14,61}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.P_out[1], turbineGovernor.P_generator) annotation (Line(
          points={{-39,57},{-25,57}},
          color={0,0,127},
          smooth=Smooth.None));

      connect(generator.f_out, powerGrid.f) annotation (Line(
          points={{-39,43},{-33,43},{-33,6},{-94,6},{-94,43},{-91,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(powerGrid.J_grid, generator.J_grid)
        annotation (Line(points={{-69,50},{-61,50}},          color={0,0,127}));
      connect(turbineGovernor.y, turbine.yGV)
        annotation (Line(points={{-3,50},{36,50},{36,-29}}, color={0,0,127}));
      connect(generator.f_out[1], turbine.f_generator) annotation (Line(points={{-39,43},{-33,43},{-33,6},{24,6},{24,-29}},
                                                         color={0,0,127}));
      connect(turbine.TurbineData[1], generator.P_turbine[1]) annotation (Line(
            points={{30,-29.6667},{30,34},{-56,34},{-56,39}},            color={0,0,
              127}));
      connect(surgeTank.a,conduit. b) annotation (Line(
          points={{-41,-40},{-49,-40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(surgeTank.b,pressureShaft. a) annotation (Line(
          points={{-19,-40},{-11,-40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(reservoir.a2_pipe,conduit. a) annotation (Line(
          points={{-74,-40},{-71,-40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(downstream.b,river. a1_pipe) annotation (Line(points={{71,-40},{74,-40}},
                                  color={0,0,255}));
      connect(waterSource.b,reservoir. a1_open) annotation (Line(points={{-91,-8},{-100,-8},{-100,-28},{-96,-28}},
                                      color={0,0,255}));
      connect(waterSink.b,river. a2_open) annotation (Line(points={{91,-8},{100,-8},{100,-28},{96,-28}},
                                       color={0,0,255}));
      connect(turbine.a, pressureShaft.b) annotation (Line(points={{19,-40},{11,-40}}, color={0,0,255}));
      connect(turbine.b, downstream.a) annotation (Line(points={{41,-40},{49,-40}}, color={0,0,255}));
      annotation (
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={1,1}), graphics={
            Rectangle(
              extent={{50,91},{90,7}},
              lineColor={215,215,215},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid,
              radius=2),           Text(
              extent={{54,13},{90,9}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="turbine power [MW]"),Text(
              extent={{57,73},{85,65}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="Main Circuit Breaker"),
                                Text(
              extent={{39,41},{103,37}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="grid power balance [MW]"),
                                Text(
              extent={{39,56},{103,52}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="generator frequency [Hz]"),
                                Text(
              extent={{38,27},{102,23}},
              lineColor={0,0,0},
              lineThickness=0.5,
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              textString="generator power [MW]"),
            Bitmap(extent={{-11,-99},{77,-56}}, fileName="modelica://FM3217_2021/Resources/Images/SundsbarmTurbine.png"),
            Bitmap(extent={{-90,-15},{20,58}}, fileName="modelica://FM3217_2021/Resources/Images/SundsbarmGenerator.png")}),
        experiment(
          StopTime=600,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"),
        Documentation(info="<html>
<p>This example illustrates a hydro power plant acting under no load and when connected to the power grid. </p>
<h4>Model experiment description</h4>
<p>When there is no load present the governor will only have the frequency error as input signal which will have the effect that the frequency of the hydro plant generator is controlled to equal the nominal frequency. This behaviour can be seen during the first 150s of simulation. </p>
<p>When 150s has passed and the frequency of the generator is synchronized to the grid frequency, the MCB is closed, the power reference is set to 45MW and new PID parameters are applied. </p>
<p>At time=350 load rejection takes place and the MCB opens once again. </p>
<h4>Simulation setup</h4>
<p>Simulate for 600s using solver Radau with a tolerance set to 1e-6.</p>
<h4>Output</h4>
<p>The most interesting variables are:</p>
<ul>
<li>generator frequency - generator.summary.f[1]</li>
<li>generated power - generator.summary.P_generator[1]</li>
</ul>
</html>",     revisions="<html>
<hr><p><font color=\"#E72614\"><b>Copyright &copy; 2004-2017, MODELON AB</b></font> <font color=\"#AFAFAF\">The use of this software component is regulated by the licensing conditions for the HydroPower Library. This copyright notice must, unaltered, accompany all components that are derived from, copied from, or by other means have their origin from the HydroPower Library. </font>
</html>"));
    end Sundsbarm;
  end Tutorial7;

  package Tutorial8 "Power flow calculations"
    model EquivalentPowerFlow
      Modelica.Mechanics.Rotational.Components.Inertia grid(J=100, w(displayUnit="Hz", start=314.15926535898)) annotation (Placement(transformation(extent={{-40,-40},{40,40}})));
      Modelica.Mechanics.Rotational.Sources.ConstantTorque generator(tau_constant=10) annotation (Placement(transformation(extent={{-90,-10},{-70,10}})));
      Modelica.Mechanics.Rotational.Sources.ConstantTorque load(tau_constant=-15) annotation (Placement(transformation(extent={{90,-10},{70,10}})));
    equation
      connect(generator.flange, grid.flange_a) annotation (Line(points={{-70,0},{-40,0}}, color={0,0,0}));
      connect(load.flange, grid.flange_b) annotation (Line(points={{70,0},{40,0}}, color={0,0,0}));
    end EquivalentPowerFlow;

    model FrequencyCalc
      Modelica.Mechanics.Rotational.Components.Inertia inertia(J=2*200e6/(2*Modelica.Constants.pi*50)^2, w(
          displayUnit="Hz",
          fixed=true,
          start=314.15926535898)) annotation (Placement(transformation(extent={{-62,70},{-42,90}})));
      Modelica.Mechanics.Rotational.Sources.ConstantTorque load(tau_constant=-10e6/(2*Modelica.Constants.pi*50)) annotation (Placement(transformation(extent={{80,70},{60,90}})));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{10,70},{30,90}})));
      Modelica.Blocks.Continuous.Integrator energy annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={12,50})));
    equation
      connect(inertia.flange_b, powerSensor.flange_a) annotation (Line(points={{-42,80},{10,80}}, color={0,0,0}));
      connect(powerSensor.flange_b, load.flange) annotation (Line(points={{30,80},{60,80}}, color={0,0,0}));
      connect(energy.u, powerSensor.power) annotation (Line(points={{12,62},{12,69}}, color={0,0,127}));
      annotation (Diagram( graphics={Text(
              extent={{-74,66},{-22,56}},
              lineColor={238,46,47},
              textString="= 200 MJ @50Hz"), Text(
              extent={{6,74},{128,22}},
              lineColor={238,46,47},
              textString="= 10 MW 
load for 1 sec
= 10 MJ")}),
        Documentation(info="<html>
<h4>Example</h4>
<p>K<sub>1</sub> = 200 MJ</p>
<p>P<sub>L</sub> = 10 MW @ 1sec</p>
<h5>Question:</h5>
<p>If the Inertia was rotating at 50 Hz in the beginning what is the speed/frequency after 1 sec?</p>
<h5>Solution:</h5>
<p><br>K<sub>2</sub> = K<sub>1</sub> - P<sub>L</sub> &middot; 1 s = 200 MJ - 10 MW &middot; 1 s = 200 MJ - 10 MJ = 190 MJ</p>
<p><br>K = 1/2 &middot; J &middot; w<sup>2</sup> = 2 &middot; J &middot; &pi;<sup>2</sup> &middot; f<sup>2</sup></p>
<p><br>K<sub>1</sub>/K<sub>2</sub> = 200 MJ / 190 MJ = f<sub>1<sup>2</sup>/f<sub>2<sup>2</sup> </p>
<blockquote>
<p><br>&rArr; f<sub>2</sub> = &radic;(50<sup>2</sup> &middot; 190/200) = <b>48.73 Hz</b></p>
</blockquote>
<h5>Tip:</h5>
<p>J = 2 &middot; K /w<sup>2</sup></p>
<p>T = P/w</p>
</html>"));
    end FrequencyCalc;

    model FrequencyCalcCorrect
      extends FrequencyCalc;
      Modelica.Mechanics.Rotational.Components.Inertia inertia1(J=2*200e6/(2*Modelica.Constants.pi*50)^2, w(
          displayUnit="Hz",
          start=314.15926535898,
          fixed=true)) annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
      Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor1 annotation (Placement(transformation(extent={{-20,-10},{0,10}})));
      Modelica.Mechanics.Rotational.Sources.Torque torque annotation (Placement(transformation(extent={{30,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant const(k=-10e6) annotation (Placement(transformation(extent={{92,0},{72,20}})));
      Modelica.Blocks.Math.Division division annotation (Placement(transformation(extent={{56,-6},{44,6}})));
      Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-50,-20})));
      Modelica.Blocks.Continuous.Integrator energy1 annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-18,-30})));
    equation
      connect(inertia1.flange_b,powerSensor1. flange_a) annotation (Line(points={{-60,0},{-20,0}},     color={0,0,0}));
      connect(powerSensor1.flange_b,torque. flange) annotation (Line(points={{0,0},{10,0}},      color={0,0,0}));
      connect(const.y,division. u1) annotation (Line(points={{71,10},{66,10},{66,3.6},{57.2,3.6}},       color={0,0,127}));
      connect(torque.tau,division. y) annotation (Line(points={{32,0},{43.4,0}},     color={0,0,127}));
      connect(inertia1.flange_b,speedSensor. flange) annotation (Line(points={{-60,0},{-50,0},{-50,-10}},     color={0,0,0}));
      connect(speedSensor.w,division. u2) annotation (Line(points={{-50,-31},{-50,-60},{66,-60},{66,-3.6},{57.2,-3.6}},   color={0,0,127}));
      connect(powerSensor1.power, energy1.u) annotation (Line(points={{-18,-11},{-18,-18}}, color={0,0,127}));
    end FrequencyCalcCorrect;

    model ElectricalPowerFlowCalc
      Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-70,0},{-50,20}})));
      Modelica.Electrical.Analog.Basic.Inductor Xs annotation (Placement(transformation(extent={{-20,60},{0,80}})));
      Modelica.Electrical.Analog.Sources.SineVoltage Ea(
        V=sqrt(2)*230,
        phase=0,
        freqHz=50) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={-60,50})));
      Modelica.Electrical.Analog.Sources.SineVoltage Vterminal(V=sqrt(2)*230, freqHz=50) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={60,50})));
      Modelica.Electrical.Analog.Sensors.PowerSensor powerSensor annotation (Placement(transformation(extent={{20,60},{40,80}})));
      Modelica.Electrical.MultiPhase.Basic.Inductor inductor(L={10,10,10}) annotation (Placement(transformation(extent={{-20,-30},{0,-10}})));
      Modelica.Electrical.MultiPhase.Sources.SineVoltage EaM(
        V=0.9*sqrt(2)*{230,230,230},
        phase=-Modelica.Electrical.MultiPhase.Functions.symmetricOrientation(3) + 0*{30,30,30}/180*Modelica.Constants.pi,
        freqHz={50,50,50}) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-60,-40})));
      Modelica.Electrical.MultiPhase.Sources.SineVoltage VterminalM(V=sqrt(2)*{230,230,230}, freqHz={50,50,50}) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=-90,
            origin={60,-40})));
      Modelica.Electrical.MultiPhase.Basic.Star star annotation (Placement(transformation(
            extent={{-10,-11},{10,11}},
            rotation=-90,
            origin={-60,-75})));
      Modelica.Electrical.Analog.Basic.Ground ground1
                                                     annotation (Placement(transformation(extent={{-70,-106},{-50,-86}})));
      Modelica.Electrical.MultiPhase.Sensors.AronSensor P annotation (Placement(transformation(extent={{10,-30},{30,-10}})));
      Modelica.Electrical.MultiPhase.Sensors.ReactivePowerSensor Q annotation (Placement(transformation(extent={{34,-30},{54,-10}})));
    equation
      connect(Xs.p, Ea.p) annotation (Line(points={{-20,70},{-60,70},{-60,60}}, color={0,0,255}));
      connect(Ea.n, ground.p) annotation (Line(points={{-60,40},{-60,36},{-60,36},{-60,20}},
                                                                           color={0,0,255}));
      connect(Vterminal.n, ground.p) annotation (Line(points={{60,40},{60,20},{-60,20}}, color={0,0,255}));
      connect(Xs.n, powerSensor.pc) annotation (Line(points={{0,70},{20,70}}, color={0,0,255}));
      connect(powerSensor.nc, Vterminal.p) annotation (Line(points={{40,70},{60,70},{60,60}}, color={0,0,255}));
      connect(powerSensor.nv, ground.p) annotation (Line(points={{30,60},{30,20},{-60,20}}, color={0,0,255}));
      connect(powerSensor.pv, Vterminal.p) annotation (Line(points={{30,80},{30,90},{60,90},{60,60}}, color={0,0,255}));
      connect(EaM.plug_p, inductor.plug_p) annotation (Line(points={{-60,-30},{-60,-20},{-20,-20}}, color={0,0,255}));
      connect(star.plug_p, EaM.plug_n) annotation (Line(points={{-60,-65},{-60,-50}}, color={0,0,255}));
      connect(VterminalM.plug_n, EaM.plug_n) annotation (Line(points={{60,-50},{60,-60},{-60,-60},{-60,-50}}, color={0,0,255}));
      connect(star.pin_n, ground1.p) annotation (Line(points={{-60,-85},{-60,-86}}, color={0,0,255}));
      connect(inductor.plug_n, P.plug_p) annotation (Line(points={{0,-20},{10,-20}}, color={0,0,255}));
      connect(P.plug_n, Q.plug_p) annotation (Line(points={{30,-20},{34,-20}}, color={0,0,255}));
      connect(Q.plug_n, VterminalM.plug_p) annotation (Line(points={{54,-20},{60,-20},{60,-30}}, color={0,0,255}));
      annotation (experiment(StopTime=0.1));
    end ElectricalPowerFlowCalc;
  end Tutorial8;

  package Tutorial9

    model SampleCalculations
      extends Modelica.Icons.Information;
      annotation (preferredView="info",Documentation(info="<html>
<h4>Sample calculations for production balance in the Power Grid</h4>
<h5>Default values</h5>
<p>P_grid = 1000 MW</p>
<p>P_1 = 0.45 * P_grid = 450 MW</p>
<p>P_1units = 20</p>
<p>P_1perUnit = P_1 / P_1units = 450 MW / 20 = 22.5 MW</p>
<p><u><b>&Delta;P</b></u> = distNoGen * P_1perUnit = -2 * 22.5 =<b> <u>45 MW</b></u></p>
<h5>Manual Droop calculations for DroopSimulation</h5>
<p>R = - (&Delta;F/F_r) / (&Delta;P/P_r)</p>
<p>&Delta;P = - P_r * (&Delta;F/F_r) / R</p>
<p>&Delta;P = - 45 MW * (-0.4043/50) / 0.1</p>
<p><u>&Delta;P = 3.639 MW (theory)</u></p>
<p><u>&Delta;P = 3.167 MW (simulation)</u></p>
</html>"));
    end SampleCalculations;

    model ResitiveLoad
      extends HydroPower.Examples.PlantConnectAndDisconnectToGrid(
        powerGrid(
          loadDiv={1,0,0},
          enableDroop=false,
          distTgen={150,1000,1e6},
          distNoGen={-1,-50,0}),
        pwr_ref(offset=45e6),
        turbineGovernor(enableDroop=false),
        generator(timeMCB_open={1e6}));

      annotation (experiment(
          StopTime=2000,
          __Dymola_NumberOfIntervals=5000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"),                                                    preferredView="info", Documentation(info="<html>
<p>The model demonstrates the behaviour of an electrical generator in combination with an electrical grid where we only have pure resistive loads present. </p>
<h4>Task</h4>
<p>Apply the following changes:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor: </li>
<li><ul>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid: </li>
<li><ul>
<li>Only <b>resistive</b> loads: <span style=\"font-family: monospace;\">loadDiv = {1,0,0}</span></li>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li>Remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>Remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s</span> from production</li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>How does the frequency behave?</li>
<li>What is the power difference &Delta;P of the generator when comparing the power produced at <span style=\"font-family: monospace;\">t=1000s</span> and the power produced by the generator at <span style=\"font-family: monospace;\">t=2000s</span>?</li>
<li>Change the power set point of the turbine governor to 45 MW and look at the frequency behaviour.</li>
</ol>
</html>"));
    end ResitiveLoad;

    model DroopSimulations
      extends ResitiveLoad(
        pwr_ref(offset=22.5e6),
        turbineGovernor(enableDroop=true, ep=0.05),
        powerGrid(P_grid=100000000, distNoGen={-10,-50,0}));
      annotation (experiment(
          StopTime=10000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"),                                                   preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Enable the droop: <span style=\"font-family: monospace;\">enableDroop = true</span> </li>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid: </li>
<li><ul>
<li>Only resistive loads: loadDiv = {1,0,0}</li>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li>Remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>Remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s </span>from production</li>
</ul></li>
</ul>
<p><br>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine governor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>What is the power difference &Delta;P of the generator when comparing the power produced at <span style=\"font-family: monospace;\">t=1000s</span> and the power produced by the generator at <span style=\"font-family: monospace;\">t=2000s</span>?</li>
<li>What is &Delta;f in this case and does &Delta;P fit with the theoretical obtained one from the calculation with the droop factor <span style=\"font-family: monospace;\">ep</span>?</li>
<li>Now re-run the simulation but set the droop for the turbine governor to <span style=\"font-family: monospace;\">ep = 0.05</span>. What difference does this make?</li>
<li>What does the <span style=\"font-family: monospace;\">DeadBand</span> setting in the turbine governor mean and what difference does it make to simulate with <span style=\"font-family: monospace;\">DeadBand=0</span>?</li>
<li>Try the simulation for a smaller power grid, e.g., <span style=\"font-family: monospace;\">P_grid = 100 MW</span>. Pay attention that you also need to adjust the number of units that are changed in order to maintain the 22.5 MW change and make the second change only 22.5 MW. How does the result differ now?</li>
</ol>
</html>"),
        experiment(StopTime=2000, Algorithm="Radau"));
    end DroopSimulations;
  end Tutorial9;

  package Tutorial10
    model LoadChanges
      extends HydroPower.Examples.PlantConnectAndDisconnectToGrid;
      annotation (preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid:</li>
<li><ul>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li><b>resistive, f and f*f dependent loads </b>(default <span style=\"font-family: monospace;\">loadDiv={0.5,0.25,0.25}</span>)</li>
<li>remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s </span>from production</li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine govenor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>Compare the behaviour of the frequency with the one obtained from model a load disttribution of <span style=\"font-family: monospace;\">loadDiv={1,0,0}</span>?</li>
<li>In case you do not see much difference run another simulation with <span style=\"font-family: monospace;\">loadDiv={0,0,1}</span> and compare again the frequency with the previous simulations.</li>
</ol>
</html>"),
        experiment(
          StopTime=2000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end LoadChanges;

    model ProductionDroop
      extends LoadChanges;
      annotation (preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid:</li>
<li><ul>
<li>enable the droop for the production: <span style=\"font-family: monospace;\">enableDroop = true</span></li>
<li><b>set the droop of the production to <span style=\"font-family: monospace;\">ep = {0.1,0.08,0.04}</span></b></li>
<li>default resistive, f and f*f dependent loads: <span style=\"font-family: monospace;\">loadDiv={0.5,0.25,0.25}</span></li>
<li>remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li>remove another 225 MW at <span style=\"font-family: monospace;\">t = 1000s </span>from production</li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine govenor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator[1]</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f[1]</span></li>
</ul>
<h4>Questions</h4>
<ol>
<li>Compare the behaviour of the frequency with the one with the droop disabled?</li>
<li>What does the power balance (<code>P_grid_tot</code>) look like compared to no droop within the power grid?</li>
</ol>
</html>"),
        experiment(
          StopTime=2000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end ProductionDroop;

    model RandomLoad
      extends LoadChanges;
                                                                   annotation (preferredView="info",
        Documentation(info="<html>
<h4>Task</h4>
<p>Set up the following:</p>
<ul>
<li>Power set-point: </li>
<li><ul>
<li><span style=\"font-family: monospace;\">pwr_ref</span>: 22.5 MW</li>
</ul></li>
<li>Turbine Governor:</li>
<li><ul>
<li>Droop: <span style=\"font-family: monospace;\">ep = 0.1 </span></li>
</ul></li>
<li>Generator: </li>
<li><ul>
<li>Connection time: connect at t = <span style=\"font-family: monospace;\">150s</span> and stay connected</li>
</ul></li>
<li>Power Grid:</li>
<li><ul>
<li>Disable the droop: <span style=\"font-family: monospace;\">enableDroop = false</span> </li>
<li>default resistive, f and f*f dependent loads: <span style=\"font-family: monospace;\">loadDiv={0.5,0.25,0.25}</span></li>
<li>remove 22.5 MW at <span style=\"font-family: monospace;\">t = 150s</span> from production</li>
<li><b>Add a random load disturbance with a time interval of <span style=\"font-family: monospace;\">h = 10s</span> after 1000 seconds of simulation.</b></li>
</ul></li>
</ul>
<p>Simulate the model for at least 2000 seconds. The long simulation time is necessary in order to let the turbine governor settle down.</p>
<p>The interesting variables are:</p>
<ul>
<li><span style=\"font-family: monospace;\">generator.summary.P_grid_tot</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.P_generator</span></li>
<li><span style=\"font-family: monospace;\">generator.summary.f</span></li>

</ul>

<h4>Questions</h4>
<ol>
<li>Look at the power balance <code>P_grid_tot</code> and try to explain the behaviour!</li>
<li>What happens when you change the interval <code>h</code> for example to <code>h = 1s</code>?</li>
<li>What happens when you enable droop in the Power Grid again? E.g., <code>ep = {0.1,0.08,0.04}</code>?</li>
</ol>
</html>"),
        experiment(
          StopTime=2000,
          Tolerance=1e-05,
          __Dymola_Algorithm="Radau"));
    end RandomLoad;
  end Tutorial10;
  
  package Tutorial11
    record Data "Data record for all simulation parameter"
      extends Modelica.Icons.Record;
      extends OpenHPL.Data;
      // General
      parameter Modelica.SIunits.Power Pn=45e6 "Nominal turbine power"
        annotation(Dialog(group="General"));
      parameter Integer np(min=2)=12
                                    "Number of poles of generator"
        annotation(Dialog(group="General"));
  
      // Mechanical System
      parameter Modelica.SIunits.Inertia J=183e3 "Turbine and generator inertia"
        annotation(Dialog(group="Mechanical System"));
      parameter Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm rpm_n=3000/np*2 "Nominal turbine speed"
        annotation(Dialog(group="Mechanical System", enable=false));
      parameter Modelica.SIunits.AngularVelocity wn = rpm_n / 60 * 2 * Modelica.Constants.pi "Nominal angular velocity"
        annotation(Dialog(group="Mechanical System", enable=false));
        parameter Modelica.SIunits.Time Ta=J*wn^2/Pn "Mechanical time constant"
        annotation(Dialog(group="Mechanical System", enable=false));
  
         // WaterWay
      parameter Modelica.SIunits.VolumeFlowRate Q_n=92 "Nominal flow rate"
        annotation (Dialog(group="WaterWay"));
      parameter Modelica.SIunits.Length H_n=50 "Nominal water head"
        annotation (Dialog(group="WaterWay"));
        parameter Modelica.SIunits.Length L[2]={200,100} "Total length of the waterway"
        annotation (Dialog(group="WaterWay"));
      parameter Modelica.SIunits.Length d[2]= {5.5,5.5} "Average pipe diameter"
        annotation (Dialog(group="WaterWay"));
  
        parameter Modelica.SIunits.Area A[2]=Modelica.Constants.pi/4*d .^ 2 "Average pipe area"
        annotation (Dialog(group="WaterWay", enable=false));
  
      parameter Modelica.SIunits.Time Tw=Q_n/(Modelica.Constants.g_n*H_n)*(L[1]/A[1]+L[2]/A[2])
        "Time constant of the waterway"
        annotation (Dialog(group="WaterWay", enable=false));
  
      annotation (preferredView="info",Documentation(info="<html>
  <p>
  The transfer function of the hydro power system can be represented by:
  </p>
  <img src=\"modelica://FM3217_2020/Resources/Images/TFsystem.png\">
  <p>For the determining the time constants of the transfer-function based model we need to determine the following two time constants:</p>
  <ul>
  <li><img src=\"modelica://FM3217_2020/Resources/Images/equations/equation-Ta.png\" alt=\"T_a = J * w_n^2 / P_n\"/></li>
  <li><img src=\"modelica://FM3217_2020/Resources/Images/equations/equation-Tw.png\" alt=\"T_w = Q_n/(g*H_n) * (L/A)\"/></li>
  </ul>
  <p>We do not worry about the controller parameters since we are going 
  to use a standard PID model from the standard libary.</p>
  <p>See also the list of parameters below.</p>
  </html>"));
    end Data;
  
    model MechanicalTF
  
        Data data                        annotation(Placement(visible = true, transformation(origin={-50,90},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.TransferFunction MechanicalSystem(a = {data.Ta, 0}, b = {1})  annotation(Placement(visible = true, transformation(origin={-10,-20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Add add annotation(Placement(visible = true, transformation(origin = {30, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant one(k = 1)  annotation(Placement(visible = true, transformation(origin={-10,-60},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Gain w(k = data.wn, y(unit="rad/s"))
                                              annotation(Placement(visible = true, transformation(origin = {70, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Trapezoid powerChange(amplitude = 0.2, falling = 50, offset = -0.1, period = 200, rising = 50, startTime = 25, width = 50)  annotation(Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Gain P(k = data.Pn)  annotation(Placement(visible = true, transformation(origin={-30,46},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Division division annotation(Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.Inertia inertia(J = data.J, w(fixed = true, start = data.wn))  annotation(Placement(visible = true, transformation(origin = {58, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Torque torque annotation(Placement(visible = true, transformation(origin = {30, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation(Placement(visible = true, transformation(origin = {80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    equation
    connect(powerChange.y, MechanicalSystem.u) annotation(Line(points={{-79,0},{-60.5,0},{-60.5,-20},{-22,-20}},          color = {0, 0, 127}));
        connect(P.u, powerChange.y) annotation(Line(points={{-42,46},{-60,46},{-60,0},{-79,0}},                    color = {0, 0, 127}));
        connect(division.u1, P.y) annotation(Line(points={{-12,46},{-19,46}},                            color = {0, 0, 127}));
        connect(speedSensor.w, division.u2) annotation(Line(points={{80,19},
                {79.75,19},{79.75,10},{-20,10},{-20,34},{-12,34}},                                                                          color = {0, 0, 127}));
        connect(inertia.flange_b, speedSensor.flange) annotation(Line(points = {{68, 40}, {80, 40}, {80, 40}, {80, 40}}));
    connect(torque.flange, inertia.flange_a) annotation(Line(points = {{40, 40}, {48, 40}}));
    connect(division.y, torque.tau) annotation(Line(points={{11,40},{18,
                40}},                                                              color = {0, 0, 127}));
    connect(add.y, w.u) annotation(Line(points = {{41, -40}, {58, -40}}, color = {0, 0, 127}));
    connect(one.y, add.u2) annotation(Line(points={{1,-60},{11,-60},{11,-46},{18,-46}},          color = {0, 0, 127}));
    connect(MechanicalSystem.y, add.u1) annotation(Line(points={{1,-20},{9.5,-20},{9.5,-34},{18,-34}},          color = {0, 0, 127}));
      annotation (experiment(
          StopTime=600,
          Tolerance=1e-05));
    end MechanicalTF;
  
    model WaterWayTF
      HydroPower.HydroSystems.PipeValve pipeValve(
          endD=data.d,
          L(displayUnit="m") = data.L[1] + data.L[2],
          ZL=data.H_n,
          m_dot_nom=data.Q_n*1000,
          dp_nom=550000)
        annotation (Placement(transformation(extent={{-2,30},{18,50}})));
      HydroPower.SinksAndSources.Fixed_pT Source_pT(paraOption=false)
        annotation (Placement(transformation(extent={{-50,30},{-30,50}})));
      HydroPower.SinksAndSources.Fixed_pT Sink_pT(paraOption=false)
        annotation (Placement(transformation(extent={{70,30},{50,50}})));
      Modelica.Blocks.Sources.Ramp Opening(
        offset=1,
        startTime=50,
        height=-0.9,
        duration=10)
        annotation (Placement(transformation(extent={{-100,-70},{-80,
                  -50}})));
      inner HydroPower.System_HPL system_HPL(
          Q_start=data.Q_n,                  constantTemperature=true,
        pipeRoughness=0,
        steadyState=true)
        annotation (Placement(transformation(extent={{-100,80},{-80,100}})));
      Modelica.Blocks.Continuous.TransferFunction waterWay(
          b={-data.Tw,1},
          a={data.Tw/2,1},
          initType=Modelica.Blocks.Types.Init.InitialOutput,
          y_start=1)
        annotation (Placement(transformation(extent={{-40,-70},{-20,-50}})));
      Modelica.Blocks.Sources.RealExpression Q_valve(y=pipeValve.summary.Q_valve)
        annotation (Placement(transformation(extent={{-58,-12},{0,12}})));
      Modelica.Blocks.Sources.RealExpression dp_valve(y=pipeValve.summary.dp_valve)
        annotation (Placement(transformation(extent={{-58,-32},{0,-10}})));
      Modelica.Blocks.Math.Product Ph
        annotation (Placement(transformation(extent={{20,-20},{40,0}})));
        Modelica.Blocks.Math.Gain Ph_TF(k=data.Pn) annotation (
            Placement(transformation(extent={{20,-70},{40,-50}})));
        Data data(Pn(displayUnit="W"))
                       annotation (Placement(transformation(extent={{
                  -60,80},{-40,100}})));
        Modelica.Blocks.Math.Feedback error annotation (Placement(transformation(extent={{60,-20},{80,0}})));
    equation
      connect(Source_pT.b, pipeValve.a) annotation (Line(
          points={{-29,40},{-3,40}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(Opening.y,pipeValve. ValveCtrl) annotation (Line(
          points={{-79,-60},{-60,-60},{-60,60},{8,60},{8,51}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(waterWay.u,Opening. y) annotation (Line(
          points={{-42,-60},{-79,-60}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pipeValve.b, Sink_pT.b)
        annotation (Line(points={{19,40},{49,40}},         color={0,0,255}));
        connect(waterWay.y, Ph_TF.u) annotation (Line(points={{-19,-60},{18,-60}},
                           color={0,0,127}));
        connect(Q_valve.y, Ph.u1) annotation (Line(points={{2.9,0},{10,0},{10,-4},{18,-4}}, color={0,0,127}));
        connect(dp_valve.y, Ph.u2) annotation (Line(points={{2.9,-21},{9.45,-21},{9.45,-16},{18,-16}}, color={0,0,127}));
        connect(Ph.y, error.u1) annotation (Line(points={{41,-10},{62,-10}}, color={0,0,127}));
        connect(Ph_TF.y, error.u2) annotation (Line(points={{41,-60},{70,-60},{70,-18}}, color={0,0,127}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)),
        experiment(StopTime=300, __Dymola_Algorithm="Radau"));
    end WaterWayTF;
  
    model TFbasedModel
  
      Modelica.Blocks.Continuous.TransferFunction mechanicalSystem(a={data.Ta,
              0})
        annotation (Placement(transformation(extent={{40,40},{60,60}})));
      Modelica.Blocks.Continuous.TransferFunction waterWay(b={-data.Tw,
              1}, a={data.Tw/2,1})
        annotation (Placement(transformation(extent={{-20,40},{0,60}})));
      Modelica.Blocks.Math.Feedback feedback
        annotation (Placement(transformation(extent={{10,40},{30,60}})));
      HydroPower.ControllersAndSensors.TurbineGovernorAnalog turbineGovernorTF(
          enableDroop=false,
        Kd_noLoad=0.05,
        Ki_noLoad=0.025,
        K_noLoad=0.2,
        tRamp=40,
          DeadBand=0,
          P_generator_nom=data.Pn)
                    annotation (Placement(transformation(extent={{-60,40},{-40,
                60}},
              rotation=0)));
      Modelica.Blocks.Sources.Constant zero(k=0)
        annotation (Placement(transformation(extent={{-92,54},{-80,66}})));
      Modelica.Blocks.Sources.BooleanConstant MCBon(k=false)
        annotation (Placement(transformation(extent={{-82,72},{-70,84}})));
      Modelica.Blocks.Sources.RealExpression turbineLosses(y=2.26e6)
        annotation (Placement(transformation(extent={{-60,16},{-30,36}})));
      Modelica.Blocks.Sources.RealExpression gridLoad(y=0e6)
        annotation (Placement(transformation(extent={{-60,4},{-30,24}})));
        Data data(
          Pn(displayUnit="W") = 45e6,
          Q_n=92,
          H_n=50,
          d={5.5,5.5},
          L={200,100},
          J=183000,
          rpm_n=500) annotation (Placement(transformation(extent={{-60,
                  80},{-40,100}})));
        Modelica.Blocks.Math.Add add(k1=1/data.Pn, k2=1/data.Pn)
          annotation (Placement(transformation(extent={{-10,10},{10,30}})));
    equation
  
      connect(waterWay.y, feedback.u1) annotation (Line(
          points={{1,50},{12,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(mechanicalSystem.u, feedback.y) annotation (Line(
          points={{38,50},{29,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(turbineGovernorTF.y, waterWay.u) annotation (Line(
          points={{-39,50},{-22,50}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(mechanicalSystem.y, turbineGovernorTF.f) annotation (Line(
          points={{61,50},{72,50},{72,70},{-66,70},{-66,43},{-61,43}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(zero.y, turbineGovernorTF.P_generator) annotation (Line(
          points={{-79.4,60},{-70,60},{-70,57},{-61,57}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(MCBon.y, turbineGovernorTF.isMCB) annotation (Line(
          points={{-69.4,78},{-50,78},{-50,61}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(turbineGovernorTF.P_reference, zero.y) annotation (Line(
          points={{-56,61},{-56,64},{-70,64},{-70,60},{-79.4,60}},
          color={0,0,127},
          smooth=Smooth.None));
        connect(turbineLosses.y, add.u1) annotation (Line(points={{
                -28.5,26},{-28.5,26},{-12,26}}, color={0,0,127}));
        connect(gridLoad.y, add.u2) annotation (Line(points={{-28.5,14},
                {-28.5,14},{-12,14}}, color={0,0,127}));
        connect(add.y, feedback.u2) annotation (Line(points={{11,20},{
                20,20},{20,42}}, color={0,0,127}));
  
      annotation (experiment(
          StopTime=600,
          Tolerance=1e-05));
    end TFbasedModel;
  
    model CleanExample
        "Hydro plant model connecting to grid at t=150s and disconnect at t=350s"
      extends Modelon.Icons.Experiment;
  
        HydroPower.HydroSystems.Reservoir reservoir1(
          L=500,
          depthIntake={0,15},
          steadyState=false,
          H_start=fill((100), (4)),
          Hmax=fill((105), (4)),
          n=4) annotation (Placement(transformation(extent={{-60,-74},{-40,-54}}, rotation=0)));
        HydroPower.HydroSystems.Reservoir reservoir2(
          L=500,
          nSegmentIntake={1,3},
          depth={70,70,70,70},
          depthIntake={0,1},
          steadyState=false,
          H_start=fill((60), (4)),
          Hmax=fill((70), (4)),
          n=4) annotation (Placement(transformation(extent={{70,-90},{90,-70}}, rotation=0)));
      HydroPower.HydroSystems.Pipe headrace(
        L=200,
        ZL=90,
        ZR=40,
        endL={5,5},
        endD={5.5,5.5},
        Q_start=0.1,
        n=5,
        p_start=540000,
        enable_dataVizPort_lower=false) annotation (Placement(
            transformation(extent={{-30,-80},{-10,-60}}, rotation=0)));
      HydroPower.MechanicalSystems.BasicTurbine turbine(
        np=32,
        H_nom=50,
        tableOnFile=true,
        LdraftTube=10,
        DavDraftTube=2,
        LscrollCase=5,
        DavScrollCasing=2,
        PUInFlowTables=true,
        QTableName="Qtab",
        Q_nom=92,
        H_start=100,
        H_start_draftTube=40,
        Ty=0.4,
        yvLim1=[-0.1, 0.1],
        yvLim2=[-0.2, 0.2],
        TurbineDataFile=Modelica.Utilities.Files.loadResource(HydroPower.TABLE_DIR
             + "TurbineDataFile.mat"),
        P_nom=45000000) annotation (Placement(transformation(extent={{-2,-80},
                  {18,-60}},
                     rotation=0)));
  
      HydroPower.HydroSystems.Pipe tailrace(
        n=4,
        endL={5,5},
        ZL=40,
        ZR=0,
        L=100,
        endD={5.5,5.5},
        Q_start=0.1,
        p_start=400000,
        enable_dataVizPort_lower=false) annotation (Placement(
            transformation(extent={{30,-80},{50,-60}}, rotation=0)));
      HydroPower.ElectricalSystems.PowerGrid powerGrid(
        startTime=1e6,
        unitsJ={122000,5.5e6,8000},
        NoLoadUnits={200,400,1000},
        distNoGen={-2,0,0},
        distTgen={150,1e6,1e6}) annotation (Placement(transformation(
              extent={{-80,-40},{-60,-20}}, rotation=0)));
  
      Modelica.Blocks.Sources.Ramp pwr_ref(
        duration=10,
        height=0,
        offset=45e6,
        startTime=1e6) annotation (Placement(transformation(extent={{-25,-16},
                {-13,-4}},
                      rotation=0)));
      HydroPower.ElectricalSystems.GeneratorAndMCB generator(
        Kdmp={0.05},
        f_start=0,
        timeMCB_open={200},
          J={183000},
          np={12},
          P_nom={45000000},
          timeMCB_close={1e6})
                          annotation (Placement(transformation(extent={{-50,-40},
                {-30,-20}},
                       rotation=0)));
      HydroPower.ControllersAndSensors.TurbineGovernorAnalog turbineGovernor(
        ep=1,
        DeadBand=0.001,
        Ki_load=0.1,
        Kd_load=0.5,
        Kd_noLoad=0.05,
        Ki_noLoad=0.025,
        K_noLoad=0.2,
        K_load=0.4,
        tRamp=40,
        P_generator_nom=generator.P_nom[1],
        enableRamp=false) annotation (Placement(transformation(extent=
               {{-10,-40},{10,-20}}, rotation=0)));
      inner HydroPower.System_HPL system_HPL(
        steadyState=true,
        pipeRoughness=0.1,
        T_start=293) annotation (Placement(transformation(extent={{
                -100,-80},{-80,-60}})));
    equation
  
      connect(powerGrid.f_grid, generator.f_grid) annotation (Line(
          points={{-59,-37},{-51,-37}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(powerGrid.P_grid_balance, generator.P_grid_balance) annotation (Line(
          points={{-59,-23},{-51,-23}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(pwr_ref.y, turbineGovernor.P_reference) annotation (Line(
          points={{-12.4,-10},{-6,-10},{-6,-19}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.f_out[1], turbineGovernor.f) annotation (Line(
          points={{-29,-37},{-11,-37}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(generator.onMCB, powerGrid.MCB) annotation (Line(
          points={{-40,-19},{-40,0},{-70,0},{-70,-19}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.onMCB[1], turbineGovernor.isMCB) annotation (Line(
          points={{-40,-19},{-40,0},{0,0},{0,-19}},
          color={255,0,255},
          smooth=Smooth.None));
      connect(generator.P_out[1], turbineGovernor.P_generator) annotation (Line(
          points={{-29,-23},{-11,-23}},
          color={0,0,127},
          smooth=Smooth.None));
  
      connect(generator.f_out, powerGrid.f) annotation (Line(
          points={{-29,-37},{-20,-37},{-20,-52},{-90,-52},{-90,-37},{
              -81,-37}},
          color={0,0,127},
          smooth=Smooth.None));
      connect(reservoir1.a2_pipe, headrace.a) annotation (Line(
          points={{-39,-70},{-31,-70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(turbine.a, headrace.b) annotation (Line(
          points={{-3,-70},{-9,-70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(turbine.b, tailrace.a) annotation (Line(
          points={{19,-70},{29,-70}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(tailrace.b, reservoir2.a1_pipe) annotation (Line(
          points={{51,-70},{60,-70},{60,-86},{69,-86}},
          color={0,0,255},
          smooth=Smooth.None));
      connect(powerGrid.J_grid, generator.J_grid)
        annotation (Line(points={{-59,-30},{-59,-30},{-51,-30}},
                                                              color={0,0,127}));
      connect(turbineGovernor.y, turbine.yGV)
        annotation (Line(points={{11,-30},{12,-30},{16,-30},{16,-60},{
                16,-59},{14,-59}},                      color={0,0,127}));
      connect(generator.f_out[1], turbine.f_generator) annotation (Line(points={{-29,-37},
                {-29,-37},{-20,-37},{-20,-38},{-20,-36},{-20,-52},{2,
                -52},{2,-59}},                           color={0,0,127}));
      connect(generator.P_turbine[1], turbine.TurbineData[1])
        annotation (Line(points={{-48,-41},{-48,-48},{8,-48},{8,-59.6667}},
                          color={0,0,127}));
      annotation (
         experiment(
          StopTime=600,
          Tolerance=1e-005,
          __Dymola_Algorithm="Radau"));
    end CleanExample;
  
    model CombinedModel
        extends TFbasedModel;
        extends CleanExample;
        annotation (experiment(
            StopTime=600,
            Tolerance=1e-05,
            __Dymola_Algorithm="Radau"));
    end CombinedModel;
  
      package OpenModelica
        extends Modelica.Icons.Package;
        record Data
          extends Modelica.Icons.Record;
          extends OpenHPL.Data;
          parameter Modelica.SIunits.Power Pn = 45e6 "Nominal Power";
          parameter Modelica.SIunits.Frequency fn = 50 "Nominal frequency";
          parameter Integer p = 12 "Number of poles of the generator";
          parameter Modelica.SIunits.AngularVelocity wn = 2*Modelica.Constants.pi * fn /(p/2)
                                                                                             "Nominal speed";
          parameter Modelica.SIunits.Inertia J = 183e3 "Inertia of gen and turbine";
          parameter Modelica.SIunits.VolumeFlowRate Qn = 92
                                                           "Nominal flow rate";
          parameter Modelica.SIunits.Height Hn = 50 "Nominal head";
          parameter Modelica.SIunits.Length L[2] = {200,100} "Lenght of the water way";
          parameter Modelica.SIunits.Area A[2] = Modelica.Constants.pi/4*d.^ 2 "Area of pipe segments";
          parameter Modelica.SIunits.Length d[2] = {5.5,5.5} "Pipe diameter";
          parameter Modelica.SIunits.Time Ta = J * wn^2/Pn "Mechanical time constant";
          parameter Modelica.SIunits.Time Tw = Qn/ (Modelica.Constants.g_n *Hn) * (L[1]/A[1] + L[2]/A[2]) "Waterway time constant";
  
        end Data;
  
        model MechanicalTF
          Modelica.Blocks.Continuous.TransferFunction mechanical(a={data.Ta,0}) annotation (
            Placement(visible = true, transformation(origin = {0, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Sources.Trapezoid powerChange(
            amplitude=0.2,
            falling=50,
            offset=-0.1,
            period=200,
            rising=50,
            startTime=25,
            width=50) annotation (
            Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Sources.Constant one(k=1) annotation (
            Placement(visible = true, transformation(origin = {0, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Add add annotation (
            Placement(visible = true, transformation(origin = {40, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Gain gain(k=data.wn) annotation (
            Placement(visible = true, transformation(origin = {70, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Division division annotation (
            Placement(visible = true, transformation(origin = {0, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Blocks.Math.Gain P(k=data.Pn) annotation (
            Placement(visible = true, transformation(origin = {-30, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (
            Placement(visible = true, transformation(origin = {80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
          Modelica.Mechanics.Rotational.Sources.Torque torque annotation (
            Placement(visible = true, transformation(origin = {30, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          Modelica.Mechanics.Rotational.Components.Inertia inertia(J=data.J, w(fixed=true, start=data.wn)) annotation (
            Placement(visible = true, transformation(origin = {58, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  FM3217_2021.Tutorial11.OpenModelica.Data data annotation(
          Placement(visible = true, transformation(origin = {-50, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
          connect(powerChange.y, mechanical.u) annotation (
            Line(points = {{-78, 0}, {-60, 0}, {-60, -52}, {-12, -52}, {-12, -52}}, color = {0, 0, 127}));
          connect(add.u1, mechanical.y) annotation (
            Line(points = {{28, -62}, {20, -62}, {20, -52}, {12, -52}}, color = {0, 0, 127}));
          connect(one.y, add.u2) annotation (
            Line(points = {{12, -80}, {20, -80}, {20, -74}, {28, -74}}, color = {0, 0, 127}));
          connect(gain.u, add.y) annotation (
            Line(points = {{58, -68}, {52, -68}, {52, -68}, {52, -68}}, color = {0, 0, 127}));
          connect(division.y, torque.tau) annotation (
            Line(points = {{11, 40}, {18, 40}}, color = {0, 0, 127}));
          connect(speedSensor.w, division.u2) annotation (
            Line(points = {{80, 19}, {79.75, 19}, {79.75, 10}, {-20, 10}, {-20, 34}, {-12, 34}}, color = {0, 0, 127}));
          connect(torque.flange, inertia.flange_a) annotation (
            Line(points = {{40, 40}, {48, 40}}));
          connect(division.u1, P.y) annotation (
            Line(points = {{-12, 46}, {-19, 46}}, color = {0, 0, 127}));
          connect(inertia.flange_b, speedSensor.flange) annotation (
            Line(points = {{68, 40}, {80, 40}, {80, 40}, {80, 40}}));
          connect(
            P.u, powerChange.y) annotation (
            Line(points = {{-42, 46}, {-60, 46}, {-60, 0}, {-78, 0}, {-78, 0}}, color = {0, 0, 127}));
          annotation (
            Icon(coordinateSystem(grid = {2, 0})));
        end MechanicalTF;
  
        model WaterWayTF
          Modelica.Blocks.Math.Gain Ph_TF(k = data.Pn) annotation (
            Placement(visible = true, transformation(extent = {{20, -70}, {40, -50}}, rotation = 0)));
          Modelica.Blocks.Continuous.TransferFunction waterWay(a = {data.Tw / 2, 1}, b = {-data.Tw, 1}, initType = Modelica.Blocks.Types.Init.InitialOutput, y_start = 1) annotation (
            Placement(visible = true, transformation(extent = {{-40, -70}, {-20, -50}}, rotation = 0)));
          Modelica.Blocks.Sources.Ramp Opening(duration = 10, height = -0.9, offset = 1, startTime = 50) annotation (
            Placement(visible = true, transformation(extent = {{-100, -70}, {-80, -50}}, rotation = 0)));
          OpenHPL.Waterway.Reservoir reservoir(H_r=10) annotation (
            Placement(visible = true, transformation(origin = {-62, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          OpenHPL.Waterway.Reservoir reservoir1(H_r=10) annotation (
            Placement(visible = true, transformation(origin = {70, 20}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
          OpenHPL.Waterway.Pipe pipe(
            D_i=data.d[1],
            H=data.Hn - 10,
            L=data.L[1]) annotation (
            Placement(visible = true, transformation(origin = {-28, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          OpenHPL.Waterway.Pipe pipe1(
            D_i=data.d[2],
            H=10,
            L=data.L[2]) annotation (
            Placement(visible = true, transformation(origin = {22, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          OpenHPL.ElectroMech.Turbines.Turbine2 turbine2(
            ConstEfficiency=true,
            H_n=data.Hn,
            Jt=data.J,
            ValveCapacity=false,
            Vdot_n=data.Qn, enableP_out = true,
            eta_h=1,
            p=data.p) annotation (
            Placement(visible = true, transformation(origin = {-2, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          inner FM3217_2020.Tutorial11.OpenModelica.Data data(
            Steady=true,
            V_0=data.Qn,
            gamma_air=1) annotation (
            Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
          connect(waterWay.u, Opening.y) annotation (
            Line(points = {{-42, -60}, {-79, -60}}, color = {0, 0, 127}));
          connect(waterWay.y, Ph_TF.u) annotation (
            Line(points = {{-19, -60}, {18, -60}}, color = {0, 0, 127}));
          connect(
            pipe.i, reservoir.o) annotation (
            Line(points = {{-38, 34}, {-45, 34}, {-45, 50}, {-52, 50}}, color = {28, 108, 200}));
          connect(
            pipe1.o, reservoir1.o) annotation (
            Line(points = {{32, 32}, {50, 32}, {50, 20}, {60, 20}}, color = {28, 108, 200}));
          connect(
            pipe.o, turbine2.i) annotation (
            Line(points = {{-18, 34}, {-12, 34}}, color = {28, 108, 200}));
          connect(
            turbine2.o, pipe1.i) annotation (
            Line(points = {{8, 34}, {11, 34}, {11, 32}, {12, 32}}, color = {28, 108, 200}));
          connect(
            turbine2.u_t, Opening.y) annotation (
            Line(points = {{-10, 46}, {-8, 46}, {-8, 70}, {-90, 70}, {-90, -18}, {-62, -18}, {-62, -60}, {-78, -60}, {-78, -60}}, color = {0, 0, 127}));
          annotation (
            Icon(coordinateSystem(grid = {2, 0})));
        end WaterWayTF;
  
        model TFbasedModel
  
          Modelica.Blocks.Continuous.TransferFunction mechanicalSystem(a={data.Ta,
                  0})
            annotation (Placement(visible = true, transformation(extent = {{50, -50}, {70, -30}}, rotation = 0)));
          Modelica.Blocks.Continuous.TransferFunction waterWay(b={-data.Tw,
                  1}, a={data.Tw/2,1})
            annotation (Placement(visible = true, transformation(extent = {{-10, -50}, {10, -30}}, rotation = 0)));
          Modelica.Blocks.Math.Feedback feedback
            annotation (Placement(visible = true, transformation(extent = {{20, -50}, {40, -30}}, rotation = 0)));
          Modelica.Blocks.Sources.RealExpression turbineLosses(y=2.26e6)
            annotation (Placement(visible = true, transformation(extent = {{-50, -74}, {-20, -54}}, rotation = 0)));
            Modelica.Blocks.Math.Add add(k1=1/data.Pn, k2=1/data.Pn)
              annotation (Placement(visible = true, transformation(extent = {{0, -80}, {20, -60}}, rotation = 0)));
        Modelica.Blocks.Continuous.PI piTF(T=10) annotation (Placement(visible=true, transformation(
                  origin={-40,-40},
                  extent={{-10,-10},{10,10}},
                  rotation=0)));
            Modelica.Blocks.Math.Gain gain(k=-1) annotation (Placement(visible = true, transformation(extent = {{-72, -44}, {-64, -36}}, rotation = 0)));
        Modelica.Blocks.Sources.Pulse gridLoad(amplitude = 10e6, period = 400, startTime = 200)  annotation (
              Placement(visible = true, transformation(origin = {-30, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
          inner FM3217_2020.Tutorial11.OpenModelica.Data data(Steady=true, V_0=data.Qn) annotation (
            Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        equation
            connect(waterWay.y, feedback.u1) annotation (
              Line(points = {{11, -40}, {22, -40}}, color = {0, 0, 127}));
            connect(mechanicalSystem.u, feedback.y) annotation (
              Line(points = {{48, -40}, {39, -40}}, color = {0, 0, 127}));
            connect(turbineLosses.y, add.u1) annotation (
              Line(points = {{-18.5, -64}, {-18.5, -64}, {-2, -64}}, color = {0, 0, 127}));
            connect(add.y, feedback.u2) annotation (
              Line(points = {{21, -70}, {30, -70}, {30, -48}}, color = {0, 0, 127}));
            connect(piTF.y, waterWay.u) annotation (
              Line(points = {{-29, -40}, {-12, -40}, {-12, -40}, {-12, -40}}, color = {0, 0, 127}));
            connect(mechanicalSystem.y, gain.u) annotation (
              Line(points = {{71, -40}, {80, -40}, {80, -20}, {-80, -20}, {-80, -40}, {-72.8, -40}}, color = {0, 0, 127}));
            connect(gain.y, piTF.u) annotation (
              Line(points = {{-63.6, -40}, {-52, -40}}, color = {0, 0, 127}));
        connect(gridLoad.y, add.u2) annotation (
              Line(points={{-19,-86},{-14,-86},{-14,-76},{-2,-76}},          color = {0, 0, 127}));
            annotation (experiment(
              StopTime=600,
              Tolerance=1e-05));
        end TFbasedModel;
  
        model CombinedModel
            extends TFbasedModel(data(
                                 Steady =      false));
        OpenHPL.Waterway.Pipe pipe1(D_i = data.d[2], H = 10, L = data.L[2]) annotation (
              Placement(visible = true, transformation(origin={30,20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OpenHPL.ElectroMech.Turbines.Turbine2 turbine2(
            ConstEfficiency=true,
            H_n=data.Hn, Jt = data.J / 2, Ploss = 1.13e6,
            ValveCapacity=false,
            Vdot_n=92, eta_h = 1, p = data.p) annotation (Placement(visible=true, transformation(
                origin={0,20},
                extent={{-10,-10},{10,10}},
                rotation=0)));
        OpenHPL.Waterway.Reservoir reservoir1(H_r = 10) annotation (
              Placement(visible = true, transformation(origin={62,20},    extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        OpenHPL.Waterway.Reservoir reservoir(H_r = 10) annotation (
              Placement(visible = true, transformation(origin={-70,20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OpenHPL.Waterway.Pipe pipe(D_i = data.d[1], H = data.Hn - 10, L = data.L[1], vertical = false) annotation (
              Placement(visible = true, transformation(origin={-30,20},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
        OpenHPL.ElectroMech.Generators.SimpleGen2 simpleGen2(Jg = data.J / 2, Ploss = 1.13e6,
            enable_f=true, enable_w = false)                                                                                          annotation (
              Placement(visible = true, transformation(origin={0,50},    extent = {{10, -10}, {-10, 10}}, rotation = 0)));
        Modelica.Blocks.Continuous.PI piReal(T = 10, k = 1) annotation (Placement(visible=true, transformation(
                  origin={-30,68},
                  extent={{-10,-10},{10,10}},
                  rotation=0)));
        Modelica.Blocks.Math.Feedback feedback2 annotation (
              Placement(visible = true, transformation(extent={{-74,58},{-54,78}},      rotation = 0)));
        Modelica.Blocks.Sources.Constant sysF(k = 50)  annotation (
              Placement(visible = true, transformation(origin={-92,68},    extent = {{-6, -6}, {6, 6}}, rotation = 0)));
        Modelica.Blocks.Math.Gain df_pu(k = 1 / 50) annotation (
              Placement(visible = true, transformation(extent={{-36, 86},{-28, 94}},      rotation = 0)));
        equation
          connect(reservoir1.o, pipe1.o) annotation (
            Line(points = {{52, 20}, {40, 20}}, color = {28, 108, 200}));
          connect(reservoir.o, pipe.i) annotation (
            Line(points = {{-60, 20}, {-40, 20}}, color = {28, 108, 200}));
          connect(pipe.o, turbine2.i) annotation (
            Line(points = {{-20, 20}, {-10, 20}}, color = {28, 108, 200}));
          connect(turbine2.o, pipe1.i) annotation (
            Line(points = {{10, 20}, {20, 20}}, color = {28, 108, 200}));
          connect(sysF.y, feedback2.u1) annotation (
            Line(points = {{-85, 68}, {-72, 68}}, color = {0, 0, 127}));
          connect(simpleGen2.f, feedback2.u2) annotation (
            Line(points = {{-11, 46}, {-64, 46}, {-64, 60}}, color = {0, 0, 127}));
          connect(feedback2.y, piReal.u) annotation (
            Line(points = {{-55, 68}, {-42, 68}}, color = {0, 0, 127}));
          connect(piReal.y, turbine2.u_t) annotation (
            Line(points = {{-19, 68}, {-8, 68}, {-8, 32}}, color = {0, 0, 127}));
          connect(df_pu.u, feedback2.y) annotation (
            Line(points = {{-37, 90}, {-46, 90}, {-46, 68}, {-55, 68}}, color = {0, 0, 127}));
          connect(turbine2.flange, simpleGen2.flange) annotation (
            Line(points = {{0, 29.8}, {0, 40.2}}, color = {0, 0, 0}));
          connect(
            simpleGen2.Pload, gridLoad.y) annotation (
            Line(points = {{0, 62}, {0, 62}, {0, 80}, {90, 80}, {90, -86}, {-18, -86}, {-18, -86}}, color = {0, 0, 127}));
          annotation (experiment(
                StopTime=600));
        end CombinedModel;
        annotation (
          Icon(coordinateSystem(grid = {2, 0})));
      end OpenModelica;
  
  end Tutorial11;
  annotation (uses(Modelica(version="3.2.3"), HydroPower(version="2.13"),
      Modelon(version="3.7")));
end FM3217_2021;
>>>>>>> Stashed changes
