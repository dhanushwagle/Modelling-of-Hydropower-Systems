within ;
package DroopExample

  model System
    extends Modelica.Icons.Example;
    PowerPlant P1(
      f_n=50,
      R=0.004,
      P_set(
        displayUnit="MW",
        start=100e6) = 50000000,
      P_n(displayUnit="MW") = 100000000)
      annotation (Placement(transformation(extent={{-80,-50},{-40,-10}})));
    Grid grid(fend=50.1) annotation (Placement(transformation(
          rotation=0,
          extent={{-20,20},{20,-20}},
          origin={0,-70})));
    PowerPlant P2(
      f_n=50,
      R=0.008,
      P_set(
        displayUnit="MW",
        start=100e6) = 50000000,
      P_n(displayUnit="MW") = 100000000)
                                     "Power plant"
      annotation (Placement(transformation(extent={{80,-50},{40,-10}})));
  //
    Real a1,a2,b1,b2,y1_100,y2_100,fsys "Linear factors for the droop line equations";

    Modelica.Blocks.Sources.Trapezoid disturbance(
      startTime=grid.startTime,
      offset=grid.fstart,
      amplitude=grid.fend - grid.fstart,
      rising=grid.duration,
      width=grid.duration,
      falling=grid.duration,
      period=6*grid.duration) annotation (Placement(transformation(
            extent={{-70,-80},{-50,-60}})));
  equation

    a1 = -P1.R*100;
    a2 = -P2.R*100;
    b1 = P1.f_n * (1 + (P1.R * P1.P_set/P1.P_n)*200);
    b2 = P2.f_n * (1 + (P2.R * P2.P_set/P2.P_n)*200);
    y1_100 = a1 * 100 + b1;
    y2_100 = a2 * 100 + b2;
    fsys = 50-(50-grid.fsys)*200;
    connect(P1.flange, grid.flange) annotation (Line(points={{-40,-30},{
            0,-30},{0,-50}},
                         color={0,0,0}));
    connect(P2.flange, grid.flange)
      annotation (Line(points={{40,-30},{0,-30},{0,-50}},
                                                       color={0,0,0}));
    connect(disturbance.y, grid.f_sys) annotation (Line(points={{-49,-70},{-24,-70}}, color={0,0,127}));
    annotation (Diagram(coordinateSystem(grid={1,1}, initialScale=0.1),
                        graphics={
          Line(
            points={{0,100},{0,0}},
            color={0,0,0},
            arrow={Arrow.Filled,Arrow.None}),
          Line(
            points={{100,0},{-100,0}},
            color={0,0,0},
            arrow={Arrow.Filled,Arrow.Filled}),
          Line(points={{-1,50},{1,50}},   color={0,0,0}),
          Line(points={{-1,60},{1,60}},   color={0,0,0}),
          Line(points={{-1,70},{1,70}},   color={0,0,0}),
          Line(points={{-1,80},{1,80}},   color={0,0,0}),
          Line(points={{-1,90},{1,90}},   color={0,0,0}),
          Line(points={{-1,40},{1,40}},   color={0,0,0}),
          Line(points={{-1,30},{1,30}},   color={0,0,0}),
          Line(points={{-1,20},{1,20}},   color={0,0,0}),
          Line(points={{-1,10},{1,10}},   color={0,0,0}),
          Text(
            extent={{-8,48},{-2,52}},
            lineColor={0,140,72},
            horizontalAlignment=TextAlignment.Right,
            textString="50.00",
            textStyle={TextStyle.Bold}),
          Text(
            extent={{-8,58},{-2,62}},
            lineColor={0,140,72},
            horizontalAlignment=TextAlignment.Right,
            textString="50.05"),
          Text(
            extent={{-8,68},{-2,72}},
            lineColor={244,125,35},
            horizontalAlignment=TextAlignment.Right,
            textString="50.10"),
          Text(
            extent={{-8,78},{-2,82}},
            lineColor={238,46,47},
            horizontalAlignment=TextAlignment.Right,
            textString="50.15"),
          Text(
            extent={{-8,88},{-2,92}},
            lineColor={162,29,33},
            horizontalAlignment=TextAlignment.Right,
            textString="50.20"),
          Text(
            extent={{-8,38},{-2,42}},
            lineColor={0,140,72},
            horizontalAlignment=TextAlignment.Right,
            textString="49.95"),
          Text(
            extent={{-8,28},{-2,32}},
            lineColor={244,125,35},
            horizontalAlignment=TextAlignment.Right,
            textString="49.90"),
          Text(
            extent={{-8,18},{-2,22}},
            lineColor={238,46,47},
            horizontalAlignment=TextAlignment.Right,
            textString="49.85"),
          Text(
            extent={{-8,8},{-2,12}},
            lineColor={162,29,33},
            horizontalAlignment=TextAlignment.Right,
            textString="49.80"),
          Text(
            extent={{-4,-3},{4,-7}},
            lineColor={0,0,0},
            textString="0"),
          Line(points={{0,0},{0,0}},       color={0,0,0}),
          Line(points={{0,1},{0,-1},{0,-1}},       color={0,0,0}),
          Line(points={{-50,1},{-50,-1},{-50,-1}},       color={0,0,0}),
          Line(points={{50,1},{50,-1},{50,-1}},       color={0,0,0}),
          Line(points={{10,1},{10,-1},{10,-1}},       color={0,0,0}),
          Line(points={{20,1},{20,-1},{20,-1}},       color={0,0,0}),
          Line(points={{30,1},{30,-1},{30,-1}},       color={0,0,0}),
          Line(points={{40,1},{40,-1},{40,-1}},       color={0,0,0}),
          Line(points={{60,1},{60,-1},{60,-1}},       color={0,0,0}),
          Line(points={{70,1},{70,-1},{70,-1}},       color={0,0,0}),
          Line(points={{80,1},{80,-1},{80,-1}},       color={0,0,0}),
          Line(points={{90,1},{90,-1},{90,-1}},       color={0,0,0}),
          Line(points={{-10,1},{-10,-1},{-10,-1}},       color={0,0,0}),
          Line(points={{-20,1},{-20,-1},{-20,-1}},       color={0,0,0}),
          Line(points={{-30,1},{-30,-1},{-30,-1}},       color={0,0,0}),
          Line(points={{-40,1},{-40,-1},{-40,-1}},       color={0,0,0}),
          Line(points={{-60,1},{-60,-1},{-60,-1}},       color={0,0,0}),
          Line(points={{-70,1},{-70,-1},{-70,-1}},       color={0,0,0}),
          Line(points={{-80,1},{-80,-1},{-80,-1}},       color={0,0,0}),
          Line(points={{-90,1},{-90,-1},{-90,-1}},       color={0,0,0}),
          Text(
            extent={{46,-3},{54,-7}},
            lineColor={0,0,0},
            textString="50"),
          Text(
            extent={{-54,-3},{-46,-7}},
            lineColor={0,0,0},
            textString="50"),
          Text(
            extent={{-100,0},{-90,-10}},
            lineColor={0,0,0},
            textString="P1
[MW]"),   Text(
            extent={{90,0},{100,-10}},
            lineColor={0,0,0},
            textString="P2
[MW]"),    Line(
            points=DynamicSelect({{50,0},{50,50}},
            {{P2.powerSensor.power/1e6,0},{P2.powerSensor.power/1e6,fsys}}),
            color={0,0,0},
            thickness=0.5,
            arrow={Arrow.None,Arrow.Half}),
          Line(points=DynamicSelect({{0,70},{-100,30}},
          {{0,b1},{-100,y1_100}}),
                color={0,0,0}),
          Line(points=DynamicSelect({{0,60},{100,40}},
          {{0,b2},{100,y2_100}}),
                color={0,0,0}),
          Line(
            points=DynamicSelect({{-50,0},{-50,50}},
            {{-P1.powerSensor.power/1e6,0},{-P1.powerSensor.power/1e6,fsys}}),
            color={0,0,0},
            thickness=0.5,
            arrow={Arrow.None,Arrow.Half}),
         Line(points=DynamicSelect({{-100,50},{100,50}},
             {{-100, fsys},{100,fsys}}),
             color={0,140,72},
             pattern=LinePattern.Dash)}),
              Icon(coordinateSystem(grid={1,1},
            initialScale=0.1)),
      experiment(StopTime=250));
  end System;

  model Grid
    parameter Modelica.SIunits.Frequency fstart=50 "Start frequency";
    parameter Modelica.SIunits.Frequency fend=49.9 "End frequency";
    parameter Modelica.SIunits.Time startTime=50 "Start time of frequency change";
    parameter Modelica.SIunits.Time duration=50 "Duration time of frequency change";
    parameter Modelica.SIunits.Inertia Jgrid=5000 "Equivalent Inertia of grid";
    Modelica.SIunits.Frequency fsys = fromHz.u "System frequency";
    Modelica.Mechanics.Rotational.Components.Inertia Grid(J=Jgrid)
      annotation (Placement(transformation(extent={{60,-10},{40,10}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed annotation (
        Placement(transformation(extent={{-10,-10},{10,10}})));
    Modelica.Blocks.Math.Gain fromHz(k=2*Modelica.Constants.pi)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=0,
          origin={-50,0})));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange
      annotation (Placement(transformation(rotation=0, extent={{-10,-110},
              {10,-90}})));
    Modelica.Blocks.Interfaces.RealInput f_sys
      "Input signal connector" annotation (Placement(transformation(
            extent={{-140,-20},{-100,20}})));
  equation
    connect(speed.flange, Grid.flange_b)
      annotation (Line(points={{10,0},{40,0}},   color={0,0,0}));
    connect(fromHz.y, speed.w_ref)
      annotation (Line(points={{-39,0},{-12,0}},
                                               color={0,0,127}));
    connect(flange, Grid.flange_a)
      annotation (Line(points={{0,-100},{80,-100},{80,0},{60,0}},
                                                  color={0,0,0}));
    connect(fromHz.u, f_sys)
      annotation (Line(points={{-62,0},{-120,0}}, color={0,0,127}));
    annotation (Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-80,80},{80,40}},
            lineColor={0,140,72},
            textString="fstart = %fstart"),
          Text(
            extent={{-80,-40},{80,-80}},
            lineColor={0,140,72},
            textString="fend = %fend"),
                               Text(
            extent={{-80,20},{80,-20}},
            lineColor={0,140,72},
            horizontalAlignment=TextAlignment.Left,
            textString="fsys = "+DynamicSelect("%fstart", String(fsys, significantDigits=4))),
          Text(
            extent={{-100,140},{100,100}},
            lineColor={0,0,255},
            textString="%name")}));
  end Grid;

  model PowerPlant
    parameter Modelica.SIunits.Power P_set(start=50e6)
      "Power setpoint of plant";
    parameter Modelica.SIunits.Power P_n(start=100e6) "Nominal power of plant";
    parameter Modelica.SIunits.Frequency f_n(start=50) "Nominal grid frequency";
    parameter Real R(start=0.1) "Droop setting for the plant governor";

    Modelica.Mechanics.Rotational.Sources.Torque torque
      annotation (Placement(transformation(extent={{40,-10},{60,10}})));
    Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={70,-30})));
    Modelica.Blocks.Math.Gain toHz(k=1/(Modelica.Constants.pi*2))
      annotation (Placement(transformation(
          extent={{10,-10},{-10,10}},
          rotation=0,
          origin={-48,-60})));
    Modelica.Blocks.Math.Division division
      annotation (Placement(transformation(extent={{10,-10},{30,10}})));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange
      annotation (Placement(transformation(extent={{110,-10},{90,10}})));
    Modelica.Blocks.Math.Division division1
      annotation (Placement(transformation(extent={{0,66},{20,86}})));
    Modelica.Blocks.Sources.RealExpression fn(y=f_n) annotation (
        Placement(transformation(extent={{-80,60},{-60,80}})));
    Modelica.Blocks.Math.Division division3
      annotation (Placement(transformation(extent={{40,50},{60,70}})));
    Modelica.Blocks.Math.Division division2
      annotation (Placement(transformation(extent={{0,30},{20,50}})));
    Modelica.Blocks.Sources.RealExpression Pn(y=P_n) annotation (
        Placement(transformation(extent={{-80,20},{-60,40}})));
    Modelica.Blocks.Sources.RealExpression Droop(y=R) annotation (
        Placement(transformation(extent={{-80,40},{-60,60}})));
    Modelica.Blocks.Math.Feedback Pstar annotation (Placement(
          transformation(extent={{-50,10},{-30,-10}})));
    Modelica.Blocks.Sources.RealExpression Pset(y=P_set) annotation (
        Placement(transformation(extent={{-80,-10},{-60,10}})));
    Modelica.Mechanics.Rotational.Sensors.PowerSensor powerSensor
      annotation (Placement(transformation(extent={{70,-10},{90,10}})));
    Modelica.Blocks.Math.Feedback deltaF annotation (Placement(
          transformation(extent={{-40,72},{-20,92}})));
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax=P_n, uMin=0)
      annotation (Placement(transformation(extent={{-22,-6},{-10,6}})));
  equation
    connect(torque.tau, division.y) annotation (Line(
        points={{38,0},{31,0}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(fn.y, division1.u2)
      annotation (Line(points={{-59,70},{-2,70}}, color={0,0,127}));
    connect(division1.y, division3.u1) annotation (Line(points={{21,76},
            {28,76},{28,66},{38,66}}, color={0,0,127}));
    connect(Pn.y, division2.u2) annotation (Line(points={{-59,30},{-20,30},
            {-20,34},{-2,34}}, color={0,0,127}));
    connect(Droop.y, division2.u1) annotation (Line(points={{-59,50},{-20,
            50},{-20,46},{-2,46}}, color={0,0,127}));
    connect(division2.y, division3.u2) annotation (Line(points={{21,40},
            {28,40},{28,54},{38,54}}, color={0,0,127}));
    connect(Pstar.u1, Pset.y)
      annotation (Line(points={{-48,0},{-59,0}}, color={0,0,127}));
    connect(torque.flange, powerSensor.flange_a)
      annotation (Line(points={{60,0},{70,0}}, color={0,0,0}));
    connect(flange, powerSensor.flange_b)
      annotation (Line(points={{100,0},{90,0}}, color={0,0,0}));
    connect(speedSensor.flange, torque.flange) annotation (Line(points={
            {70,-20},{70,0},{60,0}}, color={0,0,0}));
    connect(division3.y, Pstar.u2) annotation (Line(points={{61,60},{
            70,60},{70,20},{-40,20},{-40,8}},
                                           color={0,0,127}));
    connect(speedSensor.w, toHz.u) annotation (Line(points={{70,-41},{70,
            -60},{-36,-60}}, color={0,0,127}));
    connect(division.u2, toHz.u) annotation (Line(points={{8,-6},{0,
            -6},{0,-60},{-36,-60}},
                                  color={0,0,127}));
    connect(division1.u1, deltaF.y)
      annotation (Line(points={{-2,82},{-21,82}}, color={0,0,127}));
    connect(deltaF.u2, division1.u2) annotation (Line(points={{-30,74},{
            -30,70},{-2,70}}, color={0,0,127}));
    connect(deltaF.u1, toHz.y) annotation (Line(points={{-38,82},{-90,82},
            {-90,-60},{-59,-60}}, color={0,0,127}));
    connect(Pstar.y, limiter.u)
      annotation (Line(points={{-31,0},{-23.2,0}}, color={0,0,127}));
    connect(limiter.y, division.u1) annotation (Line(points={{-9.4,0},
            {0,0},{0,6},{8,6}}, color={0,0,127}));
    annotation (Icon(graphics={
          Rectangle(extent={{-100,100},{100,-100}}, lineColor={28,108,200},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
                               Text(
            extent={{-80,100},{80,70}},
            lineColor={238,46,47},
            horizontalAlignment=TextAlignment.Left,
            textString="P_set = %P_set"),
                               Text(
            extent={{-80,60},{84,30}},
            lineColor={238,46,47},
            horizontalAlignment=TextAlignment.Left,
            textString="P_n = %P_n"),
                               Text(
            extent={{-80,-20},{82,-50}},
            lineColor={0,140,72},
            horizontalAlignment=TextAlignment.Left,
            textString="f_n = %f_n"),
                               Text(
            extent={{-80,-60},{80,-90}},
            lineColor={244,125,35},
            horizontalAlignment=TextAlignment.Left,
            textString="R = %R "),
                               Text(
            extent={{-80,20},{80,-20}},
            lineColor={238,46,47},
            horizontalAlignment=TextAlignment.Left,
            textString="P* = "+DynamicSelect("%P_set", String(powerSensor.power, significantDigits=4))),
          Text(
            extent={{-100,140},{100,100}},
            lineColor={0,0,255},
            textString="%name")}));
  end PowerPlant;
  annotation (preferredView="info",uses(Modelica(version="3.2.3")),
    Documentation(info="<html>
<p>Droop definition: </p>
<blockquote><pre>
      &Delta;f / f<sub>n</sub>
R = - -------
      &Delta;P / P<sub>n</sub>
</pre></blockquote>
</html>"));
end DroopExample;
