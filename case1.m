function mpc = case1
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	0	110	1	1.1	0.9
	2	1	50	20	0	0	1	1	0	110	1	1.1	0.9
	3	1	50	20	0	0	1	1	0	110	1	1.1	0.9
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 0	 0	99	-99	115/110	100	1	99	0	0	0	0	0	0	0	0	0	0	0	0
	2	20	-5	99	-99	      1	100	1	99	0	0	0	0	0	0	0	0	0	0	0	0
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
R = 12/121;
X = 40.9/121;
B = 277e-6*121;
mpc.branch = [
	1	2	R	X	B	0	0	0	0	0	1	-360	360
    1	3	R	X	B	0	0	0	0	0	1	-360	360
	2	3	R	X	B	0	0	0	0	0	1	-360	360
];

%% bus location
%	bus_i	x	y
mpc.buslocation = [
    1	0	0
    2	112	0
    3	56	-74
];
