function results = my_runpf(casedata)
%MY_RUNPF  Wrapper function around MATPOWER's runpf.
%          Runs a power flow with the full AC Newton's method.
%
%   results = my_runpf(casedata)
%
%   Input:
%       casedata : either a MATPOWER case struct or a string containing
%           the name of the file with the case data
%           (see also CASEFORMAT and LOADCASE)
%   Output:
%       results : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%
%   See also RUNPF, CASEFORMAT and LOADCASE.

% Solve the power flow
define_constants;
opt = mpoption('verbose', 2, 'out.all', 0);

mpc = loadcase(casedata);
fprintf('Matlab: casedata = %s\n', casedata);
if isfield(mpc, 'run')
    if mpc.run
        results = runpf(mpc, opt);
    else
        results = ext2int(mpc);
        results.success = 1;
        fprintf('Results read from file.\n');
    end
else
    results = runpf(mpc, opt);
end
e2i = results.order.bus.e2i; % external to internal bus numbering
results = rmfield(results, {'order'});
results.branch(:, [PT QT]) = -results.branch(:, [PT QT]);

% Delete branches with zero flows at both ends
eps = 1e-4;
i1 = abs(results.branch(:, PF)) < eps;
i2 = abs(results.branch(:, PT)) < eps;
i3 = abs(results.branch(:, QF)) < eps;
i4 = abs(results.branch(:, QT)) < eps;
i = i1 & i2 & i3 & i4;
if any(i)
    disp('Branches with zeros flows at both ends (PF, PT, QF, QT)');
    disp(results.branch(i, [F_BUS T_BUS]));
    results.branch(i,:) = [];
end

% Nullify flows below given treshold eps
for col = [PF PT QF QT]
    i = abs(results.branch(:, col)) < eps;
    results.branch(i, col) = 0;
    switch col
        case PF
            i = results.branch(:, col) ~= 0;
            fprintf('min(PF) = %.6f MW\n', min(abs(results.branch(i,col))));
        case PT
            i = results.branch(:, col) ~= 0;
            fprintf('min(PT) = %.6f MW\n', min(abs(results.branch(i,col))));
        case QF
            i = results.branch(:, col) ~= 0;
            fprintf('min(QF) = %.6f Mvar\n', min(abs(results.branch(i,col))));
        case QT
            i = results.branch(:, col) ~= 0;
            fprintf('min(QT) = %.6f Mvar\n', min(abs(results.branch(i,col))));
    end
end

% PS is the active power demanded by the shunt
results.PS = results.bus(:, GS).*results.bus(:, VM).^2;
% QS is the reactive power injected by the shunt
results.QS = results.bus(:, BS).*results.bus(:, VM).^2;

% Calculate reactive power flows in branches' shunt elements
% The flows are oriented towards "from" and "to" buses at both branch ends
branch = results.branch;
stat = branch(:, BR_STATUS); % ones at in-service branches
nl = size(stat, 1);
tap = ones(nl, 1);        % default tap ratio = 1
i = find(branch(:, TAP)); % indices of non-zero tap ratios
tap(i) = branch(i, TAP);  % assign non-zero tap ratios
V = results.bus(:, VM).*exp(1j*results.bus(:,VA)/180*pi); % complex bus voltages
f = e2i(branch(:, F_BUS)); % vector of "from" buses
t = e2i(branch(:, T_BUS)); % vector of "to" buses
Yb = 1./(branch(:, BR_R) + 1j * branch(:, BR_X)); % branch series admittance
Yc = 1j*branch(:, BR_B); % branch shunt admittance
% Reactive power in shunt element at the "from" end
Yf = 1./tap.*(1./tap - 1).*Yb + Yc/2;
results.QSF = -imag(conj(Yf).*abs(V(f)).^2)*mpc.baseMVA;
% Reactive power in shunt element at the "to" end
Yt = (1 - 1./tap).*Yb + Yc/2;
results.QST = -imag(conj(Yt).*abs(V(t)).^2*mpc.baseMVA);

if ~exist('results', 'dir')
    mkdir('results');
end
save -v7 'results/results.mat' -struct results
% if exist('OCTAVE_VERSION', 'builtin')
%     % I'm in Octave
% else
%     % I'm in Matlab
%     save -v7 'results.mat' results
% end
