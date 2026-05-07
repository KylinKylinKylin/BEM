% This code determines the projectors and checks it with the existing code

% Initialize sphere
diameter = 10;
p = trisphere(144, diameter);

%plot( p, 'EdgeColor', 'b' ); hold on
%plot( p, 'EdgeColor', 'b', 'nvec', 1 ); % show outer surface normals

lambda = 900;
k0 = 2*pi/lambda;

epsz = @(w) lorenze_drude(299792458*2*pi*6.58e-16/w, 'Au', 'LD');
epsz_GaAs = @(w) GaAs_Adachi(299792458*2*pi*6.58e-16/w);

mat1 = Material(1, 1);
mat2 = Material(epsz_GaAs, 1);
mat = [mat1, mat2];

tau = BoundaryEdge(mat, p, [2,1]);

bem = galerkin.bemsolver( tau, 'relcutoff', 2, 'waitbar', 1 );

%%
% planewave excitation 
dir = [1,0,0];
pol = [0,0,1];
% if statement to stop if pol and dir are not 90deg diff

exc = galerkin.planewave(pol, dir);

% solve BEm equation
lambda = 900; %setting wavelength of light in nm
k0 = 2 * pi / lambda;
sol = bem \ exc( tau, k0 );


%% Don't need in the section
% points
n = 51;
xx = 1 * diameter * linspace( -1, 1, n);
[ x, y ] = ndgrid( xx, xx );
pt = Point( tau, [ x( : ), y( : ), 0 * x( : ) ] );

% size(pt) 2601 x 1

% evaluate fields
[e, h] = fields( sol, pt, 'relcutoff', 2, 'waitbar', 1); 
[ei, hi] = fields(exc, pt, k0); % incident field


ee_int = real( squeeze( dot( e, e, 2 ) ) );
ee = real( squeeze( dot( ei, ei, 2 ) ) );
%%

% Compute S and D operators
inout = 1;  % 1 = exterior, 2 = interior

SD_ex = eval(bem.pot, k0, inout);

S1_ex = SD_ex.SL1; % times mu
S2_ex = SD_ex.SL2; % times epsilon
CC_ex  = SD_ex.DL;

inout = 2;  % 1 = exterior, 2 = interior

SD_in = eval(bem.pot, k0, inout);

S1_in = SD_in.SL1;
S2_in = SD_in.SL2;
CC_in  = SD_in.DL;

disp("Part1 done...")
%% Making the projctors only for the electric fields
% Projecting the external field inside

I1 = eye(size(CC_in, 1));
s0 = zeros(size(CC_in));

% full calderon matrix
% tilde D [ f_E; f_H ]
tildeD_ps = [0.5 * I1 - CC_ex, 1j * k0 * S1_ex; 0.5 * I1 - CC_ex, -1j * k0 * S1_ex];
tildeD_mi = [0.5 * I1 + CC_in, -1j * k0 * S1_in; 0.5 * I1 + CC_in, 1j * k0 * S1_in];

% (1 - D_s_max) * D_i_min * g_i = pm (1 - D_s_max) * f_0
% The equation to solve.

f = exc(tau, k0);

I2 = eye(size(tildeD_ps, 1));

f_e = f.e;
f_h = f.h;

f_0 = [f_e; f_h];

% testing to see if D is orthogonal
tol = 1e-10;
cond = conj(tildeD_ps) * tildeD_mi;
isOrthogonal = all(abs(cond(:)) < tol);


if isOrthogonal == true
    disp('True');
end

% (1 − Ds+ )i = ±(1 − Ds+ )f 0 ,
rhs1 = (I2 - tildeD_ps) * f_0;
lhs1 = (I2 - tildeD_ps) * tildeD_mi;

g_i = lhs1 \ rhs1;

% -(1 − Di- )s = ±(1 − Di- )f 0 ,
rhs2 = (I2 - tildeD_mi) * f_0;
lhs2 = -(I2 - tildeD_mi) * tildeD_ps;

g_s = lhs2 \ rhs2;

writematrix(f_0, '~/BEM_1/Data/f_0.csv');

disp('Part2 done...')

%%

%Make calderon matrix
% write it out first!!!
% I = sum |l><l|
% for matrix mul between D_s, D_i
% otherwise I = sum |l'><l|

% sum f_v' (x') f_v(x) verify where these can be from the exc edge vertices