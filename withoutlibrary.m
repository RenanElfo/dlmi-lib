clc; clear;

% System Declaration
A = [0, 1; -2, -3]; size_A = size(A, 1);
B = [1; 1]; rows_B = size(B, 1); cols_B = size(B, 2);
C = [1, 1]; rows_C = size(C, 1); cols_C = size(C, 2);
D = 0; size_D = size(D, 1);
F = [A, B; zeros(cols_B, size_A+cols_B)]; size_F = size(F, 1);
G = [C, D];
T = 0.05;
Freq = 1/T;
x_0 = [1; 1];

% LMI Variables and Auxiliary Constants
setlmis([]);
n = 1;
delta_inv = Freq*2^n;
I_0 = eye(size_A, size_A+cols_B);
x_bar_0 = [x_0; zeros(cols_B, 1)];
[W, ~, sW] = lmivar(1, [size_A, 1]);
[Y, ~, sY] = lmivar(2, [size_A, cols_B]);
WY = lmivar(3, [sW, sY]);
gamma = lmivar(1, [1, 0]);
size_dlmi = 2^n+1;

%%%%%%%%
% DLMI %
%%%%%%%%

% DLMI Variable Declaration
Q = cell(1, 2^n+1);
sQ = cell(1, 2^n+1);
for k = 1:2^n+1
    [Q{k}, n_var, sQ{k}] = lmivar(1, [size(F, 1), 1]);
end

% DLMI Contraints
for k = 1:2^n
    tag = newlmi;
    lmiterm([tag, 1, 1, Q{k}], delta_inv, 1); % - P(t_{k})/delta
    lmiterm([tag, 1, 1, Q{k+1}], delta_inv, -1); % P(t_{k+1})/delta
    lmiterm([tag, 1, 1, Q{k}], F, 1, 's'); % F^T P(t_{k}) + P(t_{k}) F
    lmiterm([tag, 1, 2, Q{k}], 1, G');
    lmiterm([tag, 2, 2, 0], -1);

    tag = newlmi;
    lmiterm([tag, 1, 1, Q{k}], delta_inv, 1); % - P(t_{k})/delta
    lmiterm([tag, 1, 1, Q{k+1}], delta_inv, -1); % P(t_{k+1})/delta
    lmiterm([tag, 1, 1, Q{k+1}], F, 1, 's'); % F^T P(t_{k+1}) + P(t_{k+1}) F
    lmiterm([tag, 1, 2, Q{k+1}], 1, G');
    lmiterm([tag, 2, 2, 0], -1);
end
%%%%%%%

% LMI Constraints
tag = newlmi;
lmiterm([-tag, 1, 1, Q{2^n+1}], 1, 1); % Q(t_{2^n+1}) > 0
lmiterm([-tag, 1, 2, Q{2^n+1}], 1, I_0');
lmiterm([-tag, 2, 2, W], 1, 1);

tag = newlmi;
lmiterm([-tag, 1, 1, W], 1, 1);
lmiterm([-tag, 1, 2, WY], 1, 1);
lmiterm([-tag, 2, 2, Q{1}], 1, 1);

tag = newlmi;
lmiterm([-tag, 1, 1, gamma], 1, 1);
lmiterm([-tag, 1, 2, 0], x_bar_0');
lmiterm([-tag, 2, 2, Q{1}], 1, 1);

% Norm Minimization Problem Result
LMIs = getlmis;
c = zeros(n_var,1);
for i = 1:n_var
    c(i) = defcx(LMIs, i, gamma);
end
[copt, xopt] = mincx(LMIs, c);

copt;
Ysol = dec2mat(LMIs, xopt, Y);
Ysol2 = xopt(sY);
Wsol = dec2mat(LMIs, xopt, W);
K=Ysol'*inv(Wsol')
