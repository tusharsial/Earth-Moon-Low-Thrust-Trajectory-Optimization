function J = Hybrid_Cost_Function(s)

% Extract parameters;
t_esc = s(9);
t_capt = s(10);

% Cost function
J = t_esc + t_capt;
