phi_max = 180;
phi_min = 10;
sigma   = 10;      % motor noise SD (deg)
alpha   = 0.3;    % reward concavity
Dmax    = phi_max / phi_min;

% difficulty levels
d_names = {'easy','med','hard'};
q       = [0.9, 0.5, 0.3];    % baseline accuracy at 180°
gamma   = [1.0, 1.5, 2.5];    % tilt penalty exponent

% reward stakes
% R_low   = [0.15, 0.08, 0.05];    % per difficulty
R_low   = [0.05, 0.05, 0.05];    % per difficulty
R_high  = [0.5, 0.5, 0.5];

tilts = linspace(0,1,21);        % 21 tilt levels
phi   = phi_max - (phi_max-phi_min) * tilts;
D     = phi_max ./ phi;
Dtilde= (D - 1) / (Dmax - 1);

% motor term at 180°
p_motor_180 = normcdf(phi_max/2/sigma) - normcdf(-phi_max/2/sigma);

figure; hold on
for id = 1:numel(d_names)
    this_q     = q(id);
    this_gamma = gamma(id);
    this_Rlow  = R_low(id);
    this_Rhigh = R_high(id);

    % reward vs tilt
    g       = Dtilde .^ alpha;
    R_tilt  = this_Rlow + (this_Rhigh - this_Rlow) * g;

    % hit prob vs tilt
    p_motor = normcdf(phi/2/sigma) - normcdf(-phi/2/sigma);
    f       = p_motor / p_motor_180;
    p_hit   = this_q * (f .^ this_gamma);

    % expected reward
    EV      = p_hit .* R_tilt;

    % you can plot tilt vs EV here
    plot(tilts, EV, 'LineWidth', 2);
    xlabel('tilt (0 = wide wedge, 1 = narrow)');
    ylabel('expected reward (ml)');
end
legend(d_names);
