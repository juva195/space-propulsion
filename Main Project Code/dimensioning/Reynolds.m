function Re_t = Reynolds(D_t)


set(groot, 'defaultTextInterPreter', 'tex')
set(groot, 'defaultAxesTickLabelInterPreter', 'tex')
set(groot, 'defaultFigureColormap', turbo(256));
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', [1; 1; 1]);
set(groot, 'defaultAxesColor', 'none');
set(groot,'DefaultAxesFontSize', 15);

%% PATH GENERATION

if isunix()
    CEApath = "CEA/";
else
    CEApath = "CEA\";
end

addpath(genpath(CEApath));


ceaOutout = CEA('problem','rocket','frozen','nfz',1,'p,pa',50,'o/f', ...
2.24,'sup,ae/at',65,...
    'reactants','fuel','RP-1(L)', 'C',1.,'H',1.9423,...
    'oxid','O2(L)','output', 'tran','end');
rho_t = ceaOutout.output.froz.density(2);
v_t = ceaOutout.output.froz.sonvel(2);
mu_t = ceaOutout.output.froz.viscosity(2) * 10^-3;

Re_t = rho_t * v_t * D_t/mu_t; 
end
