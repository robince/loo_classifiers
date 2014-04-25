function I = fastinfo(R)
% special case fast info to avoid overhead of building pars
% structure

[Nvar, Ntrl, Nstm] = size(R);

pars.Nt = Ntrl*ones(1,Nstm);
pars.biasCorrNum = 2; % PT
pars.testMode = false;

pars.doHR = true;
pars.doHRS = true;
pars.doHlR = false;
pars.doHlRS = false;
pars.doHiR = false;
pars.doHiRS = false;
pars.doChiR = false;
pars.doHshR = false;
pars.doHshRS = false;


[HR, HRS] = direct_method_v6a(R, pars);

% assume equally likely stimuli
I = HR - mean(HRS);