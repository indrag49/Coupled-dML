function out = DenaturedFast
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
dydt=[kmrgd(1)^2*(1-kmrgd(1))-kmrgd(2)+par_I;
par_A*exp(par_Alpha*kmrgd(1))-par_Gamma*kmrgd(2);];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(DenaturedFast);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
jac=[ - 2*kmrgd(1)*(kmrgd(1) - 1) - kmrgd(1)^2 , -1 ; par_A*par_Alpha*exp(kmrgd(1)*par_Alpha) , -par_Gamma ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
jacp=[ 0 , 0 , 0 , 1 ; exp(kmrgd(1)*par_Alpha) , kmrgd(1)*par_A*exp(kmrgd(1)*par_Alpha) , -kmrgd(2) , 0 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
hess1=[ 2 - 6*kmrgd(1) , 0 ; par_A*par_Alpha^2*exp(kmrgd(1)*par_Alpha) , 0 ];
hess2=[ 0 , 0 ; 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
hessp1=[ 0 , 0 ; par_Alpha*exp(kmrgd(1)*par_Alpha) , 0 ];
hessp2=[ 0 , 0 ; par_A*exp(kmrgd(1)*par_Alpha) + kmrgd(1)*par_A*par_Alpha*exp(kmrgd(1)*par_Alpha) , 0 ];
hessp3=[ 0 , 0 ; 0 , -1 ];
hessp4=[ 0 , 0 ; 0 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
hessp(:,:,4) =hessp4;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
tens31=[ -6 , 0 ; par_A*par_Alpha^3*exp(kmrgd(1)*par_Alpha) , 0 ];
tens32=[ 0 , 0 ; 0 , 0 ];
tens33=[ 0 , 0 ; 0 , 0 ];
tens34=[ 0 , 0 ; 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_A,par_Alpha,par_Gamma,par_I)
