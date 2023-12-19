 v0 = 4.5;
 vA0 = 0.3.*v0;
 vV0 = 0.7.*v0;
vVU0 = 1;
Pa0 = 90;
Pv0 = 7;
xQ0 = 0;
Q0 = 5;
xR0 = 0;
R0 = 20;
xVU0 = 0;
Iphp0 = 0;
tspan = [0 60];


time_vals = 0:1:tspan(end);
SNP_infratevals = zeros(size(tspan));
for i = 1:length(time_vals)
    if time_vals(i) > 50/60 && time_vals(i) <= 10
		SNP_infratevals(i) = 30;
    elseif time_vals(i) > 10 && time_vals(i) <= 20 
		SNP_infratevals(i) =20;
	elseif time_vals(i) > 20 && time_vals(i) <= 25 
		SNP_infratevals(i) =6*(time_vals(i)-20)+20; % function
	elseif time_vals(i) > 20 && time_vals(i) <= 35
		SNP_infratevals(i) =50;
	elseif time_vals(i) > 35 && time_vals(i) <= 45 
		SNP_infratevals(i) = -4*(time_vals(i)-35)+50; % function
	elseif time_vals(i) > 45 && time_vals(i) <= 50
		SNP_infratevals(i) =10;
	elseif time_vals(i) > 50 && time_vals(i) <= tspan(end)
	    SNP_infratevals(i) = 15;
    end
end

PHP_infratevals = zeros(size(tspan));
    for i = 1:length(time_vals)
        if time_vals(i) > 0 && time_vals(i) <= 25 
		    PHP_infratevals(i) =0;
        elseif time_vals(i) > 25 && time_vals(i) <= 40
		    PHP_infratevals(i) = 2;
	    elseif time_vals(i) > 40 && time_vals(i) <= tspan(end)
		    PHP_infratevals(i) = 0;
        end
    end


ft = linspace(0,70,61);
gt = linspace(0,70,61);

y0 = [0; 0; 0; vA0; vV0; xQ0; xR0; xVU0; Iphp0];
[t,y] = ode45(@(t,y) CVODE(t,y,ft,gt, SNP_infratevals, PHP_infratevals), tspan, y0);
x1 = y(:,1);
x2 = y(:,2);
x3 = y(:,3);
vA = y(:,4);
vV = y(:,5);
xQ = y(:,6);
I = y(:,9);

[~,MAP_transposed] = CVODE(t.',y.',ft,gt, SNP_infratevals, PHP_infratevals);
MAP = MAP_transposed';


figure
colororder({'k','k'})
plot(t, MAP, 'r');
xlabel('Time [min]');
ylabel('MAP [mmHg]');
axis([0 tspan(end) 0 100])
legend('MAP')
%hold on
%yyaxis right
figure
plot(time_vals,SNP_infratevals,'k');
ylabel('PHP infusion rate [ml/h]')
axis([0 tspan(end) 0 100])
legend('PHP');
%hold off

function [dydt, Pa] = CVODE(t, y, ft,gt, SNP_infusion, PHP_infusion)
    x1 = y(1,:);
    x2 = y(2,:);
    x3 = y(3,:);
    vA = y(4,:);
    vV = y(5,:);
    xQ = y(6,:);
    xR = y(7,:);
    xVU = y(8,:);
    Iphp = y(9,:);

% ----- initaial values ------ %
    v0 = 4.5;
    vA0 = 0.3.*v0;
    vV0 = 0.7.*v0;
    R0 = 20;
    Pa0 = 90;
    Pv0 = 7;
    Q0 = 5;
    

% ----- parameters from fig5 ----- %
    Ka = 7;
    Kv = 8;
    Kq = 0.02;
    pQ = 0.5;
    zQ = 3;
    Kr = 1.5;
    tauR = 150;
    Kvu = 1;

    % PHP parameters
    kphp = 0.2; %[min^-1]
    I_sigmaR = 0.4; % [mcg/kg/min]
    lambdaR = 2; %constant
    sigma = 50; % percent of max vasoconstriction and venoconstriction effects
    k_sigma = -0.01.*sigma./log(3);
    
    %SNP parameters
    tau1 = 50; %time constant of SNP action [s]
    tau2 = 10; % time constant for flow through pulmonary circulation [s]
    tau3 = 30; % time constant for flow through systemic circulation [s]
    alpha = 0.5; % fraction of SNP recirculated
    G = 3; %plant gain [mmHg/(ml/h)]
  
    
    % -- Unstressed blood volume - 
    
    tau_vu = 7;
    I_sigmaVU = 0.4;
    lambdaVU = 1.3;
    eta_vu = 0.21;
    Dvu = eta_vu.*v0.*(Iphp.^lambdaVU./(I_sigmaVU.^lambdaVU + Iphp.^lambdaVU));
    

% ----- Interpolation ----- %
    SNP_infusion = interp1(ft, SNP_infusion, t);
    PHP_infusion = interp1(gt, PHP_infusion, t);
    
    % ----- PHP admin ----- %
    Jphp = PHP_infusion; % infusion rate
    xPHP = 1 - Iphp.^lambdaR./(I_sigmaR.^lambdaR + Iphp.^lambdaR);
    Dr = R0.*k_sigma.*log(xPHP./(2-xPHP));
    assignin('base','Dr', Dr);
% ----- Physiological model equations ----- %
    delta_vA = vA-vA0;
    delta_vV = vV - vV0;
    delta_vVU = xVU + Dvu;
    
    delta_Pv = Kv.*(delta_vV- delta_vVU);
    delta_Q = Kq.*(xQ+delta_Pv);
    delta_R = xR + Dr;
   
    delta_PaSNP = -G*x1;
    delta_PaPHP = Ka.*delta_vA;
    delta_Pa = delta_PaPHP +delta_PaSNP;
    
    Pa = delta_Pa + Pa0;
    Pv = Pv0 + Kv.*(vV-vV0);
    Q = delta_Q + Q0;   
    R = (delta_R + R0);  
   
    % ------------ diff equations ----------- %
    dx1dt = 1./tau1.* (x2-x1);
    dx2dt = 1./tau2.* (SNP_infusion-x2-(alpha.*x3));
    dx3dt = 1./tau3 .* (x2-x3);

    dvAdt = Q -((Pa-Pv)./R);
    dvVdt = (Pa-Pv)./R - Q;
    dxQdt = -pQ.*xQ + (zQ - pQ).*delta_Pv;
    dxRdt = -tauR.*xR - Kr.*delta_PaPHP;
    dxVUdt = -tau_vu*xVU + Kvu*delta_PaPHP;
    dIphpdt = -kphp.*Iphp+kphp.*Jphp;
   
    dydt = [dx1dt; dx2dt; dx3dt; dvAdt; dvVdt; dxQdt; dxRdt; dxVUdt; dIphpdt];
end

