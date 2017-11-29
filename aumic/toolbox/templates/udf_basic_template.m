function [rkin,flowin,flowout,xfeed,kla,xsat,dxdt_more] = ...
    udf_basic(t,kin,x,model)

% initialization
rkin = zeros(length(kin.kmax.num),1);
kmax = kin.kmax.num; K = kin.K.num;
nx = length(x);
flowin = 0; flowout = 0; xfeed = zeros(nx,1); vol = 1;
kla = zeros(nx,1); xsat = zeros(nx,1);
dxdt_more = zeros(nx,1);

%--------------------------------------------------------------------------
% rkin

% flowin/out and xfeed

% kla and xsat

% dxdt_more

%==========================================================================

% example:
% kmaxG1 = kmax(1); kmaxG2 = kmax(2); kmaxG3 = kmax(3);
% kmaxX1 = kmax(4); kmaxX2 = kmax(5); kmaxX3 = kmax(6); kmaxX4 = kmax(7);
% KG1 = K(1); KG2 = K(2); KG3 = K(3);
% KX1 = K(4); KX2 = K(5); KX3 = K(6); KX4 = K(7);
% KinhG = K(8); KinhX = K(9);
%
% GLC = x(4); XYL = x(8); ETH = x(9);
%
% rkin(1) = kmaxG1*GLC/(KG1+GLC)/(1+ETH/KinhG);
% rkin(2) = kmaxG2*GLC/(KG2+GLC)/(1+ETH/KinhG);
% rkin(3) = kmaxG3*GLC/(KG3+GLC)/(1+ETH/KinhG);
% rkin(4) = kmaxX1*XYL/(KX1+XYL)/(1+ETH/KinhX);
% rkin(5) = kmaxX2*XYL/(KX2+XYL)/(1+ETH/KinhX);
% rkin(6) = kmaxX3*XYL/(KX3+XYL)/(1+ETH/KinhX);
% rkin(7) = kmaxX4*XYL/(KX4+XYL)/(1+ETH/KinhX);
% 
% dilution term (step change)
% if t<20
%     dilrate = 0.5;
% else
%     dilrate = 0.8;
% end
% xfeed(4) = 200;
% 
% % dilution and feeding
% period = 5; % 5 hours
% dilrate = 0.5+0.3*sin(2*pi*t/period); % period of sin(a*t) = 2*pi/|a|
% xfeed(4) = 200;
%
% dilrate = 0.5;
% period = 10; % 10 hours
% xfeed(4) = 200+100*sin(2*pi*t/period); % period of sin(a*t) = 2*pi/|a|
% 
% mass transfer of gaseous components (CO2) 
% xsat(3) = 0.0075*1000/44; % mmol/liter
% kla(3) = 350; % 1/hour

% dxdt_more 
