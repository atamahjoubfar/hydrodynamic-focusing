clearvars
%{
Please refer to Lee's 2006 paper titled "The hydrodynamic focusing effect inside rectangular microchannels" for definition of the parameters and equations. 
This Matlab code is written by Ata Mahjoubfar (ata.m@ucla.edu).	
%}

% Channel Dimensions:
w0 = 200e-6;
% ^^^^^^^^^^^^^^^
h = 25e-6;
% ^^^^^^^^^^^^^^^
epsilon = h/w0;

% Flow rates:
% Set the pump rate based on the sheath syringe size:
pump_set_rate = 318e-9/60; % in m^3/second
% ^^^^^^^^^^^^^^^
Qs1 = pump_set_rate/2; % Flow rate for right sheath channel
Qs2 = pump_set_rate/2; % Flow rate for left sheath channel
sheath_syringe = '30_ml_BD_Plastic';
% ^^^^^^^^^^^^^^^
switch sheath_syringe
    case '1_ml_BD_Plastic'
        sheath_syringe_diameter = 4.78e-3;
    case '3_ml_BD_Plastic'
        sheath_syringe_diameter = 8.66e-3;
    case '5_ml_BD_Plastic'
        sheath_syringe_diameter = 12.06e-3;
    case '10_ml_BD_Plastic'
        sheath_syringe_diameter = 14.5e-3;
    case '20_ml_BD_Plastic'
        sheath_syringe_diameter = 19.13e-3;
    case '30_ml_BD_Plastic'
        sheath_syringe_diameter = 21.7e-3;
end
sample_syringe = '3_ml_BD_Plastic';
% ^^^^^^^^^^^^^^^
switch sample_syringe
    case '1_ml_BD_Plastic'
        sample_syringe_diameter = 4.78e-3;
    case '3_ml_BD_Plastic'
        sample_syringe_diameter = 8.66e-3;
    case '5_ml_BD_Plastic'
        sample_syringe_diameter = 12.06e-3;
    case '10_ml_BD_Plastic'
        sample_syringe_diameter = 14.5e-3;
    case '20_ml_BD_Plastic'
        sample_syringe_diameter = 19.13e-3;
    case '30_ml_BD_Plastic'
        sample_syringe_diameter = 21.7e-3;
end
Qi = pump_set_rate.*(sample_syringe_diameter/sheath_syringe_diameter).^2; % Flow rate for the sample channel

n_max = 10;
% ^^^^^^^^^^^^^^^
n_matrix = 0:n_max;
wf_start = w0/4;
% ^^^^^^^^^^^^^^^
gamma = @(wf) (1-(192*h/pi^5/wf)*sum(sinh((2*n_matrix+1)*pi*wf/2/h)./...
    cosh((2*n_matrix+1)*pi*w0/2/h)./(2*n_matrix+1).^5))/...
    (1-(192*h/pi^5/w0)*sum(tanh((2*n_matrix+1)*pi*w0/2/h)./(2*n_matrix+1).^5));
solve_function = @(wf) gamma(wf)*wf/w0 - Qi/(Qi+Qs1+Qs2);
options = optimset('Display','off');
wf = fsolve(solve_function, wf_start, options);
gamma = gamma(wf);
vf_bar = Qi/h/wf; % Average sample fluid velocity after focusing
v0_bar = (Qi+Qs1+Qs2)/w0/h; % Average total fluid velocity after focusing. This is equal to u_bar_bar!
minus_dpdx_over_mu = v0_bar./(h^2/12*(1-192*h/pi^5/w0.*sum(tanh((2.*n_matrix+1).*pi*w0/2/h)./((2*n_matrix+1).^5))));
resolution = 101;
% ^^^^^^^^^^^^^^^
if mod(resolution,2)==0
    resolution = resolution+1;
end
[y,z] = meshgrid(linspace(-w0/2,w0/2,resolution),linspace(-h/2,h/2,resolution));
u = zeros(size(y));
for n = n_matrix
    u = u + 4*h^2/pi^3*minus_dpdx_over_mu.*sum((-1).^n.*(1-cosh((2*n+1)*pi*y/h)./cosh((2*n+1)*pi*w0/2/h)).*cos((2*n+1)*pi*z/h)./((2*n+1).^3),3);
end
figure('Name','Velcity Profile', 'Renderer', 'zbuffer')
pcolor(1e6.*y,1e6.*z,u)
axis image
colormap jet
colorbar
shading interp
xlabel('y (\mum)')
ylabel('z (\mum)')
rectangle('Position',[-1e6.*wf/4,-1e6.*h/2,1e6.*wf/2,1e6.*h],'EdgeColor',[1 1 1],'LineStyle',':')

disp(['Maximum velocity is ',num2str(u((resolution+1)/2,(resolution+1)/2)),' m/s.'])
disp(['Sample flow width is ',num2str(1e6.*wf),' um.'])
