%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Naren Ganesh, Modified Jul. 9, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define lambda, k, and n matricies of the same size

wn = 10000./lambda;
d = 10e-15; % thickness of the film (um)
angle = 8; %FTIR incidence angle
theta = angle * pi/180; %convert to radians 
for i = 1:length(n)
    n_i = n(i);
    k_i = k(i);
    lambda_j = lambda(i);
%Calculate r12 for TE polarization
r12s = (cos(theta)-n(i)*cos(asin(sin(theta)/n(i))))/(cos(theta)+n(i)*cos(asin(sin(theta)/n(i))));
%Calculate r12 for TM polarization
r12p = (n(i)*cos(theta)-cos(asin(sin(theta)/n(i))))/(n(i)*cos(theta)+cos(asin(sin(theta)/n(i))));
%Calculate reflectivity averging over both lin. polarizations
rho = ((r12s)^2 + (r12p)^2)/2;
tau = exp(-((4*pi*k(i)*d)/(lambda_j)*cos(asin(sin(theta)/n(i)))));
%Reverse-propogated Reflectance
R_rp = rho*(1+((1-rho)^2*tau^2)/(1-rho^2*tau^2));
%Reverse-propogated Transmittance
T_rp = ((1-rho)^2*tau)/(1-rho^2*tau^2);
%Store Values
R_rp_values(i) = R_rp;
T_rp_values(i) = T_rp;
end

%figure;
%plot(lambda,R_rp_values)
%xlim([0,25])

% Initialize arrays to store results
n_rp_values = zeros(size(R_rp_values));
k_rp_values = zeros(size(R_rp_values));

%Step 2

for j = 1:length(R_rp_values)
    R_rp_values_j = R_rp_values(j);
    T_rp_values_j = T_rp_values(j);
    lambda_j = lambda(j);
    
    % Calculate R12
   R12_rp = (2+T_rp_values_j^2 - (1-R_rp_values_j)^2 - sqrt(((2+T_rp_values_j^2 - (1-R_rp_values_j)^2))^2 - 4*R_rp_values_j*(2-R_rp_values_j)))/(2*(2-R_rp_values_j));

%Step 3
   
    % Calculate extinction coefficient k
    k_rp = ((lambda_j / (4 * pi *d))*log((R12_rp*T_rp_values_j)/(R_rp_values_j - R12_rp)));

    % Calculate real part of refractive index n
    n_rp = (((1+R12_rp)/(1-R12_rp)) + sqrt(((1+R12_rp)/(1-R12_rp))^2 -1-k_rp^2));

    % Store the results
    n_rp_values(j) = n_rp;
    k_rp_values(j) = k_rp;
    R12_rp_values(j) = R12_rp;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constants and Coefficants
c = 3e14;%speed of light in um
w = (2*pi*c)./lambda; %Angular Frequency
w_upper = w(590); %top of measurement range
w_lower = w(1); %bottom of measurement range
I = -w./pi; %Integration Coefficant

%Refractive Index
N = real(n_rp_values)+1i*real(k_rp_values); %Define refractive index
N_r = (N-1)./(N+1); %Fresnel coefficant relation


%Step 4
beta_1 = 0; %Initialization 
beta_2 = 0;
%Two chosen points
point_1 = 586;
point_2 = 569;
omega_1 = w(point_1); %Storing angular frequency at point 1
omega_2 = w(point_2); %Storing angular frequency at point 2
%Intial Delta @ point 1 & 2
delta_1 = -1i*log((N(point_1)-1)/((N(point_1)+1)*sqrt(R12_rp_values(point_1))));
delta_2 = -1i*log((N(point_2)-1)/((N(point_2)+1)*sqrt(R12_rp_values(point_2))));

%Step 5

beta = zeros(size(w));
for z = 1:length(w)
for a = 1:length(w)
f = @(x) I(z)*(log(R12_rp_values(a))-log(R12_rp_values(z)./(w(a)^2-w(z)^2))) ;
end
    for a = 1:length(w)
        x =  a;
        %Trapaziodal sum
        beta(z) = beta(z) + ((f(x)+f(x+1))/2)*10e-15; 

    end
end

%Step 6
    syms c1 c2 
    %System of equation solver
    %Use that a + b + c == delta
eqn1 = (c1+(log(R12_rp_values(point_1))/(2*pi)))*log((omega_1-w_lower)/(omega_1+w_lower)) + (c2-(log(R12_rp_values(point_1))/(2*pi)))*log((w_upper-omega_1)/(w_upper+omega_1))+beta(point_1) == 0;
eqn2 = (c1+(log(R12_rp_values(point_2))/(2*pi)))*log((omega_2-w_lower)/(omega_2+w_lower)) + (c2-(log(R12_rp_values(point_2))/(2*pi)))*log((w_upper-omega_2)/(w_upper+omega_2))+beta(point_2) == 0;
%Two unkowns (c1 and c2) thus two equations
sol = solve([eqn1, eqn2], [c1, c2]); %Solve  
c1Sol = sol.c1; %Store in Symbolic
c2Sol = sol.c2;
ca = double(c1Sol); %Converts from Symbolic to Matlab Double
cg = double(c2Sol); 

%Store alpha and gamma values
alpha_values = zeros(size(w)); %Initialize
gamma_values = zeros(size(w));
for b = 1:length(w)
    R12_rp_values_b = R12_rp_values(b);
    w_b = w(b);
alpha = (ca+(log(R12_rp_values_b)/(2*pi)))*log((w_b-w_lower)/(w_b+w_lower));
gamma = (cg-(log(R12_rp_values_b)/(2*pi)))*log((w_upper-w_b)/(w_upper+w_b));
alpha_values(b) = alpha; 
gamma_values(b) = gamma; 

end
del_norm = alpha_values+gamma_values+beta; %calculated phase angle of fresnel coefficant through KK
KK_n = (1-R12_rp_values)./(1+R12_rp_values-2.*sqrt(R12_rp_values).*cos(del_norm));
KK_k = (2.*sqrt(R12_rp_values).*sin(del_norm))./(1+R12_rp_values-2.*sqrt(R12_rp_values).*cos(del_norm));

n_rp_first = ones(1,590);
for e = 1:526
    n_rp_first(e) = n_rp_values(e);
end
for e = 526:590
    n_rp_first(e) = inf;
end

kk_first = ones(1,590);
for e = 1:526
    kk_first(e) = inf;
end
for e = 526:569
    kk_first(e) = KK_n(e);
end
for e = 572:590
    kk_first(e) = inf;
end

n_rp_second = ones(1,590);
for e = 1:572
    n_rp_second(e) = inf;
end
for e = 570:590
    n_rp_second(e) = n_rp_values(e);
end

k_rp_first = ones(1,590);
for e = 1:522
    k_rp_first(e) = k_rp_values(e);
end
for e = 523:530
    k_rp_first(e) = inf;
end
for e = 531:590
    k_rp_first(e) = k_rp_values(e);
end

kk_k_first = ones(1,590);
for e = 1:520
    kk_k_first(e) = inf;
end
for e = 521:532
    kk_k_first(e) = KK_k(e);
end
for e = 533:590
    kk_k_first(e) = inf;
end



figure;
h = plot(flip(lambda),n_rp_first,'r',flip(lambda),kk_first,'r',flip(lambda),n_rp_second,'r',flip(lambda),n,'k');
xlabel('Wavelength, \lambda(\mum)');
ylabel('n');
legend( h, 'Nichelatti''s Method + Kramers-Kronig Analysis', 'Querry et al.', 'Location', 'NorthWest');

figure;
h2 = plot(flip(lambda),k_rp_first,'r',flip(lambda),kk_k_first,'r',flip(lambda),k,'k');
ylabel('k');
legend( h2, 'Nichelatti''s Method + Kramers-Kronig Analysis', 'Querry et al.', 'Location', 'NorthWest');


%%%%%%%%%%%%%%%%%%%%%%%%
