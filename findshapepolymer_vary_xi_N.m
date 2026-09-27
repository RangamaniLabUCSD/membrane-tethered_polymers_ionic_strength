clc;
clear;
close all;

%% membrane properties
kappa = 20;            % bending rigidity; K_BT         
sigma = 0.012;         % Surface tension: K_BT/nm^2  
lambda = sqrt(kappa/sigma);    % Unit: nm  
f0 = sqrt(kappa*sigma);     % Unit: K_BT/nm=pN
gamma = 0.0;     % line tension
c_0 = 0.00;
Rb = 1000;         
J0 = Rb;
f = 0;
%% Polymer properties
a = 10;  
a_monomer = 0.35;  
N_ionizable = a/a_monomer;
Chi = 0;   
v = (1-2*Chi)*a^3;     

R_s = 100;        %%   membrane patch radius,  Unit: nm
A_s = pi*R_s^2;   %%   membrane patch area,  Unit: nm^2
R_v = R_s/2;

alpha_amino_span = 0.0;
alpha_amino = alpha_amino_span;
Con_ion_span = 150;
Con_ion = Con_ion_span;
Phi_ion = Con_ion*6.022*10^(-4); %% unit: #/nm^3
alpha_bulk = alpha_amino*N_ionizable; %%  the degree of ionization
W = v + alpha_bulk^2/(2*Phi_ion);

N_span = 20;
% N_span = linspace(15,50,36);   
k = 0;
for N = N_span 
    N

    % xi_span = linspace(10,30,100);
    %xi_span = 10:0.5:30;
    xi_span = 30;
    j=0;
    for xi = xi_span 
        j=j+1;
        xi  
        S_R = xi^2;    
        rho_g = 1/xi^2;  
        N_p = A_s/S_R;   

        alpha = 0.001;   
        Rp = R_s/sqrt(2*(1-cos(alpha)));
        alpha_jiaodu = alpha*180/pi;
        r0 = Rp*sin(alpha);
        %% Solve equations by using bvp5c
        yeq = @(u,y,para) shapepolymer(u,y,para,sigma,kappa,f);
        ybc = @(ya,yb,para) twobc(ya,yb,para,alpha,r0,Rb);
        yinit = @(u) guess(u,Rb,alpha,r0);
        solinit = bvpinit(linspace(0,1,100),yinit,J0);
        opts = bvpset('RelTol',1e-5,'AbsTol',1e-10,'NMax',10000);
        sol = bvp5c(yeq,ybc,solinit,opts);

        k=k+1;
        alpha_c = 0.000:0.01:3.021;
        i = 0;
        E_tot = 0;
        m_alpha = alpha_c';
        r_all = cell(length(alpha_c),1);
        psi_all = cell(length(alpha_c),1);
        z_all = cell(length(alpha_c),1);
        for alpha = alpha_c
            alpha
            i = i + 1;
            if alpha == 0
                Rp = R_s/sqrt(2*(1-cos(alpha)));  %% radius of spherical cap
                R_cap(i) = Rp;
                H_brush(i) = N*S_R^(-1/3)*(W*a^2/3)^(1/3);
                Z = 0;
                z_free(i) = 0;
                Height(i) = Z+Rp*(1-cos(alpha));     %% height of the system
                E_polymer_cap(i) = 9*N/(2*kappa)*(R_s/xi)^2*(a^2/3/xi^2)^(2/3);
                Eb_cap(i) = 4*(1-cos(alpha));     % bending energy of spherical cap part
                Eb_free(i) = 0;   % bending energy of free part
                Et_cap(i) = 1/2*sigma*R_s^2*(1-cos(alpha))/kappa;  % surface tension energy of spherical cap part
                Elt_cap(i)= gamma*R_s/(lambda*f0)*sqrt(2*(1+cos(alpha)));  % line tension energy of spherical cap part
                Et_free(i) = 0;  % surface tension energy of free part
                E_force_cap(i) = -f*A_s*sqrt((1-cos(alpha))/2)/kappa;  
                E_force_free(i) = 0;
                E_force_free_Z(i)  = -f*Z/kappa;
                E_tot_force(i) = E_force_cap(i)+E_force_free(i);  
                E_tot_cap(i) = Eb_cap(i)+Et_cap(i)+Elt_cap(i)+E_force_cap(i);  %total energy of spherical cap part
                E_mem_cap(i) = 4*(sqrt(1-cos(alpha))-sqrt(2)*R_s*c_0/4)^2 + R_s^2/(2*lambda^2)*(1-cos(alpha)) + gamma*R_s/(lambda*f0)*sqrt(2*(1+cos(alpha)));
                E_cap(i) = 9*N/(2*kappa)*(R_s/xi)^2*(a^2/3/xi^2)^(2/3) + 4*(sqrt(1-cos(alpha))-sqrt(2)*R_s*c_0/4)^2 + R_s^2/(2*lambda^2)*(1-cos(alpha))...
                    + gamma*R_s/(lambda*f0)*sqrt(2*(1+cos(alpha)));
                E_tot1(i) = E_cap(i);  %total free energy
            else
                Rp = R_s/sqrt(2*(1-cos(alpha)));  %% radius of spherical cap
                R_cap(i) = Rp;
                H_flat = N*S_R^(-1/3)*(W*a^2/3)^(1/3);
                H_brush(i) = Rp*(1+5/3*H_flat/Rp)^(3/5)-Rp;
                %% calculate the energy of the spherical cap part
                E_polymer_cap(i) = 9/(2*kappa)*(R_s/xi)^2*(R_s/sqrt(2*(1-cos(alpha))))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+5*N/(3*R_s/sqrt(2*(1-cos(alpha))))*(v*a^2/3/xi^2)^(1/3))^(1/5)-1);
                E_mem_cap(i) = 4*(sqrt(1-cos(alpha))-sqrt(2)*R_s*c_0/4)^2 + R_s^2/(2*lambda^2)*(1-cos(alpha)) + gamma*R_s/(lambda*f0)*sqrt(2*(1+cos(alpha)));
                E_cap(i) = 9/(2*kappa)*(R_s/xi)^2*(R_s/sqrt(2*(1-cos(alpha))))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+5*N/(3*R_s/sqrt(2*(1-cos(alpha))))*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                    + 4*(sqrt(1-cos(alpha))-sqrt(2)*R_s*c_0/4)^2 + R_s^2/(2*lambda^2)*(1-cos(alpha)) + gamma*R_s/(lambda*f0)*sqrt(2*(1+cos(alpha)));
                %% calculate the energy of the free part
                alpha_xaixs(i) = alpha;
                r0 = Rp*sin(alpha);
                para = sol.parameters;
                yeq = @(u,y,para) shapepolymer(u,y,para,sigma,kappa,f);
                ybc = @(ya,yb,para) twobc(ya,yb,para,alpha,r0,Rb);
                sol = bvp5c(yeq,ybc,sol,opts);
                u = sol.x;
                psi = sol.y(1,:);
                dpsi = sol.y(2,:);
                J = sol.y(3,:);
                r = sol.y(4,:);
                z = sol.y(5,:);
                err_alpha(i) = sol.stats.maxerr;
                r_all{i} = sol.y(4,:);
                z_all{i} = sol.y(5,:);
                psi_all{i} = sol.y(1,:);
                psi_0 = psi(:,1);
                dpsi_0 = dpsi(:,1);
                R_0 = r(:,1);
                Z = z(:,1);
                z_free(i) = Z;
                Height(i) = Z+Rp*(1-cos(alpha));     %% height of the system
                Eb_cap(i) = 4*(1-cos(alpha));     % bending energy of spherical cap part
                Eb_free(i) = trapz(u,(dpsi./J+sin(psi)./r).^2.*r.*J);   % bending energy of free part
                Et_cap(i) = sigma*Rp^2*(1-cos(alpha))^2/kappa;  % surface tension energy of spherical cap part
                Elt_cap(i)= gamma*R_s/(lambda*f0)*sqrt(2*(1+cos(alpha)));  % line tension energy of spherical cap part
                Et_free(i) = trapz(u,2*sigma*r.*J.*(1-cos(psi))/kappa);  % surface tension energy of free part
                E_force_cap(i) = -f*Rp*(1-cos(alpha))/kappa;   
                E_force_free(i) = -trapz(u,f*J.*sin(psi)/kappa); 
                E_force_free_Z(i)  = -f*Z/kappa;
                E_tot_force(i) = E_force_cap(i)+E_force_free(i); 

                E_tot_cap(i) = Eb_cap(i)+Et_cap(i)+Elt_cap(i)+E_force_cap(i);  %total energy of spherical cap part
                E_tot_free(i) = Eb_free(i)+Et_free(i)+E_force_free(i);   %total energy of free part
                E_tot1(i) = E_cap(i) + Eb_free(i) + Et_free(i) + E_force_cap(i) + E_force_free(i);  %total free energy
                E_tot2(i) = E_polymer_cap(i) + Eb_cap(i) + Et_cap(i) + Elt_cap(i) + Eb_free(i) + Et_free(i);  
            end
        end

        E_tot1 = E_tot1';
        [E_tot_min,index] = min(E_tot1);
        alpha_min_rad(k) = m_alpha(index);
        alpha_min(k) = m_alpha(index)*180/pi;
        Height_min(k) = Height(index);
        h_eq0(k) = z_free(index);
        z_eq0 = z_free(index);
        CofS = z_free(index)-R_cap(index)*cos(m_alpha(index));
        x_all{k} = r_all{index};
        y_all{k} = z_all{index};
        jiajiao_all{k} = psi_all{index};
        E_all{k}=E_tot1;
        alpha_all{k} = alpha_c;
        WrappingDegree(k,:) = [xi,rho_g,N_p,N,m_alpha(index),m_alpha(index)*180/pi,H_brush(index),...
            Height(index),z_eq0,CofS,E_tot_min,kappa,sigma,gamma,c_0,R_s,f,index];
        phaseDiagram(k,:) = [xi,rho_g,N_p,N,m_alpha(index)/(alpha_c(end)+0.01),m_alpha(index),m_alpha(index)*180/pi,H_brush(index),index];
        leg{k}=strcat('\xi=',num2str(xi),'nm');
    end
end
x = r_all{index};
y = z_all{index};
h = y(:,1);
alpha = alpha_c(index);
Rp = R_cap(index);
alpha_nd = alpha_c/(alpha_c(end)+0.01);
alpha_nd = alpha_nd';

subplot(2,2,1)
memColor = [1.0  0.5  0.0];
coatColor = [0.54 0.17 0.89];
theta2 = linspace(pi/2-alpha,pi/2+alpha,1001);
cir_x2 = Rp * cos(theta2);
cir_y2 = Rp * sin(theta2) - Rp * cos(alpha) + h;
plot(x,y,'Color', memColor,'Linewidth',4)
hold on;
plot(-x,y,'Color', memColor,'Linewidth',4)
hold on
plot(cir_x2, cir_y2,'Color', coatColor, 'LineWidth',4)
xlabel('r(nm)')
ylabel('z(nm)')
xlim([-200 200])
ylim([-10 110])
axis equal

subplot(2,2,2)
for k=1:length(xi_span)
    plot(alpha_all{k},E_all{k},'linewidth',2)
    hold on
end
    hold off
    legend(leg)