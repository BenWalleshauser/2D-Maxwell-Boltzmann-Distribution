%Ben Walleshauser
%10/15/22

%% Initialization
clc
clear
close all

%Number of particles
N = 16;

%Temporal domain
h = 0.01;
t = 0:h:5000;

%Size of grid (dimensionless units of r/sigma)
Lx = 30;   
Ly = 30;

%Particle Positions (all start on LHS)
X = zeros(N,2,length(t));
X(:,:,1) = [repelem([4/5:1.2*4/5:1.2*16/5],4)' repmat([2:1.5*2:1.5*8],[1 4])'];


%Particle Velocities
vmax = 0;
V = zeros(N,2,length(t));
V(:,:,1) = sqrt(2)*vmax.*(rand(N,2)-0.5);

%Particle Accelerations
A = zeros(N,2,length(t));

%[x,y]=Periodic(-8.1,-8.1,Lx,Ly)
%% Simualtion

%Initial Accelerations
k = 1;
for i = 1:N
    
    for j = 1:N

        %Skip if particle index same
        if i == j
            continue
        end

        del_x = X(j,1,k)-X(i,1,k);    
        del_x = Return_dX(del_x,Lx,Ly);

        del_y = X(j,2,k)-X(i,2,k);    
        del_y = Return_dX(del_y,Lx,Ly);

        r = sqrt(del_x^2 + del_y^2);
        Fmag = -24*(1/r^7 -2/r^13);

        %Newtons 2nd Law (m = unity)
        A(i,1,k) = A(i,1,k) + Fmag*(-del_x/r);
        A(i,2,k) = A(i,2,k) + Fmag*(-del_y/r);

    end
end

%Verlet algorithm
for k = 1:length(t)-1
    
    %Updated positions
    for i = 1:N
        X(i,1,k+1) = X(i,1,k) + V(i,1,k)*h + 0.5*A(i,1,k)*h^2;
        X(i,2,k+1) = X(i,2,k) + V(i,2,k)*h + 0.5*A(i,2,k)*h^2;
        
        [X(i,1,k+1), X(i,2,k+1)] = Periodic(X(i,1,k+1),X(i,2,k+1),Lx,Ly);
    end
    
    %Updated accelerations
    for i = 1:N
        
        %Effect of jth particle
        for j = 1:N
            
            %Skip if particle index same
            if i == j
                continue
            end
            
            del_x = X(j,1,k+1)-X(i,1,k+1);    
            del_x = Return_dX(del_x,Lx,Ly);

            del_y = X(j,2,k+1)-X(i,2,k+1);    
            del_y = Return_dX(del_y,Lx,Ly);
            
            r = sqrt(del_x^2 + del_y^2);
            Fmag = -24*(1/r^7 -2/r^13);
            
            %Newtons 2nd Law (m = unity)
            A(i,1,k+1) = A(i,1,k+1) + Fmag*(-del_x/r);
            A(i,2,k+1) = A(i,2,k+1) + Fmag*(-del_y/r);
            
        end             
    end
    
    %Updated Velocities
    for i = 1:N
        V(i,1,k+1) = V(i,1,k) + 0.5*(A(i,1,k) + A(i,1,k+1))*h;
        V(i,2,k+1) = V(i,2,k) + 0.5*(A(i,2,k) + A(i,2,k+1))*h;
    end
end

%% Magnitude

%Initial configuration is at equilubrium
twarmup = 100;
ti_start = twarmup/h+1;

V_eq = V(:,:,ti_start:end);

%Computing velocity magnitudes
for i = 1:length(V_eq(1,1,:))
   V_eq_mag(:,i) = sqrt(V_eq(:,1,i).^2 + V_eq(:,2,i).^2);
end

%Calculating kbT
kbT = (1/N).*sum(0.5.*mean(V_eq_mag(:,i).^2,2))

delV = 0.1*sqrt(kbT);

Veq_max = max(max(V_eq_mag));
bin_lower = 0:delV:Veq_max;
bin_upper = bin_lower + delV;
V_eq_mag_vec = V_eq_mag(:);
probk = zeros(length(bin_lower),1);
hist=histogram(V_eq_mag_vec,bin_lower);
probk = hist.Values;

%Normazlizing
a = 1/(delV*sum(probk));
probk_normalized = a.*probk;

%Averages
mean(probk_normalized)
avg_v = mean(V_eq_mag_vec);

%Theoretical
v = 0:0.01:Veq_max;
f = (1/(kbT)).*v.*exp(-v.^2/(2*kbT));

%Plotting Normalized Histogram
figure(3)
hold on
histogram(V_eq_mag_vec,bin_lower,'Normalization','pdf')
plot(v,f)
xlabel('Velocity Magnitude V')
ylabel('P(V)')
title('Probability Distribution of Velocity Magnitudes')
grid on
legend('Numerical','Theoretical')
hold off


%%
function dXnew = Return_dX(del_x,Lx,Ly)
        if abs(del_x) > 0.5*Lx
           if del_x > 0
               dXnew = -Lx + del_x;
           elseif del_x < 0
               dXnew = Lx + del_x;
           end 
        else
            dXnew = del_x;
        end
end

function [xnew, ynew] = Periodic(xold,yold,Lx,Ly)
        %if diverge...
        if abs(xold) > 1e5 | abs(yold) > 1e5
            disp('Model Divergence!')
            return
        end
        
        while xold > Lx
            xold = xold-Lx;
        end
        while yold > Ly
            yold = yold-Ly;
        end        
        while xold < 0
            xold = xold+Lx;
        end
        while yold < 0
            yold = yold+Ly;
        end
        xnew = xold;
        ynew = yold;
end
