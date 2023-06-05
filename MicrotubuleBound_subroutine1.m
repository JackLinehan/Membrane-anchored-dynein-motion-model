function [rho_pull_vector] = MicrotubuleBound_subroutine1(theta_pull,length_data, r)
%% Microtubule Pulling Simulation 
[local_count,~] = size(length_data); 
%dt= 0.027; 
dt = 1; 
for j = 1:local_count % number of replicates 
%% 
dot_product = []; ds = []; x=[]; y=[]; z=[]; 
theta = theta_pull; % altitude 

theta = deg2rad(theta);

phi = -180:.1:180; %azimuth

phi = deg2rad(phi); 

motor_v = .02; % size of complex %0.19

theta_prob = makedist('Normal','mu',0,'sigma',max(theta).*0.34); 

phi_prob = makedist('Uniform',min(phi),max(phi)); 

% x_c = motor_v.*sin(deg2rad(90)).*cos(phi); 
% 
% y_c = motor_v.*sin(deg2rad(90)).*sin(phi);  
% hold on 
% scatter(x_c,y_c,'r'); 
% hold off 

runs = 2; 
%r = .09; %.15

theta_keep(1,1) = abs(random(theta_prob));
phi_keep(1,1) = random(phi_prob);

while (runs <= length_data(j,1)) 
    
            %x = r*sin(theta_keep(runs-1,1))*cos(phi_keep(runs-1,1)); 
            theta_loop = abs(random(theta_prob));
            phi_loop = random(phi_prob); 
            
            theta_keep(runs,1) = theta_loop + (0-theta_keep(runs-1,1)) + theta_keep(runs-1,1); 
            %phi_keep(runs,1) = acos(x/(r*sin(theta_keep(runs,1)))) + phi_keep(runs-1,1); 
            phi_keep(runs,1) = phi_loop + phi_keep(runs-1,1); 
            
            %theta_keep(runs,1) = theta_loop*dt + theta_keep(runs-1,1); 
            %phi_keep(runs,1) = phi_loop*dt + phi_keep(runs-1,1); 
            
            while (theta_keep(runs,1) > theta)
                       
                theta_keep(runs,1) = theta - 0.017; 
            end 
            
            while (theta_keep(runs,1) < -theta) 
                
                theta_keep(runs,1) = -theta + 0.017;
                
            end 
            
            
             
            %ds(runs-1,1) = sqrt(r.^2.*sin(theta_keep(runs-1,1)).^2.*phi_loop.^2 + r.^2.*theta_loop.^2); 
%                      x(runs-1,1) = r*sin(theta_loop)*cos(phi_loop); 
%                      y(runs-1,1) = r*sin(theta_loop)*sin(phi_loop); 
%                      z(runs-1,1) = r*cos(theta_loop); 
                     %dot_product(runs-1,1) = dot([x y z],[1 1 0]); 
                     
%                      hold on 
%                      scatter3(x,y,z,'*'); 
%                      hold off
%                      pause(.1); 
    runs = runs + 1; 
    
end 

%% Map Cartesian Coordinates to Fluorescent Intensity 

            x = r*sin(theta_keep).*cos(phi_keep); 
            y = r*sin(theta_keep).*sin(phi_keep); 
            z = r*cos(theta_keep); 
            
%             hold on 
%             scatter(x_c,y_c,'.');  
%             scatter3(x,y,z,'*'); 
%             hold off 


  rho = sqrt(diff(x).^2 + diff(y).^2);  

    x = []; y = []; z = []; ds=[]; dot_product = []; phi_keep = []; theta_keep = []; 

rho_pull_vector{j,1} = rho; 

end 
%% End of Simulation 

end 
