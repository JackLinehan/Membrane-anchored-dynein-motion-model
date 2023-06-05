tic(); 


%% Prophase measured data 
%measured_steps = 20; 
       % Anterior 
        load('cull_traj_less5_dynein_steps.mat')
%         idx = []; 
%         [idx,~] = find(steps<=0.14); 
        [h1, step_edges] = histcounts(steps); 
        x = h1./sum(h1); 
        measured_probability_anterior_prophase = x; 
%         AP_E = sum(measured_probability_anterior_prophase.*step_edges(1:end-1)); 
%         AP_Q = prctile(steps,75) - prctile(steps,25); 
%         AP_std = std(log(steps))^2; 
       measured_steps = step_edges; 
        
        % Posterior 
        load('posterior_prophase_culled_6step_trajectories.mat')
%           idx = []; 
%         [idx,~] = find(steps<=0.14); 
         [h1, step_edges] = histcounts(steps, measured_steps); 
        x = h1./sum(h1); 
        measured_probability_posterior_prophase = x; 
%         PP_E = sum(measured_probability_posterior_prophase.*step_edges(1:end-1)); 
%         PP_Q = prctile(steps,75) - prctile(steps,25); 
%         PP_std = std(sqrt(steps))^2; 
%% Anaphase measured data 
load('Anaphase_Anterior_Steps.mat')
r_local = cell2mat(r); clear vars r
%   idx = []; 
%         [idx,~] = find(r<=0.14); 
        [h1, step_edges] = histcounts(r_local, measured_steps); 
        x = h1./sum(h1); 
        measured_probability_anterior_anaphase = x; 
%         AA_E = sum(measured_probability_anterior_anaphase.*step_edges(1:end-1)); 
%         AA_Q = prctile(steps,75) - prctile(steps,25); 
%         AA_std = std(sqrt(steps))^2; 
load('posterior_anaphase_steps.mat')
%   idx = []; 
%         [idx,~] = find(steps<=0.14); 
        [h1, step_edges] = histcounts(steps, measured_steps); 
        x = h1./sum(h1); 
        measured_probability_posterior_anaphase = x; 
%         PA_E = sum(measured_probability_posterior_anaphase.*step_edges(1:end-1)); 
%         PA_Q = prctile(steps,75) - prctile(steps,25); 
%         PA_std = std(log(steps))^2; 
        


%% Simulation Runs 

fg_length = 0.31:0.01:.35; 
rho_analysis = zeros(3,length(fg_length)); 

for outer_loop = 1:length(fg_length) 
stable_sim = 1; 
while(stable_sim <=1) 
%% Anterior Prophase 
theta_bound = 10:5:90; 
%load('posterior_prophase_lengths.mat')
load('Aprophase_lengths.mat');
trajectory_length = trajectory_length'; 
for theta_pull_specification = 1:length(theta_bound)
        
    for dynein_pop = 1:length(trajectory_length)
    
        rho_data_vector1 = MicrotubuleBound_subroutine1(theta_bound(theta_pull_specification),trajectory_length{dynein_pop,1}, fg_length(1,outer_loop)); 
        
        rho_data_vector{theta_pull_specification,1}{dynein_pop,1} = cell2mat(rho_data_vector1); 
        
        clear vars rho_data_vector1
    end 
   
end
toc(); 
   
rho_data_vector = cellfun(@(x) cell2mat(x),rho_data_vector,'UniformOutput',false); 
%% Posterior Prophase 

theta_bound = 10:5:90; 
load('posterior_prophase_lengths.mat')
%trajectory_length = trajectory_length'; 
for theta_pull_specification = 1:length(theta_bound)
        
    for dynein_pop = 1:length(trajectory_length)
    
        rho_data_vector1pp = MicrotubuleBound_Slave1(theta_bound(theta_pull_specification),trajectory_length{dynein_pop,1},fg_length(1,outer_loop)); 
        
        rho_data_vector_pp{theta_pull_specification,1}{dynein_pop,1} = cell2mat(rho_data_vector1pp); 
        
        clear vars rho_data_vector1pp
    end 
   
end


rho_data_vector_pp = cellfun(@(x) cell2mat(x),rho_data_vector_pp,'UniformOutput',false); 

%% Anterior Anaphase 

theta_bound = 10:5:90; 
load('anterior_anaphase_length.mat')
%trajectory_length = trajectory_length'; 
for theta_pull_specification = 1:length(theta_bound)
        
    for dynein_pop = 1:length(trajectory_length)
    
        rho_data_vector1aa = MicrotubuleBound_Slave1(theta_bound(theta_pull_specification),trajectory_length{dynein_pop,1},fg_length(1,outer_loop)); 
        
        rho_data_vector_aa{theta_pull_specification,1}{dynein_pop,1} = cell2mat(rho_data_vector1aa); 
        
        clear vars rho_data_vector1aa
    end 
   
end


rho_data_vector_aa = cellfun(@(x) cell2mat(x),rho_data_vector_aa,'UniformOutput',false); 

%% Posterior Anaphase 

theta_bound = 10:5:90; 
load('posterior_anaphase_length.mat')
%trajectory_length = trajectory_length'; 
for theta_pull_specification = 1:length(theta_bound)
        
    for dynein_pop = 1:length(trajectory_length)
    
        rho_data_vector1pa = MicrotubuleBound_Slave1(theta_bound(theta_pull_specification),trajectory_length{dynein_pop,1},fg_length(1,outer_loop)); 
        
        rho_data_vector_pa{theta_pull_specification,1}{dynein_pop,1} = cell2mat(rho_data_vector1pa); 
        
        clear vars rho_data_vector1pa
    end 
   
end


rho_data_vector_pa = cellfun(@(x) cell2mat(x),rho_data_vector_pa,'UniformOutput',false); 

        
%% Metrics 

for j = 1:length(theta_bound) % through theta bound
              
        for k = 1:j % through theta pull vector 
            
                
                        rho_theory_all = []; rho_theory_all_pp = []; rho_theory_all_aa = []; 
                        rho_theory_all_pa = []; 

         %% Floppy Dynein 
                rho_theory_all = [rho_theory_all; rho_data_vector{j,1}]; 
                
                rho_theory_all_pp = [rho_theory_all_pp; rho_data_vector_pp{j,1}]; 
                
                rho_theory_all_aa = [rho_theory_all_aa; rho_data_vector_aa{j,1}]; 
                
                rho_theory_all_pa = [rho_theory_all_pa; rho_data_vector_pa{j,1}]; 
              
        %% Pulling 
                rho_theory_all = [rho_theory_all; rho_data_vector{k,1}]; 
                
                rho_theory_all_pp = [rho_theory_all_pp; rho_data_vector_pp{k,1}]; 
                
                rho_theory_all_aa = [rho_theory_all_aa; rho_data_vector_aa{k,1}]; 
                
                rho_theory_all_pa = [rho_theory_all_pa; rho_data_vector_pa{k,1}]; 
           
%% Combined Theoretical Values 
                %rho_theory_all = rho_theory_all; 
        [x, ~] = histcounts(rho_theory_all,step_edges); 
        x_theory = x./sum(x); 
        
        %rho_theory_pp = rho_theory_pp; 
        [x_pp,~] = histcounts(rho_theory_all_pp,step_edges); 
        x_theory_pp = x_pp./sum(x_pp); 
        
        [x_aa,~] = histcounts(rho_theory_all_aa,step_edges); 
        x_theory_aa = x_aa./sum(x_aa); 
        
        [x_pa,~] = histcounts(rho_theory_all_pa,step_edges); 
        x_theory_pa = x_pa./sum(x_pa); 
        
%         subplot(1,2,1); 
%         hold on 
%         plot(step_edges(2:end),x_theory,'r'); 
%         plot(step_edges(2:end),measured_probability_anterior_prophase,'k'); 
%         hold off 
%         title('anterior prophase'); 
%         
%         subplot(1,2,2); 
%         hold on 
%         plot(step_edges(2:end),x_theory_pp,'r'); 
%         plot(step_edges(2:end),measured_probability_posterior_prophase,'k'); 
%         title('posterior prophase'); 
%         pause(); 
%         figure(clf);
        
        % anterior prophase 
        AP_error = sum((measured_probability_anterior_prophase - x_theory).^2);%/AP_std;
        %AP_error = rmse(x_theory, measured_probability_anterior_prophase); 
        % posterior prophase 
        PP_error = sum((measured_probability_posterior_prophase - x_theory_pp).^2);%/PP_std;
        %PP_error = rmse(x_theory, measured_probability_posterior_prophase);
        % anterior anaphase 
        AA_error = sum((measured_probability_anterior_anaphase - x_theory_aa).^2);%/AA_std;%/AA_Q;
        %AA_error = rmse(x_theory, measured_probability_anterior_anaphase); 
        % posterior anaphase
        %PA_error = rmse(x_theory, measured_probability_posterior_anaphase); 
        PA_error = sum((measured_probability_posterior_anaphase - x_theory_pa).^2);%/PA_std;%/PA_Q;
        
        metric{j}{k} = [AP_error, PP_error, AA_error, PA_error];  
        
        rho_theory_all = []; rho_theory_pp = []; 
        
%                     [prob_AP,edges1] = histcounts(metric{j}{k}{free_N,1}{m,1}(:,1)); 
%                     edges1 = movmean(edges1,2); edges1 = edges1(2:end); 
%                     [prob_PP,edges2] = histcounts(metric{j}{k}{free_N,1}{m,1}(:,2)); 
%                     edges2 = movmean(edges2,2); edges2 = edges2(2:end); 
%                     [prob_AA,edges3] = histcounts(metric{j}{k}{free_N,1}{m,1}(:,3)); 
%                     edges3 = movmean(edges3,2); edges3 = edges3(2:end); 
%                     [prob_PA, edges4] = histcounts(metric{j}{k}{free_N,1}{m,1}(:,4)); 
%                     edges4 = movmean(edges4,2); edges4 = edges4(2:end); 
%                     
%                     AP_E = sum(prob_AP.*edges1); PP_E = sum(prob_PP.*edges2); 
%                     AA_E = sum(prob_AA.*edges3); PA_E = sum(prob_PA.*edges4); 
%                     
%                     metric_var{j}{k}{free_N,1}{m} = [sum(prob_AP.*((edges1-AP_E).^2)) sum((prob_PP.*((edges2-PP_E).^2))) ...
%                         sum(prob_AA.*((edges3-AA_E).^2)) sum(prob_PA.*((edges4-PA_E).^2))]; 
                    
                    metric_mean{j}{k} = mean([AP_error, PP_error, AA_error, PA_error]); 
                    %metric_mean{j}{k}{m} = mean(metric{j}{k}{m,1}); 
        
                
            local_metric = []; 
            
        end 
        
end 

% for j = 1:length(metric_mean) 
%     keep_minimums = [];
%     for k = 1:length(metric_mean{j}) 
%          keep_error_param = []; 
%         for free_N = 1:length(metric_mean{j}{k})
%             
%                 
%                 %a = cellfun(@cell2mat, metric{j}, 'UniformOutput', false); 
%                 a = cell2mat(metric_mean{1,j}{1,k}{free_N,1}); % pull theta free set
% 
%                 [value,location] = min(a); 
%                 keep_error_param(1,free_N) = theta_bound(j); % theta free
%                 keep_error_param(2,free_N) = theta_bound(k); % theta interacting 
%                 disp('Maximum Angle of Deflection for Swivel Free Model is'); 
%                 disp(num2str(theta_bound(j)));  
%         
%                 disp('Maximum Angle of Deflecton for Swivel Interacting Model is'); 
%                 disp(num2str(theta_bound(k))); 
%             
%                 keep_error_param(3,free_N) = dynein_number_eval(free_N); 
%                 keep_error_param(4,free_N) = dynein_number_eval(location); 
%                 keep_error_param(5,free_N) = a(location); 
%         
%         end 
%         
%            b = keep_error_param(5,:); 
%           [value, location] = min(b); 
%           keep_minimums(:,k) = keep_error_param(:,location); 
%         
%     end 
%     
%     c = keep_minimums(5,:); 
%     [value,location] = min(c); 
%     all_keep_minimimums(:,j) = keep_minimums(:,location); 
%    
% end 
clear vars all_keep_minimums 
for j = 1:length(metric_mean) 
    
    a = cell2mat(metric_mean{j}); 
    
    [value,loc] = min(a); 
    
    all_keep_minimums(:,j) = [theta_bound(j); theta_bound(loc); value]; 
    
    clear vars value loc
    
end 

[val,loc] = min(all_keep_minimums(end,:)); 
        disp('%%%%%%%%%%%%%%%%%%%'); 
        disp(all_keep_minimums(:,loc)); 
        disp('%%%%%%%%%%%%%%%%%%')
         
        theta_free(stable_sim,1) = all_keep_minimums(1,loc); 
        theta_interacting(stable_sim,1) = all_keep_minimums(2,loc); 
        error_keeper(stable_sim,1) = all_keep_minimums(3,loc); 
       % clear vars rho_data_vector rho_data_vector1 metric metric_mean
       % clear vars rho_data_vector_pp rho_data_vector_pa rho_data_vector_aa
       clear vars rho_data_vector rho_data_vector_pp rho_data_vector_aa rho_data_vector_pa 
       clear vars rho_theory_all rho_theory_all_pp rho_theory_all_pa rho_theory_all_aa
       clear vars x_theory x_theory_pp x_theory_aa x_theory_pa 
       clear vars AP_error PP_error AA_error PA_error metric 
       clear vars all_keep_minimums 
       
       stable_sim = stable_sim + 1; 
end 

       rho_analysis(:,outer_loop) = [mean(theta_free); mean(theta_interacting); mean(error_keeper)]; 
       
end 
    
        
        
        