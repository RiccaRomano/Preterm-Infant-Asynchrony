%% Analyze the results from Morris screening
clear; clc; close all;
%%
% Load in the Morris results
load('Morris_results_VE.mat');
% For plotting, it can be useful to put all the parameter variable names
% here
Names={'TLC','RV','p_{l,opn}','p_{l,ran}','\beta',...
    '\gamma','k','v_{rc,str}','v_{ab,str}','p_{rc,sl}',...
    'RR','f','T','p_{frac}','A_{pia}',...
    'A_{mus}','R_{sm}','R_{sd}','K_s','R_{um}',...
    'K_u','I_u','p_{c,max}','p_{c,ran}','K_c',...
    'v_{c,max}','C_{ve}','R_{ve}','R_{rc}','R_{ab}'};


% Define these two values or load them from your Morris results
p = length(par_ids);  % NUMBER OF PARAMETERS (typically the size of 'par_ids' if you have a subset)
r = size(d,1);        % NUMBER OF MORRIS SAMPLES (first dimension of 'd' cell)

% If you specified a cell of names, you can use this line to help with
% plotting
par_names = Names(par_ids);
%% If you have multiple outputs, pick which outputs you want to use
%  when computing your morris indices, e.g., d{i}(:,ids);
ids = [1]; 

% Construct a vector for average (mu) and sample variance (s2)
% Note that this will store the values for all parameters (first index) AND
% all outputs considered in the sensitivity analysis (second index)
mu_all     = zeros(length(par_ids),length(ids));
s2_all     = zeros(length(par_ids),length(ids));
mustar_all = zeros(length(par_ids),length(ids));
output_counter = 1; %Index for storing each outputs sensitivity
% Now loop over all the outputs you want Morris effects on (e.g., pressure,
% flow)
for ii=ids
    mu_star = zeros(1,p); mu = zeros(1,p);
    s2 = zeros(1,p);
    
    % NOTE: if using a PDE or DAE system, you might have system crashes, so
    % be sure to account for this. r_counter incremenets the number of
    % successful solutions
    r_counter = 0; 

    % Inner loop over number of parameters (p) and number of samples (r)
    for i=1:p
        for j=1:r
            if ~isempty(d{j,i}) % Only get non-crash solutions
                % If you have a time dependent output, specify how you want
                % to construct d. Here is an example using the 2-norm

%                 mu_star(i) = mu_star(i) + norm(abs(d{j,i}(ii,:)),2); 
%                 mu(i)      = mu(i) + norm(d{j,i}(ii,:),2);

                % If your output is static, use this (can also replace with
                % max, min, mean of signal)
                mu_star(i) = mu_star(i) + abs(d{j,i}(ii,:)); 
                mu(i)      = mu(i) + d{j,i}(ii,:);

                % Increment for successful solutions
                r_counter = r_counter+1;
            end
        end
        % Now calculate the mean from the above loop
        mu_star(i) = (1./r_counter).*mu_star(i); 
        mu(i)      = (1./r_counter).*mu(i);
    end

    for i=1:p
        for j=1:r
            if ~isempty(d{j,i})
                % If you have a time dependent output, specify how you want
                % to construct d. Here is an example using the 2-norm
%                 s2(i) = s2(i) + (norm(d{j,i}(ii,:),2) - mu(i)).^2;
                
                % For scalar, see above as before
                s2(i) = s2(i) + (d{j,i}(ii,:) - mu(i)).^2;
            end
        end
        % Calculate the sample variance
        s2(i) = s2(i)./(r_counter-1);
    end

    % This is a scalar representation of parameter importance (See Colebank
    % and Chesler 2022 for discussion on this index.
    morris_rank = sqrt(mu_star.^2+s2);
    
    %Calculate average morris rank value
    avgmor=mean(morris_rank);
    median_morris=median(morris_rank);

    %Find and save influential parameters
    idxparVE=find(morris_rank>avgmor);
    infparVE=morris_rank(idxparVE);
    infnamesVE=par_names(idxparVE);
    idxparVE=[];
    for l=1:length(infnamesVE)
        index=find(strcmp(Names,infnamesVE{l}));
        idxparVE=[idxparVE index];
    end
    infnamesVE=Names(idxparVE);
    currentFile = sprintf('Inf_Pars_VE.mat');
    save(currentFile,"idxparVE","infnamesVE","infparVE");
    save VEMorris.mat
    
    % Now plot mu* vs s2
    figure(ii);
    set(gcf,'Position',[1 1 500 500]);
    loglog(mu_star,s2,'.','MarkerSize',12);
    % If you provided parameter names in the top of the script use this
    text(mu_star*1.1,s2.*1.1,par_names,'FontSize',12);
    set(gca,'FontSize',20);
    ylabel('Morris Variance');
    xlabel('Morris Mean');
    title('Minute Volume Morris Results');
    grid on;
    saveas(gcf,'VE_Morris_Results.fig');

    % Plot ranking
    fig=figure(1000+ii);
    set(gcf,'Position',[501 50 1350 500]);
    [rank_val,rank_loc] = sort(morris_rank,'descend');
    semilogy(rank_val,'ok','LineWidth',3,'MarkerSize',8);
    tex=annotation('textbox','String','B');
    tex.Parent=fig.CurrentAxes;
    tex.Position=[0 .2 3 0.1];
    tex.FontSize=24;
    tex.FontWeight='bold';
    tex.LineStyle="none";

    yline(avgmor,'LineWidth',2);
    xticks(1:length(par_ids));
    xticklabels(par_names(rank_loc)); % If you provided Names, use this line
    xtickangle(-30);
    set(gca,'FontSize',20);
    ylabel('Morris Indices for $\dot{V}_E$','interpreter','latex',"FontSize",30);
    grid on;
    saveas(gcf,'VE_Morris_Indices.fig');
    saveas(gcf,'VE_Morris_Indices.eps');
    saveas(gcf,'VE_Morris_Indices.svg');
    % Store the sensitivity for output ii
    mu_all(:,output_counter) = mu;
    s2_all(:,output_counter) = s2;
    mustar_all(:,output_counter) = mu_star;
    output_counter = output_counter+1;
    
end