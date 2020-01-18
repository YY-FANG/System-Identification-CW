clear;

%% Initialization
N = 1000; batch =100; alpha_s = 0.01; alpha = 10; y = zeros(N,1); 

a = 1/2; lambda = 3; mu = 1; y_0 = 0; w_0 = 0;

theta_1 = zeros(batch,1); theta_2 = zeros(2,1,batch);                      % parameters of batch LS 

theta_S_1 = zeros(batch,1); theta_S_2 = zeros(2,1,batch);                  % parameters of recursive LS (1st form)

theta_R_1 = zeros(batch,1); theta_R_2 = zeros(2,1,batch);                  % parameters of recursive LS (2nd form)

theta_V_1 = zeros(batch,1); theta_V_2 = zeros(2,1,batch);                  % parameters of recursive LS (3rd form)

theta_V_1s = zeros(batch,1); theta_V_2s = zeros(2,1,batch);                % parameters of recursive LS (3rd form) 
                                                                           % with small alpha

%% Batch LS Initialization
varphi_1 = zeros(N,1); varphi_2 = zeros(2,1,N);                            % known past quantities

phiphi_AR2 = zeros(2,2,N); phiy_AR2 = zeros(2,1,N);

%% Recursive LS (1st form) Initialization 
theta_S1 = zeros(N,1); K_S1 = zeros(N,1);     

epsilon_S1 = zeros(N,1); S1 = zeros(N,1);                                  % AR1

theta_S2 = zeros(2,1,N); K_S2 = zeros(2,1,N);

epsilon_S2 = zeros(N,1); S2 = zeros(2,2,N);                                % AR2

%% Recursive LS (2nd form) Initialization 
theta_R1 = zeros(N,1); K_R1 = zeros(N,1); 

epsilon_R1 = zeros(N,1); R1 = zeros(N,1);                                  % AR1

theta_R2 = zeros(2,1,N); K_R2 = zeros(2,1,N); 

epsilon_R2 = zeros(N,1); R2 = zeros(2,2,N);                                % AR2

%% Recursive LS (3rd form) Initialization
V1_0 = alpha*1; theta_V1_0 = 0;

theta_V1 = zeros(N,1); K_V1 = zeros(N,1); 

epsilon_V1 = zeros(N,1); V1 = zeros(N,1);                                  % AR1

V2_0 = alpha*eye(2,2); theta_V2_0 = [0;0];

theta_V2 = zeros(2,1,N); K_V2 = zeros(2,1,N); 

epsilon_V2 = zeros(N,1); V2 = zeros(2,2,N);                                % AR2

%% Recursive LS (3rd form) Initialization with small alpha
V1_0s = alpha_s*1; theta_V1_0s = 0;

theta_V1s = zeros(N,1); K_V1s = zeros(N,1); 

epsilon_V1s = zeros(N,1); V1s = zeros(N,1);                                  % AR1

V2_0s = alpha_s*eye(2,2); theta_V2_0s = [0;0];

theta_V2s = zeros(2,1,N); K_V2s = zeros(2,1,N); 

epsilon_V2s = zeros(N,1); V2s = zeros(2,2,N);                                % AR2

%% Generation
for j = 1:batch
    e = lambda.*randn(N,1);                                                % generate the white noise of e
    w = mu.*randn(N,1);                                                    % generate the white noise of w
    
    for i = 1:N        
%% Batch LS
        if i == 1
            y(i) = a*y_0+e(i)+w(i)-a*w_0;
            varphi_1(i) = y_0; varphi_2(:,:,i) = [y_0;y_0];
        elseif i == 2
            y(i) = a*y(i-1)+e(i)+w(i)-a*w(i-1);
            varphi_1(i) = y(i-1); varphi_2(:,:,i) = [y(i-1);y_0];
        elseif i > 2
            y(i) = a*y(i-1)+e(i)+w(i)-a*w(i-1);
            varphi_1(i) = y(i-1); varphi_2(:,:,i) = [y(i-1);y(i-2)];
        end
        phiphi_AR2(:,:,i) = varphi_2(:,:,i)*varphi_2(:,:,i)';
        phiy_AR2(:,:,i) = varphi_2(:,:,i)*y(i);
        
%% Recursive LS (1st form)       
        if i>1
            S1(i) = S1(i-1)+varphi_1(i)*varphi_1(i)';
            epsilon_S1(i) = y(i)-varphi_1(i)'*theta_S1(i-1);
            K_S1(i) = S1(i)\varphi_1(i);
            theta_S1(i) = theta_S1(i-1)+K_S1(i)*epsilon_S1(i);             
        end                                                                % AR1 
        
        if i>2                                   
            S2(:,:,i) = S2(:,:,i-1)+varphi_2(:,:,i)*varphi_2(:,:,i)';
            epsilon_S2(i) = y(i)-varphi_2(:,:,i)'*theta_S2(:,:,i-1);
            K_S2(:,:,i) = S2(:,:,i)\varphi_2(:,:,i);
            theta_S2(:,:,i) = theta_S2(:,:,i-1)+K_S2(:,:,i)*epsilon_S2(i);      
        end                                                                % AR2 
        
%% Recursive LS (2nd form)        
        if i>1
            R1(i) = R1(i-1)+1/i*(varphi_1(i)*varphi_1(i)'-R1(i-1));
            epsilon_R1(i) = y(i)-varphi_1(i)'*theta_R1(i-1);
            K_R1(i) = 1/i*(R1(i)\varphi_1(i));
            theta_R1(i) = theta_R1(i-1)+K_R1(i)*epsilon_R1(i);
        end                                                                % AR1
            
        if i>2
            R2(:,:,i) = R2(:,:,i-1)+1/i*(varphi_2(:,:,i)*varphi_2(:,:,i)'-R2(:,:,i-1));
            epsilon_R2(i) = y(i)-varphi_2(:,:,i)'*theta_R2(:,:,i-1);
            K_R2(:,:,i) = 1/i*(R2(:,:,i)\varphi_2(:,:,i));
            theta_R2(:,:,i) = theta_R2(:,:,i-1)+K_R2(:,:,i)*epsilon_R2(i);
        end                                                                % AR2         

%% Recursive LS (3rd form) 
        if i == 1
            beta1 = 1+varphi_1(i)'*V1_0*varphi_1(i);
            V1(i) = V1_0-1/beta1*V1_0*varphi_1(i)*varphi_1(i)'*V1_0;
            epsilon_V1(i) = y(i)-varphi_1(i)'*theta_V1_0;
            K_V1(i) = V1(i)*varphi_1(i);
            theta_V1(i) = theta_V1_0+K_V1(i)*epsilon_V1(i);                % AR1
            
            beta2 = 1+varphi_2(:,:,i)'*V2_0*varphi_2(:,:,i);
            V2(:,:,i) = V2_0-1/beta2*V2_0*varphi_2(:,:,i)*varphi_2(:,:,i)'*V2_0;
            epsilon_V2(i) = y(i)-varphi_2(:,:,i)'*theta_V2_0;
            K_V2(:,:,i) = V2(:,:,i)*varphi_2(:,:,i);
            theta_V2(:,:,i) = theta_V2_0+K_V2(:,:,i)*epsilon_V2(i);        % AR2
        else
            beta1 = 1+varphi_1(i)'*V1(i-1)*varphi_1(i);
            V1(i) = V1(i-1)-1/beta1*V1(i-1)*varphi_1(i)*varphi_1(i)'*V1(i-1);
            epsilon_V1(i) = y(i)-varphi_1(i)'*theta_V1(i-1);
            K_V1(i) = V1(i)*varphi_1(i);
            theta_V1(i) = theta_V1(i-1)+K_V1(i)*epsilon_V1(i);             % AR1
            
            beta2 = 1+varphi_2(:,:,i)'*V2(:,:,i-1)*varphi_2(:,:,i);
            V2(:,:,i) = V2(:,:,i-1)-1/beta2*V2(:,:,i-1)*varphi_2(:,:,i)*varphi_2(:,:,i)'*V2(:,:,i-1);
            epsilon_V2(i) = y(i)-varphi_2(:,:,i)'*theta_V2(:,:,i-1);
            K_V2(:,:,i) = V2(:,:,i)*varphi_2(:,:,i);
            theta_V2(:,:,i) = theta_V2(:,:,i-1)+K_V2(:,:,i)*epsilon_V2(i); % AR2
        end
        
%% Recursive LS (3rd form) with small alpha
        if i == 1
            beta1s = 1+varphi_1(i)'*V1_0s*varphi_1(i);
            V1s(i) = V1_0s-1/beta1s*V1_0s*varphi_1(i)*varphi_1(i)'*V1_0s;
            epsilon_V1s(i) = y(i)-varphi_1(i)'*theta_V1_0s;
            K_V1s(i) = V1s(i)*varphi_1(i);
            theta_V1s(i) = theta_V1_0s+K_V1s(i)*epsilon_V1s(i);            % AR1
            
            beta2s = 1+varphi_2(:,:,i)'*V2_0s*varphi_2(:,:,i);
            V2s(:,:,i) = V2_0s-1/beta2s*V2_0s*varphi_2(:,:,i)*varphi_2(:,:,i)'*V2_0s;
            epsilon_V2s(i) = y(i)-varphi_2(:,:,i)'*theta_V2_0s;
            K_V2s(:,:,i) = V2s(:,:,i)*varphi_2(:,:,i);
            theta_V2s(:,:,i) = theta_V2_0s+K_V2s(:,:,i)*epsilon_V2s(i);    % AR2
        else
            beta1s = 1+varphi_1(i)'*V1s(i-1)*varphi_1(i);
            V1s(i) = V1s(i-1)-1/beta1s*V1s(i-1)*varphi_1(i)*varphi_1(i)'*V1s(i-1);
            epsilon_V1s(i) = y(i)-varphi_1(i)'*theta_V1s(i-1);
            K_V1s(i) = V1s(i)*varphi_1(i);
            theta_V1s(i) = theta_V1s(i-1)+K_V1s(i)*epsilon_V1s(i);         % AR1
            
            beta2s = 1+varphi_2(:,:,i)'*V2s(:,:,i-1)*varphi_2(:,:,i);
            V2s(:,:,i) = V2s(:,:,i-1)-1/beta2s*V2s(:,:,i-1)*varphi_2(:,:,i)*varphi_2(:,:,i)'*V2s(:,:,i-1);
            epsilon_V2s(i) = y(i)-varphi_2(:,:,i)'*theta_V2s(:,:,i-1);
            K_V2s(:,:,i) = V2s(:,:,i)*varphi_2(:,:,i);
            theta_V2s(:,:,i) = theta_V2s(:,:,i-1)+K_V2s(:,:,i)*epsilon_V2s(i); % AR2
        end        

    end
    
%% Least-square Solution: Batch LS
    theta_1(j) = (sum(varphi_1.*varphi_1))\sum(varphi_1.*y);               % AR(1)
    theta_2(:,:,j) = (sum(phiphi_AR2,3))\sum(phiy_AR2,3);                  % AR(2)

%% Least-square Solution: Recursive LS (1st form)
    theta_S_1(j) = theta_S1(N);                                            % AR(1)
    theta_S_2(:,:,j) = theta_S2(:,:,N);                                    % AR(2)

%% Least-square Solution: Recursive LS (2nd form)
    theta_R_1(j) = theta_R1(N);                                            % AR(1)
    theta_R_2(:,:,j) = theta_R2(:,:,N);                                    % AR(2)    
    
%% Least-square Solution: Recursive LS (3rd form)
    theta_V_1(j) = theta_V1(N);                                            % AR(1)
    theta_V_2(:,:,j) = theta_V2(:,:,N);                                    % AR(2) 
    
%% Least-square Solution: Recursive LS (3rd form) with small alpha
    theta_V_1s(j) = theta_V1s(N);                                          % AR(1)
    theta_V_2s(:,:,j) = theta_V2s(:,:,N);                                  % AR(2)      
    
end

%% Empirical Mean of Parameters
theta_1_mean = mean(theta_1); theta_2_mean = mean(theta_2,3);              % the empirical mean of batch LS

theta_S_1_mean = mean(theta_S_1); theta_S_2_mean = mean(theta_S_2,3);      % the empirical mean of recursive LS (1st form)

theta_R_1_mean = mean(theta_R_1); theta_R_2_mean = mean(theta_R_2,3);      % the empirical mean of recursive LS (2nd form)

theta_V_1_mean = mean(theta_V_1); theta_V_2_mean = mean(theta_V_2,3);      % the empirical mean of recursive LS (3rd form)

theta_V_1s_mean = mean(theta_V_1s); theta_V_2s_mean = mean(theta_V_2s,3);  % the empirical mean of recursive LS (3rd form) 
                                                                           % with small alpha
AR1 = theta_1_mean
AR1_S = theta_S_1_mean
AR1_R = theta_R_1_mean
AR1_V = theta_V_1_mean
AR1_Vs = theta_V_1s_mean

AR2 = theta_2_mean
AR2_S = theta_S_2_mean
AR2_R = theta_R_2_mean
AR2_V = theta_V_2_mean
AR2_Vs = theta_V_2s_mean

%% Plot for Comparison
AR1_theory = zeros(N,1); AR1_theory(:) = 0.4615; x = 1:N;
R_bar = zeros(N,1); R_bar(:) = 13; 
V_bar = zeros(N,1); V_bar(:) = 0; 
figure(1)
subplot(3,1,1);
plot(x,AR1_theory,'Color','#0072BD'); hold on;
plot(x,theta_S1,'Color','#D95319');
legend('theoretical','recursive LS (S)');
xlabel('N'); ylabel('Value of Parameter');
axis([1 1000 -1 1.5]); 
subplot(3,1,2)
plot(x,AR1_theory,'Color','#0072BD'); hold on;
plot(x,theta_R1,'Color','#77AC30');
legend('theoretical','recursive LS (R)');
xlabel('N'); ylabel('Value of Parameter');
axis([1 1000 -1 1.5]); 
subplot(3,1,3)
plot(x,AR1_theory,'Color','#0072BD'); hold on;
plot(x,theta_V1,'Color','#A2142F'); 
legend('theoretical','recursive LS (V)');
xlabel('N'); ylabel('Value of Parameter');
axis([1 1000 -1 1.5]); 

figure(2)
subplot(2,1,1);
plot(x,AR1_theory,'Color','#0072BD'); hold on;
plot(x,theta_V1,'Color','#A2142F'); 
legend('theoretical','recursive LS (V)');
xlabel('N'); ylabel('Value of Parameter');
axis([1 1000 -1 1.5]); 
subplot(2,1,2);
plot(x,AR1_theory,'Color','#0072BD'); hold on;
plot(x,theta_V1s,'Color','#7E2F8E'); 
legend('theoretical','recursive LS (V) with small alpha');
xlabel('N'); ylabel('Value of Parameter');
axis([1 1000 -1 1.5]); 

figure(3)
plot(x,S1,'Color','#D95319');
xlabel('N'); ylabel('S(t)');

figure(4)
plot(x,R_bar,'Color','#0072BD'); hold on;
plot(x,R1,'Color','#77AC30');
xlabel('N'); ylabel('R(t)');
axis([1 1000 0 16]); 

figure(5)
plot(x,V1,'Color','#A2142F');
xlabel('N'); ylabel('V(t)');
axis([1 1000 0 0.01]); 
