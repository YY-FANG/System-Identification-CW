clear;

%% Initialization
N = 2000; batch = 100; y = zeros(N,1); e = zeros(N,1);

a = 1/2; lambda = 3; mu = 1; y_0 = 0; w_0 = 0;

varphi_1 = zeros(N,1); varphi_2 = zeros(2,1,N);                            % known past quantities

theta_1 = zeros(batch,1); theta_2 = zeros(2,1,batch);                      % parameters of batch LS

phiphi_AR2 = zeros(2,2,N); phiy_AR2 = zeros(2,1,N);

theta_1_ = []; theta_2_ = []; color = 0;
 
%% Generation
for k = 1:1                                                                % number of attempts when batch=1
    for j = 1:batch     
        if color == 1
            for i = 1:N
                if i == 1
                    e(i) = randn(1);
                else
                    e(i) = -0.5*e(i-1)+randn(1);                           
                end
            end                                                            % generate the color noise of e
        else
            e = lambda.*randn(N,1);                                        % generate the white noise of e
        end 
        w = mu.*randn(N,1);                                                % generate the white noise of w
        for i = 1:N
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
        end
    
%% Least-square Solution                                                   
        theta_1(j) = (sum(varphi_1.*varphi_1))\sum(varphi_1.*y);           % AR(1)
        theta_2(:,:,j) = (sum(phiphi_AR2,3))\sum(phiy_AR2,3);              % AR(2)
                                                                           % sovle Ax=b by (A\b)
    end           

%% Empirical Mean and Variances of Parameters
    theta_1_mean = mean(theta_1); theta_2_mean = mean(theta_2,3);          % the empirical mean

    theta_1_ = [theta_1_ theta_1_mean]; theta_2_ = [theta_2_ theta_2_mean]; 

    sigma_11 = mean((theta_2(1,:,:)-theta_2_mean(1,:)).^2);

    sigma_12 = mean((theta_2(1,:,:)-theta_2_mean(1,:)).*(theta_2(2,:,:)-theta_2_mean(2,:)));

    sigma_21 = mean((theta_2(1,:,:)-theta_2_mean(1,:)).*(theta_2(2,:,:)-theta_2_mean(2,:)));

    sigma_22 = mean((theta_2(2,:,:)-theta_2_mean(2,:)).^2);

    theta_1_var = var(theta_1);  theta_2_var = [sigma_11 sigma_12; sigma_21 sigma_22];      % the variances  
end

AR1 = mean(theta_1_) 
AR2 = mean(theta_2_,2)