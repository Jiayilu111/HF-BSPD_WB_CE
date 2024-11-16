function [alpha_pattern, sin_theta_pattern] = generate_BiLinear(f, fc, alpha_list, sin_theta_list)
%%Minimum error index
    eta = f ./ fc;
    M = length(eta);
    S = length(alpha_list);
    D = length(sin_theta_list);
    alpha_pattern = zeros(S, M);
    sin_theta_pattern = zeros(D, M);
    %%
    alpha_value = alpha_list * eta;
%     sin_theta_value = mod(sin_theta_list ./ eta + 1, 2 ) - 1;
    sin_theta_value = sin_theta_list * eta ;
    for i_s = 1:S
       for i_m = 1:M
          temp = abs(alpha_value(i_s, i_m) - alpha_list);
          [~, alpha_pattern(i_s, i_m)] = min(temp); % approximates the optimal sample
       end
    end
    
    for i_d = 1:D
       for i_m = 1:M
          temp = abs(sin_theta_value(i_d, i_m) - sin_theta_list);
          [~, sin_theta_pattern(i_d, i_m)] = min(temp);
       end
    end
    
end

