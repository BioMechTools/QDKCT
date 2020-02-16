function [k_list] = convhull_list(k,n_length)
%convhull_list is aimed to convert the k from convhull
%Input:
%   k: k from convhull
%   n_length: total size of the vertices
n_len_k = length(k);
k_list = cell(n_len_k-1,1);
% min_k = min(k);
% max_k = max(k);
if k(3)-k(2) > 0
    for num = 1:n_len_k-2
        k_list{num} = k(num):1:k(num+1);
    end
    k_list{n_len_k-1} = [k(n_len_k-1):n_length 1:k(1)];  
else
    for num = 2:n_len_k-1
        k_list{num} = k(num+1):1:k(num);
    end
    k_list{1} = [k(2):n_length 1:k(1)]; 
end
	
end

