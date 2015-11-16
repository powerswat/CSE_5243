function kMinHash(tot_mat, tot_vec_lbl)

%% Generate an input matrix (# of shingles * # of docs)
input_mat = tot_mat';

%% (Opt 1): Approx. Linaer Permutation Hashing
% a = randperm(10,1);
% b = randperm(10,1)+1;
% N = size(input_mat,1);
% tmp_ps = primes(N+100);
% p = tmp_ps(max(find(tmp_ps > N)));
% h_ab = mod(mod(a.*input_mat(:,1)+ b, p), N);

%% (Opt 2): Linear Permutation Hashing


end