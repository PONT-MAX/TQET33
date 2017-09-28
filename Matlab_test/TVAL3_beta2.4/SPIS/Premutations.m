%% Code generated permutation of the Walsh Hadamar Matrix

% How to scramble the cols for three different sizes. Once run the same
% premutations will be used for the whole project.

N = 64^2;
col_perm = randperm(N); % column permutations allowable
save('Constants/col_perm_64.mat','col_perm');

N = 128^2;
col_perm = randperm(N); % column permutations allowable
save('Constants/col_perm_128.mat','col_perm');

N = 256^2;
col_perm = randperm(N); % column permutations allowable
save('Constants/col_perm_256.mat','col_perm');

%%
N = 512^2;
col_perm = randperm(N); % column permutations allowable
save('Constants/col_perm_512.mat','col_perm');

%%

% How to scramble the rows for three different sizes. Once run the same
% premutations will be used for the whole project. User have to chose how
% many rows to use N*ratio = M from this premutations.

N = 64^2;
row_perm = randperm(N-1) + 1; % Cant have row 1
save('Constants/row_perm_64','row_perm');

N = 128^2;
row_perm = randperm(N-1) + 1; % Cant have row 1
save('Constants/row_perm_128','row_perm');

N = 256^2;
row_perm = randperm(N-1) + 1; % Cant have row 1
save('Constants/row_perm_256','row_perm');

%%
N = 512^2;
row_perm = randperm(N-1) + 1; % Cant have row 1
save('Constants/row_perm_512','row_perm');