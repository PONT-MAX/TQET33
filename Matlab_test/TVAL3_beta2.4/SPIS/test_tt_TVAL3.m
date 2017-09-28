clear all; close all;
% problem size
sidelength = 128;
N = sidelength^2;
fullscreen = get(0,'ScreenSize')*0.7;

% set the optional paramaters
clear opts
opts.mu = 2^8; opts.beta = 2^7;
opts.maxcnt = 150; opts.tol_inn = 1e-3;
opts.tol = 1E-6; opts.maxit = 750;
opts.mu0 = 2^4; opts.beta0 = 2^0;    % trigger continuation scheme
opts.nonneg = true; opts.isreal = true

% permutation to WH-matrix
load(['row_perm_', int2str(sidelength)]);
load(['col_perm_' int2str(sidelength)]);

% 


for mu = [2^6 2^8 2^10]
    opts.mu = mu;
    for i = 0:15
        % Importera SWIR-bild
        Im0 = imread(['../im/Shapes/',int2str(i),'.PNG']);
        Im0 = double(Im0(:,:,1));
        for ratio = [0.1 0.2 0.3]
            M = round(N*ratio);
            picks = row_perm(1:M);
            
            % Create Sensing matrix
            A_1 = double(fWHtrans(eye(N))*N);
            % Sensing matrix 1
            A_1 = (A_1 + 1)/2;
            A_1 = A_1(:,col_perm);
            % Sensing matrix 2
            A_2 = abs(A_1 -1); % Representing "Negative/compliment" to A_1
            
            % Effected by M
            A_1 = A_1(picks,:);
            A_2 = A_2(picks,:);
            %
            
            
            %
            sigma = 0.000007;  % noise std
            
            
            % First obseravation (one sensor)
            b_1 = A_1*Im0(:);
            %
            % add noise
            bavg = mean(abs(b_1));
            noise = randn(M,1);
            b_1_n_3e3 = b_1 + sigma*bavg*noise*100;
            
            
            % Second sensor observation
            b_2 = A_2*Im0(:);
            
            % Add noise
            noise = randn(M,1); % New noise
            bavg = mean(abs(b_2));
            b_2_n_3e3 = b_2 + sigma*bavg*noise*100;
            
            % Put signals together
            b_2 = b_1 - b_2; % Two sensors togheter simutaniously
            b_2_n_3e3 = b_1_n_3e3 - b_2_n_3e3;
            
            A_11 = A_1 - A_2; % Create combined sensing matrix
            
            
            %imwrite(uint8(Im0),['sim/', int2str(i), '_ns.PNG'], 'PNG')
            %imwrite(uint8(Im1),['sim/', int2str(i), '_bs.PNG'], 'PNG')
            % Reconstruction
            
            
            
            % A [0 1] Noiseless
            [estIm, out] = TVAL3(A_1,b_1/sqrt(N),sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro');
            imwrite(uint8(estIm),['sim/tt/', int2str(i), '_01_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_0e0_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0) ,'\nM = ', num2str(ratio),...
                '\nNoise: 0\nre_er = ', num2str(re_er) ,'\nDist: [0,1]' ])
            
            % A [0 1] n_3e3
            [estIm, out] = TVAL3(A_1,b_1_n_3e3,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro');
            imwrite(uint8(estIm),['sim/tt/', int2str(i), '_01_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
                '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [0,1]' ])
            
                        
            % A [-1 1] Noiseless
            [estIm, out] = TVAL3(A_11,b_2,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro');
            imwrite(uint8(estIm),['sim/tt/', int2str(i), '_11_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_0e0_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0) ,'\nM = ', num2str(ratio),...
                '\nNoise: 0\nre_er = ', num2str(re_er) ,'\nDist: [-1,1]' ])
            
            
            % A [-1 1] n_3e3
            [estIm, out] = TVAL3(A_11,b_2_n_3e3,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro');
            imwrite(uint8(estIm),['sim/tt/', int2str(i), '_11_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
                '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [-1,1]' ])
            
            
            
        end
    end
end


%%
% plotting
figure('Name','Barbara','Position',...
    [fullscreen(1) fullscreen(2) fullscreen(3) fullscreen(4)]);
colormap(gray);

subplot(1,3,1);
imshow(Im0,[]);
title(sprintf('Original %dx%d Barbara',sidelength,sidelength),'fontsize',16);
subplot(1,3,2);
imshow(estIm,[]);
title(sprintf('Recovered with %2.0f%% measurements',ratio*100),'fontsize',16);
xlabel(sprintf('Noise level: %2d%%  \n Rel-Err: %4.2f%%,   CPU time: %4.2fs',100*sigma,100*re_er,t),'fontsize',14);
subplot(1,3,3);
imshow(estImd,[]);
title(sprintf('Recovered with %2.0f%% measurements Double',ratio*100),'fontsize',16);
xlabel(sprintf('Noise level: %2d%%  \n Rel-Err: %4.2f%%,   CPU time: %4.2fs',100*sigma,100*re_erd,t),'fontsize',14);

%% Im2

            % Beskärd
            
            % First obseravation (one sensor)
            b_1 = A_1*Im1(:);
            %
            % add noise
            bavg = mean(abs(b_1));
            noise = randn(M,1);
            b_1_n_3e3 = b_1 + sigma*bavg*noise*100;
            
            
            % Second sensor observation
            b_2 = A_2*Im1(:);
            
            % Add noise
            noise = randn(M,1); % New noise
            bavg = mean(abs(b_2));
            b_2_n_3e3 = b_2 + sigma*bavg*noise*100;
            
            % Put signals together
            b_2 = b_1 - b_2; % Two sensors togheter simutaniously
            b_2_n_3e3 = b_1_n_3e3 - b_2_n_3e3;
            
            % A [0 1] Noiseless
            [estIm, out] = TVAL3(A_1,b_1/sqrt(N),sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im1,'fro')/norm(Im1,'fro');
            imwrite(uint8(estIm),['sim/', int2str(i), '_01_bs_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_0e0_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0) ,'\nM = ', num2str(ratio),...
                '\nNoise: 0\nre_er = ', num2str(re_er) ,'\nDist: [0,1]\nNedskalad' ])
            
            
            % A [0 1] n_3e3
            [estIm, out] = TVAL3(A_1,b_1_n_3e3,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im1,'fro')/norm(Im1,'fro');
            imwrite(uint8(estIm),['sim/', int2str(i), '_01_bs_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
                '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [0,1]\nNedskalad' ])
            
            
            
            % A [-1 1] Noiseless
            [estIm, out] = TVAL3(A_11,b_2,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im1,'fro')/norm(Im1,'fro');
            imwrite(uint8(estIm),['sim/', int2str(i), '_11_bs_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_0e0_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0) ,'\nM = ', num2str(ratio),...
                '\nNoise: 0\nre_er = ', num2str(re_er) ,'\nDist: [-1,1]\nNedskalad' ])
            
            
            % A [-1 1] n_3e3
            [estIm, out] = TVAL3(A_11,b_2_n_3e3,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im1,'fro')/norm(Im1,'fro');
            imwrite(uint8(estIm),['sim/', int2str(i), '_11_bs_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
                '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [-1,1]\nNedskalad' ])

