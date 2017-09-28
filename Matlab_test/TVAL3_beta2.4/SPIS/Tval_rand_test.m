clear all; %close all;
% problem size
sidelength = 256;
N = sidelength^2;
fullscreen = get(0,'ScreenSize')*0.7;

% set the optional paramaters
clear opts
opts.mu = 2^5; opts.beta = 2^7;
opts.maxcnt = 50; opts.tol_inn = 1e-4;
opts.tol = 1E-6; opts.maxit = 1324;
opts.mu0 = 2^4; opts.beta0 = 2^0;    % trigger continuation scheme
opts.nonneg = true; opts.isreal = true


% 
% Create Sensing matrix
%A = randi(2,N,N,'single')-1;
ratio = 0.2;
M = round(N*0.2);
A_1 = randi(2,M,N,'single')-1;

for mu = [2^9]
    opts.mu = mu;
    for i = [90 91 93]
        % Importera SWIR-bild
        Im0 = imread(['../im/Shapes/',int2str(i),'.PNG']);
        Im0 = single(Im0(:,:,1));
        %[ Im0, ~ ] = get_swir_im( i, 256 );
        %Im0 = single(Im0(:,:,1));
            
            % Effected by M
            %A_1 = single(A(1:M,:));
            %
            % First obseravation (one sensor)
            b_1 = A_1*Im0(:)./sqrt(N);
            %
            % add noise
            sigma = 0.001;  % noise std
            bavg = mean(abs(b_1));
            noise = randn(M,1);
            b_1_n_3e3 = b_1 + sigma*bavg*noise;
            
 
            
            % Reconstruction
            % A [0 1] Noiseless
            [estIm, out] = TVAL3(A_1,b_1,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro');
            figure(1); imshow(estIm, [0 255])
            pause(0.1);
            imwrite(uint8(estIm),['sim/256/', int2str(i), '_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_0e0_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0) ,'\nM = ', num2str(ratio),...
                '\nNoise: 0\nre_er = ', num2str(re_er) ,'\nDist: [0,1]\nNedskalad' ])
            
            % A [0 1] n_3e3
            [estIm, out] = TVAL3(A_1,b_1_n_3e3,sidelength,sidelength,opts);
            estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
            re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro');
            figure(2); imshow(estIm, [0 255])
            pause(0.1);
            imwrite(uint8(estIm),['sim/256/', int2str(i), '_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
                'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
                '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [0,1]\nNedskalad' ])
            
            
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

