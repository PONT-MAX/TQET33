%clear all; %close all;

% set the optional paramaters
clear opts
opts.mu = 2^9; opts.beta = 2^7;
opts.maxcnt = 300; opts.tol_inn = 1e-5; %5
opts.tol = 1E-10; opts.maxit = 2500; % 10
opts.mu0 = 2^4; opts.beta0 = 2^0;    % trigger continuation scheme
opts.nonneg = true; opts.isreal = true

% Im size
for sidelength = [256 512]
    
    N = sidelength^2;
    
    % permutation to WH-matrix
    load(['row_perm_', int2str(sidelength)]);
    load(['col_perm_' int2str(sidelength)]);
    
    for ratio = [0.05 0.1 0.15 0.2 0.25 0.3]
        % Ratio
        M = round(N*ratio);
        picks = row_perm(1:M);
        
        % Initilize Sensing matrix
        A = @(x,mode) dfA(x,picks,col_perm,mode);
        %A_1 = Phi(N,ratio);
        %A_1 = randi(2,M,N) -1 ;
        
        for i = 1:21
            % Importera SWIR-bild
            [ Im0, ~ ] = get_swir_im( i, sidelength );
            Im0 = double(Im0(:));
            
            % First obseravation (one sensor)
            %b_1 = A_1*Im0(:);%./sqrt(N);
            
            
            b_1 = Phi_512_s(ratio,Im0);
            b_1 = noise_b(b_1,5);
            % Second sensor observation
            %b_2 = abs(A_1-1)*Im0(:);%./sqrt(N);
            %b_2 = noise_b(b_2);
            
            % Put 1 -1 signal together
            %b_2 = b_1 - b_2;
            %A_2 = A_1;
            %A_2(A_2 < 1) = -1;
            
            % A [0 1] Noiseless
            for nl = 1:6
                % [0 1]
                [estIm, out] = TVAL3(A,b_1(:,nl),sidelength,sidelength,opts);
                estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
                figure(i); imshow(estIm, [0 255]);
                %imwrite(uint8(estIm),['res_WH/',int2str(sidelength),'/', int2str(i), '_0_',...
                %    num2str(ratio*100),'_' int2str(nl), '.PNG']);
                
            end
        end
    end
end

% Im size
for sidelength = [128 256]
    
    N = sidelength^2;
    
    % permutation to WH-matrix
    load(['row_perm_', int2str(sidelength)]);
    load(['col_perm_' int2str(sidelength)]);
    
    for ratio = [0.1 0.2 0.3]
        % Ratio
        M = round(N*ratio);
        picks = row_perm(1:M);
        
        % Initilize Sensing matrix
        A = @(x,mode) dfA(x,picks,col_perm,mode);
        A_1 = Phi(N,ratio);
        %A_1 = randi(2,M,N) -1 ;
        
        for i = 0:21
            % Importera SWIR-bild
            [ Im0, ~ ] = get_swir_im( i, sidelength );
            Im0 = double(Im0(:));
            
            %imwrite(uint8(Im0),['sim/', int2str(i), '_ns.PNG'], 'PNG')
            
            % First obseravation (one sensor)
            b_1 = A_1*Im0(:);%./sqrt(N);
            
            
            %b_1 = Phi_512_s(ratio,Im0)
            b_1 = noise_b(b_1,3);
            % Second sensor observation
            %b_2 = abs(A_1-1)*Im0(:);%./sqrt(N);
            %b_2 = noise_b(b_2);
            
            % Put 1 -1 signal together
            %b_2 = b_1 - b_2;
            %A_2 = A_1;
            %A_2(A_2 < 1) = -1;
            
            % A [0 1] Noiseless
            for nl = 1:4
                % [0 1]
                [estIm, out] = TVAL3(A,b_1(:,nl),sidelength,sidelength,opts);
                estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
                %figure(i); imshow(estIm, [0 255]);
                imwrite(uint8(estIm),['res_WH/',int2str(sidelength),'/', int2str(i), '_0_',...
                    num2str(ratio*100),'_' int2str(nl), '.PNG']);
                
            end
        end
    end
end

opts.maxcnt = 300; opts.tol_inn = 1e-3; %5
opts.tol = 1E-6; opts.maxit = 2500; % 10
% Im size
for sidelength = [128 256]
    
    N = sidelength^2;
    
    % permutation to WH-matrix
    %load(['row_perm_', int2str(sidelength)]);
    %load(['col_perm_' int2str(sidelength)]);
    
    for ratio = [0.1 0.2 0.3]
        % Ratio
        M = round(N*ratio);
        %picks = row_perm(1:M);
        
        % Initilize Sensing matrix
        %A = @(x,mode) dfA(x,picks,col_perm,mode);
        %A_1 = Phi(N,ratio);
        A_1 = randi(2,M,N) -1 ;
        
        for i = 0:21
            % Importera SWIR-bild
            [ Im0, ~ ] = get_swir_im( i, sidelength );
            Im0 = double(Im0(:));
            
            %imwrite(uint8(Im0),['sim/', int2str(i), '_ns.PNG'], 'PNG')
            
            % First obseravation (one sensor)
            b_1 = A_1*Im0(:);%./sqrt(N);
            
            
            %b_1 = Phi_512_s(ratio,Im0)
            b_1 = noise_b(b_1,3);
            % Second sensor observation
            b_2 = abs(A_1-1)*Im0(:);%./sqrt(N);
            b_2 = noise_b(b_2,3);
            
            % Put 1 -1 signal together
            b_2 = b_1 - b_2;
            A_2 = A_1;
            A_2(A_2 < 1) = -1;
            
            % A [0 1] Noiseless
            for nl = 1:4
                % [0 1]
                [estIm, out] = TVAL3(A_1,b_1(:,nl),sidelength,sidelength,opts);
                estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
                %figure(i); imshow(estIm, [0 255]);
                imwrite(uint8(estIm),['res_WH/r/',int2str(sidelength),'/', int2str(i), '_0_',...
                    num2str(ratio*100),'_' int2str(nl), '.PNG']);
                
                [estIm, out] = TVAL3(A_2,b_2(:,nl),sidelength,sidelength,opts);
                estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
                %figure(i); imshow(estIm, [0 255]);
                imwrite(uint8(estIm),['res_WH/r/',int2str(sidelength),'/', int2str(i), '_1_',...
                    num2str(ratio*100),'_' int2str(nl), '.PNG']);
                
            end
        end
    end
end


%%

% M = round(N*0.2);
%A_1 = randi(2,M,N,'single') - 1;
%A_2 = abs(A_1 - 1);
%A_11 = A_1 - A_2; % Create combined sensing matrix
for i = [1]
    % Importera SWIR-bild
    [ Im0, ~ ] = get_swir_im( i, sidelength );
    Im0 = double(Im0);
    %Im0 = imread(['../im/Shapes/',int2str(i),'.PNG']);
    %Im0 = single(Im0(:,:,1));
    for ratio = [0.1]
        M = round(N*ratio);
        % Create Sensing matrix
        
        % Sensing matrix 1
        
        % Sensing matrix 2
        %A_2 = abs(A_1 -1); % Representing "Negative/compliment" to A_1
        
        % Effected by M
        %A_1 = A_1(1:M,:);
        %A_2 = A_2(picks,:);
        %
        %
        %sigma = 0.0001;  % noise std
        
        
        % First obseravation (one sensor)
        b_1 = A_1*Im0(:)./sqrt(N);
        %
        % add noise
        %bavg = mean(abs(b_1));
        %noise = randn(M,1);
        %b_1_n_3e3 = b_1 + sigma*bavg*noise;
        
        % Second sensor observation
        b_2 = abs(A_1-1)*Im0(:)./sqrt(N);
        
        % Add noise
        %noise = randn(M,1); % New noise
        %bavg = mean(abs(b_2));
        %b_2_n_3e3 = b_2 + sigma*bavg*noise;
        
        % Put signals together
        b_2 = b_1 - b_2; % Two sensors togheter simutaniously
        %b_2_n_3e3 = b_1_n_3e3 - b_2_n_3e3;
        
        %A_11 = A_1 - A_2; % Create combined sensing matrix
        
        
        %imwrite(uint8(Im0),['sim/', int2str(i), '_ns.PNG'], 'PNG')
        %imwrite(uint8(Im1),['sim/', int2str(i), '_bs.PNG'], 'PNG')
        % Reconstruction
        
        
        
        % A [0 1] Noiseless
        [estIm, out] = TVAL3(A,b_1,sidelength,sidelength,opts);
        estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
        re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro')*100
        figure(1); imshow(estIm, [0 255])
        %imwrite(uint8(estIm),['sim/', int2str(i), '_01_ns_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_0e0_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
        %    'Description', ['Image: ', int2str(0) ,'\nM = ', num2str(ratio),...
        %    '\nNoise: 0\nre_er = ', num2str(re_er) ,'\nDist: [0,1]\nNedskalad' ])
        
        % A [0 1] n_3e3
        %[estIm, out] = TVAL3(A_2,b_2,sidelength,sidelength,opts);
        %estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
        %re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro')*100
        %figure(2); imshow(estIm, [0 255])
        %imwrite(uint8(estIm),['sim/', int2str(i), '_01_ns_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
        %    'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
        %    '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [0,1]\nNedskalad' ])
        
        
        % A [-1 1] n_3e3
        [estIm, out] = TVAL3(A,b_2,sidelength,sidelength,opts);
        estIm = estIm - min(estIm(:)); estIm = imadjust(estIm./max(estIm(:)))*255;
        re_er = norm(estIm-Im0,'fro')/norm(Im0,'fro')*100
        figure(4); imshow(estIm, [0 255])
        %imwrite(uint8(estIm),['sim/', int2str(i), '_11_ns_r_', num2str(ratio*100), '_mu_', int2str(opts.mu), '_3e3_reer_', num2str(round(re_er*100)), '.PNG'], 'PNG','Author', 'A. Brorsson',...
        %    'Description', ['Image: ', int2str(0),  '\nM = ', num2str(ratio),...
        %    '\nNoise: 3e-3\nre_er = ', num2str(re_er) ,'\nDist: [-1,1]\nNedskalad' ])
        
        
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

