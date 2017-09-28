
% Extraherar signalen, kompenserar f�r variationer i ljus
% rekonstruerar med TVAL3 som �r explicit gjord f�r WH och compressed 
% imaging (men kan anv�ndas f�r tex slump�ssiga m�tmatriser ).
% K�r lie enkel bildbehandling sist

% ----------
% TVAL3 users guide: http://www.caam.rice.edu/~optimization/L1/TVAL3/v.beta/User_Guide_beta2.4.pdf
% Viktigaste �r: 
% opts.mu (ju h�gre desto er detaljrik och (mer brus)) 
% och 
% opts.beta (ungef�r som opts.mu)

% opts.tol (L�gre v�rde --> L�ngsammare, borde ge b�ttre bild) (H�gre v�rde --> snabbare)

% F�rslag fr�n Andreas:
% opts.mu - minst 2^10 (�ka med potenser)
% opts.beta - unders�k
% opts.tol - H�ll l�gt f�r b�ttre bildkvalitet

% opts.nonneg = true - sant f�r bilder
% opts.isreal = true - s�ker en bild som inte �r komplex
% ---------

% F�r att ta bort tex r�relse fr�n ett objekt under m�tsignalen:
% Ta bort samples i y_dn_f[a:b]
% och motsvarande samples i picks[a:b]. picks = row_perm[1:M]

clear all
clear opts
% set the optional paramaters of TVAL3
opts.mu = 2^13; opts.beta = 2^5; % 2^5 , 2^7
opts.maxcnt = 10; opts.tol_inn = 1e-5; %5 300
opts.tol = 1E-10; opts.maxit = 1000; % WH = 4500
opts.mu0 = 2^10; opts.beta0 = 2^5;    % trigger continuation scheme
opts.nonneg = true; opts.isreal = true
%
for nr = [24]%[30:38 42:46 60:65 70:84 92:93 95:96]
    nr
    sidelength = 512;
    N = sidelength^2;
    ratio = 0.3
    M = round(N*ratio/48)*48 % Antal m�tmatriser. Must match number of Measurement matrix in sampling
    f_s = 69120; % Sampling freq in A/D (Andreas valde 69120 som ger 48 samples per m�tmatris vid 60Hz och 24 m�tmatris per frame)
    y_dn_f = signal_process_48([int2str(nr),'_', int2str(sidelength),'.tdms'], sidelength, ratio, f_s, false);
    
%     M = round(N*ratio/3)*3
%     f_s = 18000;
%     y_dn_f = signal_process([int2str(nr),'_', int2str(sidelength),'.tdms'], sidelength, ratio, f_s, true);
%     
    for ratio_i =  1
        ratio_i
        % kompensera f�r ljus variationer moving mean
        y_dn = y_dn_f(1:end*ratio_i);
        x = 1:length(y_dn);
        dt_mm = movmean(y_dn,75);
        y_final = y_dn - dt_mm;
        y_final = y_final - min(y_final);
        %
        
        
        figure(6)
        plot(x,y_dn,'r',x,y_final-0.007,'b',x,dt_mm,'g')
        pause(0.1)
        
        % Tar bort r�relser (m�tv�rden som ligger v�ldigt h�gt/l�gt)
        %--------------------------------------------
        figure(10)
        m = x./x.*mean(y_final);
        
        s = std(y_final)*3.8;
        s1 = m + s;
        s2 = m - s;
        plot(x,y_final,'b',x,m,'g',x,s1,'r',x,s2,'r')
        pause(0.1)
        %
        b = [];
        for i = 1:length(y_final)
            if (y_final(i) > s1) | (y_final(i) < s2)
                b = [b min(length(y_final),max(1,i-50:i+50))];
            end
        end
        c = unique(b);
        d = setxor(x,c);
        y_ff = y_final(d);
        figure(11); plot(y_ff); pause(0.1)
        
        load(['row_perm_', int2str(sidelength)]);
        load(['col_perm_' int2str(sidelength)]);
        %picks = row_perm(1:M);
        picks = row_perm(d);
        % -------------------------------------------
        
        % Initilize Sensing matrix
        A = @(x,mode) dfA(x,picks,col_perm,mode);
        %
        
        [estIm, out] = TVAL3(A,y_ff,sidelength,sidelength,opts);
        estIm = estIm - min(estIm(:));
        estImm = medfilt2(estIm,[5 5]);
        estImm = flipdim(estImm ,1);           %# vertical flip
        estImNorm = uint16(imadjust(estImm./max(estImm(:)))*2^16);
        figure(15); imshow(estImNorm, [0 2^16]);
        
        imwrite(estImNorm,['./im/', int2str(nr),'_', int2str(sidelength),'_13_2_m',int2str(ratio*ratio_i*100),'.PNG'],'BitDepth',16)
    end
end
%%
figure(13); imshow(estImm, [0 0.2]);