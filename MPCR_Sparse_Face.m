function MPCR_Sparse_Face
clear all
close all
clc

N=20; %number of patients
M=28; %Number of photos per patient

lambda=0.1;

[L,patient_names,patient_names_key]=load_library(N,M);

r=randperm(size(L,2));

L=L(:,r);
patient_names_key=patient_names_key(r);

% for i=1:size(L,2)
%
%     imagesc(reshape(L(:,i),120,128))
%     colormap(gray)
%     patient_names_key(i)
%     patient_names{patient_names_key(i)}
%     pause
%
% end

L1=L;

p1=[];

t=[];

for k=1:size(L,2)%randperm(size(L,2))%
    
    L=L1;
    
    y=L(:,k);
    L(:,k)=0;
    
    
    a=LCA(y, L, lambda);
    
    
    b=[];
    
    for j=1:N
        
        b=[b sum(abs(a(find(patient_names_key==j))))];
        
    end
    
    [b1,b2]=max(b);
    
    
    figure(1)
    subplot(511)
    bar(b)
    
    p=[b2 patient_names_key(k)]
    p1=[p1; p];
    
    subplot(512)
    plot(p1)
     
    %['Class:' patient_names{patient_names_key(k)} '  Test:' patient_names{patient_names_key(b2)}]

    subplot(513)
    imagesc(reshape(y,120,128))
    colormap(gray)

    subplot(514)
    bar(abs(a))
    
    subplot(515)
    t=[t b2==patient_names_key(k)];
    hist(t,2);
    sum(t)/length(t)
    
    drawnow()
     
end

end



function [L,patient_names,patient_names_key]=load_library(N,M)

cd faceinput

patient_names={'an2i','at33','boland','bpm','ch4f','cheyer','choon','danieln','glickman','karyadi','kawamura','kk49','megak','mitchell','night','phoebe','saavik','steffi','sz24','tammo'};
patient_names_key=[];%ceil((1:(M*N))/M);

L=[]; %Library

for i =1:size(patient_names,2)
    
    patient_names{i};
    dr1=dir([patient_names{i} '*open.pgm']);
    dr2=dir([patient_names{i} '*sunglasses.pgm']);
    
    f1={dr1.name, dr2.name}; % get only filenames to cell
    
    D=[]; %Dictionary
    
    for j=1:M%length(f1) % for each image
        
        a1=f1{j};
        
        b1=im2double(imread(a1));
        
        b1=b1(1:end)';
        
        b1 = b1 - min(b1(:));
        b1 = b1 / max(b1(:));
        

        D=[D b1];
        
        patient_names_key=[patient_names_key i];
        
    end
    
    D = bsxfun(@minus,D,mean(D)); %remove mean
    fX = fft(fft(D,[],2),[],3); %fourier transform of the images
    spectr = sqrt(mean(abs(fX).^2)); %Mean spectrum
    D = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); %whitened L
    
    L=[L D];
    
    
end


% L = bsxfun(@minus,L,mean(L)); %remove mean
fX = fft(fft(L,[],2),[],3); %fourier transform of the images
spectr = sqrt(mean(abs(fX).^2)); %Mean spectrum
L = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); %whitened L


end



function [a, u] = LCA(y, D, lambda)


t=.01;
h=.00001;

d = h/t;
u = zeros(size(D,2),1);


for i=1:100
    
    a=u.*(abs(u) > lambda);
    
    u =   u + d * ( D' * ( y - D*a ) - u - a  ) ;
    
    
end




end



