%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------%
%
% Machine Perception and Cognitive Robotics Laboratory
%
%     Center for Complex Systems and Brain Sciences
%
%              Florida Atlantic University
%
%------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------%
function MPCR_Sparse_Face_Camera

lambda=0.1;

L = [];

N = 0;

key=[];

name_key={};

cam = webcam(1);
preview(cam)

for loop=1:100
    
    while size(L,2)<200
        [L,N,key,name_key]=new_person(L,N,key,name_key,cam);
    end
    
    y=take_photo(cam);
    y = y - min(y(:));
    y = y / max(y(:));
    y=whiten(y);
   
    y=y(:);
 
    a=LCA(y, L, lambda);
   
        bar(a)
        pause

    b=[];
    
    for j=1:N
        
        b=[b sum(abs(a(find(key==j))))];
        
    end
    
    [b1,b2]=max(b);
    
    
    bar(b)
    drawnow()
     
end


clear cam

end



function y=take_photo(cam)

y=downscale(downscale(im2double(rgb2gray(snapshot(cam)))));

end



function [L,N,key,name_key]=new_person(L,N,key,name_key,cam)

name=input('Please Enter Name:')

name_key=[name_key name];



D=[];
Dk=[];

N=N+1;

for k = 1:100
    
    k
    p =take_photo(cam);
    
    
    
    p = p - min(p(:));
    p = p / max(p(:));
    
    
    D=[D p(:)];
    Dk=[Dk N];
    
    
end



D = bsxfun(@minus,D,mean(D)); %remove mean
fX = fft(fft(D,[],2),[],3); %fourier transform of the images
spectr = sqrt(mean(abs(fX).^2)); %Mean spectrum
D = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); %whitened L


L=[L D];
key=[key Dk];


% L = bsxfun(@minus,L,mean(L)); %remove mean
fX = fft(fft(L,[],2),[],3); %fourier transform of the images
spectr = sqrt(mean(abs(fX).^2)); %Mean spectrum
L = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); %whitened L


disp('Patient Added')

end


function y=whiten(y)

y = bsxfun(@minus,y,mean(y)); %remove mean
fX = fft(fft(y,[],2),[],3); %fourier transform of the images
spectr = sqrt(mean(abs(fX).^2)); %Mean spectrum
y = ifft(ifft(bsxfun(@times,fX,1./spectr),[],2),[],3); %whitened L

end


function y=upscale(x)

y=zeros(1,2*size(x,2));
y(1:2:end-1)=x;
y(2:2:end)=x;

end


function y=downscale(x)

y=x(1:2:end-1,1:2:end-1);

end



function [a, u] = LCA(y, D, lambda)


t=0.01;
h=0.000001;

d = h/t;
u = zeros(size(D,2),1);


for i=1:100
    
    
    a = ( u - sign(u).*(lambda) ) .* ( abs(u) > (lambda) );
    %     a=u.*(abs(u) > lambda);
    
    u =   u + d * ( D' * ( y - D*a ) - u - a  ) ;
    
    
end




end




















