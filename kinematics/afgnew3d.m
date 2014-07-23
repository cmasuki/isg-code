function[odata]=afgnew3d(idata,v);
% function AFGNEW; 
% calculates numerical derivative for position- and velocity data
% constants chosen for best frequency response around 50-100Hz sample frequency
% v = sample frequency (Hz)
% c (in program) = differentiation constants
% DirkJan, april 1991
% rewritten version of afgnew, can be used for both two-dimensional and
% three-dimensional arrays. In the latter case, time should be the last
% dimension (i.e 'k' in [i,j,k])

if ndims(idata)==3,
    [I,J,K]=size(idata);

    c=[-1.5 -4 4 1.5]/14*v;
    c=diag(c);

    for i=1:I,
            for j=1:J
        tussen(:,1)=[0;0;squeeze(idata(i,j,1:K-2))];
        tussen(:,2)=[0;squeeze(idata(i,j,1:K-1))];
        tussen(:,3)=[squeeze(idata(i,j,2:K));0];
        tussen(:,4)=[squeeze(idata(i,j,3:K));0;0];

        tussen=tussen*c;
        odata(i,j,:)=sum(tussen');
    end
    end


    odata(:,:,1:2)=zeros(I,J,2);
    odata(:,:,K-1:K)=zeros(I,J,2);
end

if ndims(idata)==2,
    
    [m,n]=size(idata);
    
    gekanteld=0;
    
    if m==3, 
        gekanteld=1;
        idata=idata'; % uses lying vectors.....
        [m,n]=size(idata);
    end

    c=[-1.5 -4 4 1.5]/14*v;
    c=diag(c);

    for i=1:n,
        tussen(:,1)=[0;0;idata(1:m-2,i)];
        tussen(:,2)=[0;idata(1:m-1,i)];
        tussen(:,3)=[idata(2:m,i);0];
        tussen(:,4)=[idata(3:m,i);0;0];

        tussen=tussen*c;
        odata(:,i)=sum(tussen')';
    end

    odata(1:2,:)=zeros(2,n);
    odata(m-1:m,:)=zeros(2,n);
    
    if gekanteld==1, 
        odata=odata'; %back to standing vectors
    end

end
