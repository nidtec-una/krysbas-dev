function miter = pdrule(m,minitial,mmin,res,iter,mstep,mmax, alpha, delta) %implementacion Cabral

 ap=0;
 ad=0;

if iter >3
    
mj= m + ceil(alpha*(res(iter,:)/res(iter-1,:)) + delta*((res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))));
%mj= m + round(alpha*(res(iter,:)/res(iter-1,:)) + delta*((res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))));
    
 ap=res(iter,:)/res(iter-1,:);
 ad=(res(iter,:) - res(iter-2,:))/(2*res(iter-1,:));

elseif iter >2
    mj= m + ceil(alpha*(res(iter,:)/res(iter-1,:)));
    %    mj= m + round(alpha*(res(iter,:)/res(iter-1,:)));

else
    mj=minitial;
end

% ap=res(iter,:)/res(iter-1,:)
% if iter >3
% ad=(res(iter,:) - res(iter-2,:))/(2*res(iter-1,:))
%end


% if mj <= 0
%      mj=abs(mj) + mstep;
% end
if mj < mmin
    minitial = minitial + mstep;
    mj=minitial;
end

if mj > mmax
%    miter=minitial;
%     miter=miter-mstep;
    %mj=mmax-mstep;
    mj=mmax;
end

miter= [mj minitial] ;
miter= [mj minitial ap ad] ;

    

    

