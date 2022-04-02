function [q] = grid_hexagonal(N_molecole,l)
    a=l*cos(pi/3);
    b=l*sin(pi/3);
    N= ceil(sqrt(N_molecole));
    AA=length(l:2*a:2*a*N);
    
    q(1,1)=0;
    q(2,1)=0;
    
    i=1;
    first=0;
    step=0;
    while i<N+1
        q=[q(1,:),l:2*a:2*a*N;
           q(2,:), step*ones(1, AA)];
        if first==0
            i=i+1;
            first=1;
            step=step+l;
        else
            i=i+3;
            first=0;
            step=step+2*b+l;
        end
    end

    q(:,1)=[];
    i=3;
    first=0;
    step=l+b;
    while i<N+1
        q=[q(1,:),l+a:2*a:2*a*N+a;
           q(2,:), step*ones(1, AA)];
        if first==0
            i=i+1;
            first=1;
            step=step+l;
        else
            i=i+3;
            first=0;
            step=step+2*b+l;
        end
    end

end