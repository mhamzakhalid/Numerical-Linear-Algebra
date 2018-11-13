function [U rhs X Y] = geometryPoisson(N,f,g,task)
    h = 1/N;
    U = zeros(N+1);
    %grid
    x = 0:h:1;
    y = 0:h:1;
    [X,Y] = meshgrid(x,y);
    switch task
        case 'CG'
            %Boundary on U
            U(:,1) = g(X(:,1),Y(:,1));
            U(:,end) = g(X(:,end),Y(:,end));
            U(1,:) = g(X(1,:)',Y(1,:)');
            U(end,:) = g(X(end,:)',Y(end,:)');
            rnd = rand(N-1);
            U(2:N,2:N) = rnd(1:end,1:end);
            %Right Hand side
            rhs = f(X,Y);
        case {'MGV', 'PCG'}
            %Boundary Conditions
            U(:,1) = g(X(:,1),Y(:,1));
            rnd = rand(N-1);
            U(2:N,2:N) = rnd(1:end,1:end);
            %%Right Hand side
            rhs = zeros(N+1);
            rhs(:,:) = f(X,Y);
        otherwise 
            disp('Enter Correct Case CG or MGV')
    end
end