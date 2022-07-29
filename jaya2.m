clear all, clc
pi = 3.1415926535;
d = 42; %dimension
n = 250; %population size
for i = 1:(d/2)
    lb(i) = 0;     
    ub(i) = pi;                     %lower bound
end
for i = (d/2+1):d
    lb(i) = -80;
    ub(i) = 80;                   %upper bound   
end 
          
%lb = [-100 -100 -100]; 
%ub = [100 100 100]; 
y0 = [1648 -1.6 18148 0 940];


itermax = 30;

for i = 1:n
    for j = 1:d
        u(i,j) = lb(j) + rand.*(ub(j) - lb(j));
    end
end

    for j = 1:(d/2)
        u1(:,j) = u(:,j);
    end
    for j = (d/2+1):d
        u2(:,j-d/2) = u(:,j);
    end

tspan1 = 0:15:300;
tspan2 = 0:0.01:300;
    for j = (d/2+1):d
       hf = 1500;
       vf = 5;
    end
[T2, Y2] = MoonShipModel1(tspan2, y0, u1(1,1), u2(1,1));
%for i = 1:n
    for j = 1:(d/2)
         [T1, Y1] = MoonShipModel(tspan1, y0, u1(:,j), u2(:,j), d/2);
    J = ffunc1(T1, [Y1(j,3) Y1(j,1)], hf, vf); % + רענאפ
 % ך-וםע רענאפא
    end
%end
for iter = 1:itermax
            for j = 1:d/2
    %Finding the Xbest position
    [fmin, minind] = min(J); 
    Xbest = u2(minind,j);
     %Finding the Xworst position
    [fmax, maxind] = max(J);
    Xworst = u2(maxind,j);
    %for i = 1:n
        X = u2(:,j);
        Xnew = X + rand.*(Xbest - abs(X)) - rand.*(Xworst - abs(X));
        % Checking the bounds
        Xnew = max(Xnew, lb(j)); % Preserve lower bounds
        Xnew = min(Xnew, ub(j)); % Preserve upper bounds
        
        Fnew = ffunc1(T1, [Y1(:,3) Y1(:,1)], hf, vf);
        % Greedy selection
            if Fnew < J(j)
                u2(:,j) = Xnew;
                J(:,j) = Fnew;
            end
        end
    %end
        %Memorize the best solution
    [optval, optind] = min(J);
    BestFx(iter) = optval;
    BestX(iter,:) = u2(optind,:);
    % Iteration info
    disp(['Iteration ' num2str(iter) ...
    ': Best Cost = ' num2str(BestFx(iter))]);
    % Plot
    plot(BestFx, 'LineWidth', 0.5);
    xlabel('Iteration number');
    ylabel('Fitness value');
    %title('Convergence by iteration');
    grid on;
        
end
