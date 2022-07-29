clear all, clc
d = 2; %dimension
n = 250; %population size
lb = [0 -80]; %lower bound
ub = [pi 80]; %upper bound
itermax = 20;
dt = 15;

hf = 150; % termination height, м
Vf = 5; % termination speed, м/s

for i = 1:n
    for j = 1:d
        pos(i,j) = lb(j) + rand.*(ub(j) - lb(j));
    end
end
time = 0:dt:dt*itermax;
%Initialize gbest
iter = 1;
while iter <= itermax
    for j = 1:n
        [T, Y] = MoonShipModel(time, [1648 -1.6 18148 0 940], pos(j, 1), pos(j, 2), 0);
        fun(j, 1) = Y(itermax+1, 1);
        fun(j, 2) = Y(itermax+1, 3);
        fx(j, 1) = sqrt((fun(j, 2) - hf).^2);
        if fun(j, 2) <= 0
            fx(j, 1) = fx(j, 1) .* 10^9;
        end
    end
    [fminvalue, ind] = min(fx);
    gbest = pos(ind,:);
    Fgbest = fminvalue;
    a = 2 - 2*iter/itermax;
    i = 1;
    for i = 1:n
        X = pos(i,:);
        pos1 = pos;
        % THE best
        A1 = 2.*a.*rand(1,d) - a;
        C1 = 2.*rand(1,d);
        [alpha, alphaind] = min(fx);
        Max_Value = max(fx);
        alphapos = pos1(alphaind,:);
        Dalpha = abs(C1.*alphapos - X);
        X1 = alphapos - A1.*Dalpha;
        %pos1(alphaind,:) = [];
        %fx1 = fun(pos1);
        fx1 = fx;
        fx1(alphaind,:) = Max_Value;
        
        % 2nd best
        [beta, betaind] = min(fx1);
        betapos = pos1(betaind,:);
        A2 = 2.*a.*rand(1,d) - a;
        C2 = 2.*rand(1,d);
        Dbeta = abs(C2.*betapos - X);
        X2 = betapos - A2.*Dbeta;
        %pos1(betaind,:) = [];
        %fx1 = fun(pos1);
%         fx1(alphaind,:) = [];
        fx1(betaind,:) = Max_Value;
        
        % 3rd best
        [delta, deltaind] = min(fx1);
        deltapos = pos1(deltaind,:);
        A3 = 2.*a.*rand(1,d) - a;
        C3 = 2.*rand(1,d);
        Ddelta = abs(C3.*deltapos - X);
        X3 = deltapos - A3.*Ddelta;
        
        Xnew = (X1 + X2 + X3) ./ 3;
        %Checking the bounds
        Xnew = max(Xnew, lb); % Preserve lower bounds
        Xnew = min(Xnew, ub); % Preserve upper bounds
        %Fnew = fun(Xnew)
        %time = dt*(iter-1):dt:dt*iter+1;
        [T, Y] = MoonShipModel(time, [1648 -1.6 18148 0 940], Xnew(1, 1), Xnew(1, 2), 0);
        fun1_new(1, 1) = Y(itermax, 1);
        fun1_new(1, 2) = Y(itermax, 3);
        Fnew = sqrt((fun1_new(1, 1) - Vf).^2 + (fun1_new(1, 2) - hf).^2);
        if fun1_new(1, 2) <= 0
            Fnew = Fnew .* 10^9;
        end
        
        % Greedy selection
        if Fnew < fx(i)
            pos(i,:) = Xnew;
            fx(i,:) = Fnew;
        end
    end
    % Udpate gbest
    [fmin,find] = min(fx);
    if fmin < Fgbest
        Fgbest = fmin;
        gbest = pos(find,:);
    end
    u1_best(1, iter) = gbest(1, 1);
    u2_best(1, iter) = gbest(1, 2);
    %Memorize the best solution
    [optval, optind] = min(fx);
    BestFx(iter) = optval;
    BestX(iter,:) = pos(optind,:);
    % Iteration info
    disp(['Iteration ' num2str(iter) ...
    ': Best Cost = ' num2str(BestFx(iter))]);
    iter = iter + 1;
end

% Plot
figure(1); 
plot(BestFx, 'LineWidth', 1.5);
xlabel('Номер итерации');
ylabel('Значение функционала J');
title('Минимизация функционала');
grid on;

%figure(2);
%plot(T, Y(:, 1));
%grid on
%legend('V(t)')
%xlabel('Время, сек');
%ylabel('Модуль скорости движения космического аппарата, м/с');

%figure(3);
%plot(T, Y(:, 2));
%grid on
%legend('\theta(t)')
%xlabel('Время, сек');
%ylabel('Угол наклона траектории относительно гравитационной вертикали, рад');

%figure(4);
%plot(T, Y(:, 3));
%grid on
%legend('h(t)')
%xlabel('Время, сек');
%ylabel('Высота полёта космического аппарата, м');

%figure(5);
%plot(T, Y(:, 4));
%grid on
%legend('\phi(t)')
%xlabel('Время, сек');
%ylabel('Угол дальности вдоль поверхности, рад');
%grid on

%figure(6);
%plot(T, Y(:, 5));
%grid on
%legend('m(t)')
%xlabel('Время, сек');
%ylabel('Масса космического аппарата, кг');
