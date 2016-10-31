%   check_domain: Checks whether the solution to the time 't' Bellman
%   results in transitioning outside the collocation domain in time 't+1'
%   and reports the number of quadrature points that result in evaluations
%   of the continuation value outside the collocation domain.

k_jumps_up = 0;
k_jumps_down = 0;
M_jumps_up = 0;
M_jumps_down = 0;
Mu_jumps_up = 0;
Mu_jumps_down = 0;
Ml_jumps_up = 0;
Ml_jumps_down = 0;
T_jumps_up = 0;
T_jumps_down = 0;
Tocean_jumps_up = 0;
Tocean_jumps_down = 0;
mean_jumps_up = 0;
mean_jumps_down = 0;
var_jumps_up = 0;
var_jumps_down = 0;
g = zeros(size(s,1),size(s,2));


value_variables = struct(...
    'betaeffi',VarsExog.betaeff,'Bi',VarsExog.Bs(space.year/Params.timestep+1),...
    'Psii',VarsExog.Psi(space.year/Params.timestep+1),'eminti',VarsExog.emint(space.year/Params.timestep+1),...
    'EFi',VarsExog.EF(space.year/Params.timestep+1),'EFi_next',VarsExog.EF_next(space.year/Params.timestep+1),...
    'Ati',VarsExog.At(space.year/Params.timestep+1),'Lti',VarsExog.Lt(space.year/Params.timestep+1),...
    'Atni',VarsExog.Atn(space.year/Params.timestep+1),'Ltni',VarsExog.Ltn(space.year/Params.timestep+1),...
    'sigmai',VarsExog.sigma(space.year/Params.timestep+1),...
    'eff_adjusti',VarsExog.eff_adjust(space.year/Params.timestep+1),'eff_adjustni',VarsExog.eff_adjustn(space.year/Params.timestep+1));

for ii = 1:size(q,1)
    
    % Conditional on selected controls from the Bellman maximization,
    % determine what our next state will be at both the simulated value and
    % at all of our quadrature points so that we get a sense of where the
    % value function is being evaluated
    [g,maxo(ii,:),mino(ii,:)] = transitions(s(ii,:),NaN,q(ii,1),q(ii,2),...
        value_variables,Params,LogicalOpts,space,feedback_nodes(:,ii));
    
    if g(1,4)<Params.smin(4) %check M
        M_jumps_down = M_jumps_down +1;
        outbound(end+1,:) = [4 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'M' 'down' 'it:' '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,4)>Params.smax(4)
        M_jumps_up = M_jumps_up +1;
        outbound(end+1,:) = [4 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'M' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    if g(1,3)<Params.smin(3) %check M
        Mu_jumps_down = Mu_jumps_down +1;
        outbound(end+1,:) = [3 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'Mu' 'down' 'it:' '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,3)>Params.smax(3)
        Mu_jumps_up = Mu_jumps_up +1;
        outbound(end+1,:) = [3 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'Mu' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    if g(1,2)<Params.smin(2) %check M
        Ml_jumps_down = Ml_jumps_down +1;
        outbound(end+1,:) = [2 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'Ml' 'down' 'it:' '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,2)>Params.smax(2)
        Ml_jumps_up = Ml_jumps_up +1;
        outbound(end+1,:) = [2 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'Ml' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    if g(1,1)<Params.smin(1) %check k
        k_jumps_down = k_jumps_down+1;
        outbound(end+1,:) = [1 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'k' 'down' 'it:' '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,1)>Params.smax(1)
        k_jumps_up = k_jumps_up +1;
        outbound(end+1,:) = [1 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'k' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    
    if g(1,5)<Params.smin(5) %check T actual cs
        T_jumps_down = T_jumps_down+1;
        outbound(end+1,:) = [5 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'T actual' 'down' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    if g(1,5)>Params.smax(5) % check T actual cs
        T_jumps_up = T_jumps_up +1;
        outbound(end+1,:) = [5 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'T actual' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:6),5) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    % Ignore if within machine precision: rounding error from temperature
    % transitions
    if g(1,6)<Params.smin(6) && abs(g(1,7)-Params.smin(6)) > eps %check To
        Tocean_jumps_down = Tocean_jumps_down+1;
        outbound(end+1,:) = [6 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'To' 'down' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:6),5) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,6)>Params.smax(6)
        Tocean_jumps_up = Tocean_jumps_up +1;
        outbound(end+1,:) = [6 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'To' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:6),5) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    if g(1,7)<Params.smin(7) %check mean
        mean_jumps_down = mean_jumps_down+1;
        outbound(end+1,:) = [7 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'mean' 'down' 'it:' '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,7)>Params.smax(7)
        mean_jumps_up = mean_jumps_up +1;
        outbound(end+1,:) = [7 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'mean' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
    
    if g(1,8)<Params.smin(8) %check variance
        var_jumps_down = var_jumps_down+1;
        outbound(end+1,:) = [8 -1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'k' 'down' 'it:' '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    elseif g(1,8)>Params.smax(8)
        var_jumps_up = var_jumps_up +1;
        outbound(end+1,:) = [8 1 0 s(ii,:) g(1,:)];
        outboundstring(end+1,:)= {'k' 'up' 'it:'  '|' 'node' num2str(ii) 'unc:' 'none' '|' 'jump: [' mat2str(s(ii,1:8),4) '] to ['  mat2str(g(1,1:8),4) ']'};
    end
end

% Write out-of-domain transitions to a text file
it=1;
if size(outbound,1)>0
    logical_wroteoutbounds = 1;
    dlmwrite([outboundsfile,num2str(it),'.txt'],outbound,'newline','pc');
    fid = fopen([outboundsfile,num2str(period),'.txt'], 'wt');
    [stringrows, stringcols] = size(outboundstring);
    for ii=1:stringrows
        fprintf(fid, '%s\t ',outboundstring{ii,1:end-1}); %comma alternative: '%s, ' instead
        fprintf(fid,'%s\n',outboundstring{ii,end});
    end
    fclose(fid);
end

% Report maximum and minimum temperature and mean belief evaluations at our
% quadrature points
disp(['Max temp eval: ' num2str(max(maxo(:,1))) ' min temp eval: ' num2str(min(mino(:,1)))]);
disp(['Max exp eval: ' num2str(max(maxo(:,2))) ' min exp eval: ' num2str(min(mino(:,2)))]);


% Report how many collocation nodes had a quadrature points evaluated outside
% the collocation domain
disp(['TEMP: ' num2str(sum(maxo(:,1) > Params.smax(5))) ' quad nodes jumped above ' num2str(Params.smax(5)) ' '...
    num2str(sum(mino(:,1) < Params.smin(5))) ' quad nodes jumped below ' ...
    num2str(Params.smin(5)) '.']);
disp(['EXP: ' num2str(sum(maxo(:,2) > Params.smax(7))) ' quad nodes jumped above ' num2str(Params.smax(7)) ' '...
    num2str(sum(mino(:,2) < Params.smin(7))) ' quad nodes jumped below ' ...
    num2str(Params.smin(7)) '.']);

% Report how many collocation points actually transitioned outside the
% collocation domain
disp(['Total jumps out of the state space this period: k [low,high] =[',num2str(k_jumps_down),' ',num2str(k_jumps_up),...
    '], M =[',num2str(M_jumps_down),' ',num2str(M_jumps_up),'], Mu =[',num2str(Mu_jumps_down),...
    ' ',num2str(Mu_jumps_up),'], Ml =[',num2str(Ml_jumps_down),...
    ' ',num2str(Ml_jumps_up),'], T =[',num2str(T_jumps_down),' ',num2str(T_jumps_up),...
    '], To  =[',num2str(Tocean_jumps_down),' ',num2str(Tocean_jumps_up), '], mean =[',num2str(mean_jumps_down),' ',num2str(mean_jumps_up),...
    '], var  =[',num2str(var_jumps_down),' ',num2str(var_jumps_up), '],  (at year ',num2str(period*Params.timestep),').']);