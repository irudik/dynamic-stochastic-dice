function [q,v,c] = function_iterate(outboundsfile,s,v,q,...
    VarsExog,LogicalOpts,Params,space,space_terminal,feedback_nodes,feedback_weights,period)

%   function_iterate: Loops over the collocation grid points and maximizes
%   the right hand side of the Bellman at some time 't'. Checks whether the
%   solution results in transitioning outside the domain next period, and
%   reports the number of quadrature points evaluated outside the domain.

if LogicalOpts.outbounds == 1
    outbound = [];
    outboundstring = {};    %outbounds_nonzero = 0;
end

% Store temporary structure of exogenous trajectories at this year
value_variables = struct(...
    'betaeffi',VarsExog.betaeff,'Bi',VarsExog.Bs(period+1),...
    'Psii',VarsExog.Psi(period+1),'eminti',VarsExog.emint(period+1),...
    'EFi',VarsExog.EF(period+1),'EFi_next',VarsExog.EF_next(period+1),...
    'Ati',VarsExog.At(period+1),'Lti',VarsExog.Lt(period+1),...
    'Atni',VarsExog.Atn(period+1),'Ltni',VarsExog.Ltn(period+1),'sigmai',VarsExog.sigma(period+1),...
    'eff_adjusti',VarsExog.eff_adjust(period+1),'eff_adjustni',VarsExog.eff_adjustn(period+1));

loopcount = tic;

parfor i = 1:space.num_nodes
    % Pull this collocation node's quadrature nodes
    feedback = feedback_nodes(:,i);
    % Get the 1 dimensional weights to do a 2 dimensional expectation over
    % the feedback and noise (if we want to allow for a truncated feedback
    % distribution)
    weights = feedback_weights(:,i)*Params.tw';
    weights = reshape(weights,numel(weights),1);
    
    % Maximize the Bellman at node 'i'
    [qnew vnew] = maxbell(s(i,:),q(i,:),...
        value_variables,Params,LogicalOpts,space,space_terminal,feedback,weights);
    q(i,:) = qnew;
    v(i,1) = -vnew;
end


toc(loopcount)

% Check if we transition outside the grid
if LogicalOpts.outbounds
    run check_domain
end

% Update coefficients
c = smol_fit(space,v);