function [model,dynamic,accept,sample,LnL] = HMC(model,I,epsilon)


dynamic.q              = model.B(model.sEL);
dynamic.p              = normrnd(0,1,numel(dynamic.q),1);
[dynamic.dU,model.LnL] = dLnL(model);

for i = 1:I
    [model,dynamic,accept(i)] = nUTs(model,7,epsilon);
    if isequal(accept(end),1)
        sample(i,:) = model.B(model.sEL)';
    else
        sample(i,:) = NaN(1,numel(model.sEL));
    end        
    N(i)                      = size(dynamic.q,2)
    LnL(i)                    = model.LnL;
end

function [model,dynamic,accept] = nUTs(model,K,epsilon)
dynamic.q              = model.B(model.sEL);
dynamic.p              = normrnd(0,1,numel(dynamic.q),1);
[dynamic.dU,model.LnL] = dLnL(model);
F                      = model.LnL - dynamic.p'*dynamic.p/2;
proposal               = model;
S                      = 1;
for k = 1:K
    s(k)               = sign(rand - .5);
    if ~isequal(s(end),S) && k > 1
        dynamic.p  = flip(dynamic.p,2);
        dynamic.q  = flip(dynamic.q,2);
        dynamic.dU = flip(dynamic.dU,2);
        S          = s(end);
    end
    [dynamIC,proposAL] = LeApFrOG(proposal,dynamic,2^(k - 1),s(end)*epsilon);
    stop(1)            = min(s(end)*(dynamic.q(:,end) - dynamic.q(:,1))'*dynamic.p(:,[1 end])) < 0;
    stop(2)            = (LnLikelihood(proposal) - dynamic.p(:,end)'*dynamic.p(:,end)/2) - F < -1000;
    
    if isequal(max(stop),1) || isequal(proposAL.LnL,-Inf) || isnan(proposAL.LnL) 
        break
    else
        dynamic    = dynamIC;
        proposal   = proposAL;
    end
end

alpha      = min(1,exp((proposal.LnL - dynamic.p(:,end)'*dynamic.p(:,end)/2))/exp(F));
if rand() < alpha  
    model  = proposal;
    accept = 1;
else
    accept = 0;
end

function epsilon = FindReasonableEpsilon(model,dynamic)
epsilon         = .1;
[dynamIC,modEL] = LeApFrOG(model,dynamic,1,epsilon);
alpha           = exp(modEL.LnL - model.LnL - 1/2*(dynamIC.p(:,end)'*dynamIC.p(:,end) - dynamic.p(:,end)'*dynamic.p(:,end)));
a               = 2*(alpha > 0.5) - 1;
N               = 1;

while alpha^a > 2^(-a)
    epsilon         = epsilon*(2^a);
    [dynamIC,modEL] = LeApFrOG(model,dynamic,1,epsilon);
    alpha           = exp(modEL.LnL - model.LnL - 1/2*(dynamIC.p(:,end)'*dynamIC.p(:,end) - dynamic.p(:,end)'*dynamic.p(:,end)));
    N               = N + 1;
end

function [dynamic,model] = LeApFrOG(model,dynamic,L,epsilon) % Leapfrog
sEL             = model.sEL;
for i = 1:L
    dynamic.p(:,end + 1)              = dynamic.p(:,end) + epsilon*dynamic.dU(:,end)/2;
    dynamic.q(:,end + 1)              = dynamic.q(:,end) + epsilon*dynamic.p(:,end);
    model.B(sEL)                      = dynamic.q(:,end);
    [dynamic.dU(:,end + 1),model.LnL] = dLnL(model);
    dynamic.p(:,end)                  = dynamic.p(:,end) + epsilon*dynamic.dU(:,end)/2;
end

function lnL = LnLikelihood(model) % Likelihood
if isequal(model.type,'NB1')
    alpha        = exp(model.B(end))';
    Beta         = mat2cell(model.B(1:end - 1),[size(model.events,2) size(model.X{1},2)],1);
    AeXB         = alpha.*exp(Beta{1}' + cell2mat(model.X)*kron(eye(size(model.events,2)),Beta{2})).*model.exposure;
    C            = gamma(model.events + 1./alpha)./(gamma(model.events + 1).*gamma(1./alpha));
    p            = C.*((1./(1 + AeXB)).^(1./alpha)).*((AeXB./(1 + AeXB)).^model.events);
    lnL          = sum(sum(model.w.*log(p)));
elseif isequal(model.type,'NB2')
    delta        = exp(model.B(end))';
    Beta         = mat2cell(model.B(1:end - 1),[size(model.events,2) size(model.X{1},2)],1);
    mu           = exp(Beta{1}' + cell2mat(model.X)*kron(eye(size(model.events,2)),Beta{2})).*model.exposure;
    p            = (mu.^model.events).*((1/(1 + delta)).^(mu/delta + model.events));
    lnL          = sum(sum(model.w.*log(p)));
else
    Beta         = mat2cell(model.B,[size(model.events,2) size(model.X{1},2)],1);
    eXB          = exp(Beta{1}' + cell2mat(model.X)*kron(eye(size(model.events,2)),Beta{2})).*model.exposure;
    p            = (1./factorial(model.events)).*(eXB.^model.events).*exp(-eXB);
    lnL          = sum(sum(model.w.*log(p)));
end

function [D,LnL] = dLnL(model) % gradient
LnL             = LnLikelihood(model);
d               = 10^-5;
sEL             = model.sEL;
for i = 1:numel(sEL)
    Bi           = model.B;
    Bi(sEL(i),:) = Bi(sEL(i),:) + d;
    model_i      = model;
    model_i.B    = Bi;
    D(i,1)       = (LnLikelihood(model_i) - LnL)/d;
end