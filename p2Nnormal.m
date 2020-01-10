% This script plots function f(p,N) = c/p*(1-(1-p)^N) for the paper on
% atatcking UTS-CAPTCHA.
% The results revealed that using this function does not help in practice
% due to two facts: 1) when N is small the function has a more converged
% form (very flat) which is very sensitive to noise; 2) when N is large the
% function is very sensitive to small change of p and any numeric algorithm
% will become very slow in solving the equation so practically unusable.
% A better approach seems to be using an estimate of M_MN (based on the
% number of observed images) and estimate p=c/M_MN.
% 
% Shujun Li @ www.hooklee.com 2017

% Create a number of different styles for plotting ROC curves.
colors = {'b','r','g','m','c','y','k'};
markers = {'*','+','o','s'};
i = 0;
style_number = numel(colors) * numel(markers);
styles = cell(1, style_number);
for c=1:numel(colors)
    for s=1:numel(markers)
        i = i + 1;
        styles{i} = sprintf('-%s%s', colors{c}, markers{s});
    end
end

c = 22;
fun = @(p,N) c./p.*(1-(1-p).^N);
Ns = [100:100:1000 2000:1000:5000];
p = 0.001:0.0001:0.01;
y = cell(1,numel(Ns));
figure;
legend_info = cell(1,numel(Ns));
for i=1:numel(Ns)
    y{i} = fun(p, Ns(i));
    loglog(x, y{i}, styles{mod(i-1,style_number)+1});
    % plot(p, y{i}, styles{mod(i-1,style_number)+1});
    legend_info{i} = ['$N=' num2str(Ns(i)) '$'];
    hold on;
end
xlabel('$p$', 'interpreter', 'latex');
ylabel('$\overline{N_{\mathrm{normal}}}$', 'interpreter', 'latex');
h = legend(legend_info, 'location', 'northeast');
set(h,'Interpreter','latex');
