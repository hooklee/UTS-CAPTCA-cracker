% This script estimates the key parameters p, p_t and M_MN of the
% UTS-CAPTCHA scheme.
% 
% Shujun Li @ www.hooklee.com 2017

% set paratemers and variables for the simulated UTS-CAPTCHA service.
header;

% Maximum number of challenges observed
N_C = 2000;

% Ground truth value of p.
p_gt = c/M_MN;
% Declare a real symbolic variable for solving p.
syms p;
% Ground truth value of p_t.
p_t_gt = zeros(1,N_C);

% The counters for all images (in real world case only observed images can
% be counted, so this array should actually be a dynamic one; here we use a
% static array to simplify the code).
counters = zeros(1,M_MN);
% This is a separate array of counters working with vpasolve when N goes
% beyond vpasolve_N_max.
counters_vpasolve = zeros(1,M_MN);
% Estimated parameters given i observed challenges
% Estimated p based on mean frequency of all observed images minus any
% candidate trap images.
p_mean = zeros(1,N_C);
% Estimated p based on median frequency of all observed images minus any
% candidate trap images.
p_median = zeros(1,N_C);
% Estimated p based on median frequency of all observed images.
p_median2 = zeros(1,N_C);
% Estimated p based on numerically solving the equation
% c(1-(1-p)^N)) = -mean(N_normal)
p_vpasolve = zeros(1,N_C);
% Set the step of conducting vpasolve and displaying results to show
% progress of the parameter estimation.
N_C_step = 10;
% Set the precision of vpasolve to balance speed and precision.
digits(6);
% vpasolve function becomes very slow when N becomes large. So we use a
% number of N-challenge sessions to better estimate mean(N_normal) so we
% can still use vpasolve with a smaller N.
vpasolve_N_max = 200;
% Estimated p_t based on maximum frequency of all observed images.
p_t_max = zeros(1,N_C);
% Estimated p_t based on k-means (using the cluster with higher frequency).
p_t_kmeans = zeros(1,N_C);
% The number of unique images observed so far
N_unique = zeros(1,N_C);
% This is a separate array of N_unique working with vpasolve when N goes
% beyond vpasolve_N_max.
N_unique_vpasolve = zeros(1,N_C);

% The initial accuracy of a bot for recognising a single image.
p_bot0 = 0.8;
% Lables of all M_MN images known to the bot with p_bots wrong lables.
labels_bot = labels_truth;
error_indices = randperm(M_MN, round((1-p_bot0)*M_MN));
labels_bot(error_indices) = ~labels_bot(error_indices);
% Record the number of passed challenges.
challenges_passed = 0;
% Store a separate counter and the index of the first challenge each
% candidate trap image appeared.
trap_counters = zeros(1,M_MN);
% Set the first index to be Inf for images which have not been
% identified as an candidate.
trap_1st_indices = Inf(1,M_MN);

% If an estimate of p_t is needed, the attacker will need to make an
% non-empty TI. This can be done by randomly responding to all challenges.
for i=1:N_C
    % Generate a random challenge with the trap image numer and neutral labels.
    [C, t, valid_labels] = generate_challenge(TI, M, MN, c, n_max, t_max, 1);
    R_truth = (C<=M); % The groud truth labels are for those small indices.
    % Increase the counter of all observed images by 1.
    counters(C) = counters(C) + 1;
    counters_vpasolve(C) = counters_vpasolve(C) + 1;
    % Increase the counter of candidate trap images by 1.
    C_trap = C(~isinf(trap_1st_indices(C)));
    trap_counters(C_trap) = trap_counters(C_trap) + 1;

    % If the estimated value of p_t is much (>10 times) larger than that of
    % p, stop responding to challenges to avoid TI keeping enlarging.
    if p_t_max(i)<=10*p_mean(i)
        % The bot is trying to respond to the challenge following its labels.
        R_bot = labels_bot(C);
        % The UTS-CAPTCHA service checks if the responses are all correct.
        if isequal(R_bot(valid_labels), R_truth(valid_labels))
            fprintf('Challenge %d: passed!\n', i);
            % Before it is set the value is Inf. After it is set it is
            % something smaller than the current i. So min() will work.
            trap_1st_indices(C) = min(trap_1st_indices(C), i);
            challenges_passed = challenges_passed + 1;
            % The UTS-CAPTCHA service updates TI.
            % Any mistmatches must be for neutral images.
            trap_images = C(R_bot~=R_truth);
            TI = cat(2, TI, trap_images);
            if ~isempty(TI)
                % Update the ground truth value of p_t.
                p_t_gt(i:end) = (1+min(numel(TI),t_max))/2/numel(TI);
            end
        end
    end
    
    % N_unique is a lower bound of M_MN and an approximate of N_normal.
    N_unique(i) = sum(counters>0);
    N_unique_vpasolve(i) = sum(counters_vpasolve>0);
    % Consider images which are have appeared at least once and not a
    % candidate trap image.
    counters_i = counters(counters>0 & trap_counters==0);
    p_mean(i) = mean(counters_i) / i;
    % An open question:
    % When median frequency is used, there is a clear and stable
    % pseudo-periodic pattern appearing in the results. The source of this
    % pattern remains unclear, but may be proved mathematically.
    p_median(i) = median(counters_i) / i;
    p_median2(i) = median(counters(counters>0)) / i;
    % Estimate p using vpasolve() function.
    if mod(i,vpasolve_N_max)==0
        % N_unique_mean = mean(N_unique_vpasolve(vpasolve_N_max:vpasolve_N_max:i));
        N_unique_mean = median(N_unique_vpasolve(vpasolve_N_max:vpasolve_N_max:i));
        p_temp = vpasolve(c*(1-(1-p)^vpasolve_N_max)/p == N_unique_mean, p, [0 1]);
    elseif (i<=vpasolve_N_max && mod(i,N_C_step)==0)
        p_temp = vpasolve(c*(1-(1-p)^i)/p == N_unique_vpasolve(i), p, [0 1]);
    else
        p_temp = p_vpasolve(i); % Take the last (current) value.
    end
    % Use the current value for all future ones as vpasolve will not run
    % for every i.
    p_vpasolve(i:end) = min(p_temp);
    % Estimate p_t based on the separate counters for candidate trap images.
    trap_counters_i = trap_counters(trap_counters>0);
    % Only estimate when there are some candidate trap images.
    if ~isempty(trap_counters_i)
        i_trap = i - trap_1st_indices(trap_counters>0);
        p_t_max(i) = max(trap_counters_i./i_trap);
        % Use k-means to estimate the value of p_t.
        if numel(i_trap)>2
            [~,ps] = kmeans(trap_counters_i'./i_trap', 2, 'Start', [p_mean(i);p_t_max(i)], 'EmptyAction', 'singleton');
            p_t_kmeans(i) = max(ps);
        else
            % Just one or two candidate trap images with at least one
            % occurrence, so no need to run k-means.
            p_t_kmeans(i) = max(trap_counters_i./i_trap);
        end
    end
    
    % Reset the counters for vpasolve every vpasolve_N_max challenges after
    % the they are used.
    if mod(i,vpasolve_N_max)==0
        counters_vpasolve = zeros(1,M_MN);
    end
    
    % Display some results to show the script is still running.
    if mod(i,N_C_step)==0
        fprintf('Challenge %d: M_MN = %d (N_unique = %d); p = %g (%g, %g); p_t = %g (%g, %g)\n', ...
            i, M_MN, N_unique(i), p_gt, p_mean(i), p_vpasolve(i), p_t_gt(i), p_t_max(i), p_t_kmeans(i));
    end
end

close all;

figure;
semilogy(1:N_C, p_mean, 'b-', 'LineWidth', 3);
hold on;
semilogy(1:N_C, p_median, 'c-', 'LineWidth', 3);
semilogy(1:N_C, p_median2, 'y-', 'LineWidth', 3);
semilogy(1:N_C, p_vpasolve, 'm-', 'LineWidth', 3);
semilogy(1:N_C, c./N_unique, 'g-', 'LineWidth', 2);
line([1 N_C], [p_gt p_gt], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
legend({'Mean frequency', 'Median frequency (non-trap images only)', ...
    'Median frequency (all images)', 'Estimated using vpasolve', ...
    '$c/N_{\mathrm{unique}}$', 'Ground truth'}, 'Interpreter', 'latex', 'FontSize', 12);
axis tight;
set(gca,'XLim',[0 N_C]);
grid on;
xlabel('Number of challenges observed', 'Interpreter', 'latex');
ylabel('Estimated values of $p$', 'Interpreter', 'latex');

figure;
semilogy(1:N_C, round(c./p_mean), 'b-', 'LineWidth', 3);
hold on;
semilogy(1:N_C, round(c./p_median), 'g-', 'LineWidth', 3);
semilogy(1:N_C, round(c./p_median2), 'c-', 'LineWidth', 3);
semilogy(1:N_C, round(c./p_vpasolve), 'm-', 'LineWidth', 3);
semilogy(1:N_C, N_unique, 'y-', 'LineWidth', 2);
line([1 N_C], [M_MN M_MN], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 2);
legend({'$N_{\mathrm{unique}}$', '$c/p_{\mathrm{mean}}$', '$c/p_{\mathrm{median}}$', '$c/p_{\mathrm{median2}}$', ...
    '$c/p_{\mathrm{vpasolve}}$', 'Ground truth'}, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Number of challenges observed', 'Interpreter', 'latex');
ylabel('Estimated values of $|\mathrm{M}\cup\mathrm{MN}|$', 'Interpreter', 'latex');
axis tight;
grid on;

figure;
plot(1:N_C, p_t_max, 'b-', 'LineWidth', 3);
hold on;
plot(1:N_C, p_t_kmeans, 'c-', 'LineWidth', 3);
plot(1:N_C, p_t_gt, 'r--', 'LineWidth', 2);
legend({'Maximum frequency', 'Estimated via $k$-means', ...
    'Ground truth'}, 'Interpreter', 'latex', 'FontSize', 12);
axis tight;
set(gca,'XLim',[0 N_C]);
grid on;
xlabel('Number of challenges observed', 'Interpreter', 'latex');
ylabel('Estimated values of $p_t$', 'Interpreter', 'latex');
