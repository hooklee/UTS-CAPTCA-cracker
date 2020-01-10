% This script simulates the k-means based attack on the
% UTS-CAPTCHA scheme.
% 
% Shujun Li @ www.hooklee.com 2017

% set paratemers and variables for the simulated UTS-CAPTCHA service.
header;

% Assume p is known (a separate script will be used to show how p can be
% estimated from N_normal).
p = c/M_MN;
% The following code is for estimating p from N_normal.
% syms p; p = vpasolve(c*(1-(1-p)^N)/p == N_normal, p, [-Inf Inf]);
% p_t = (1+t_max)/2/numel(TI);

% The initial accuracy of a bot for recognising a single image.
% p_bot0 = 0.6;
p_bot0 = 0.8;
% p_bot0 = 0.825;
% p_bot0 = 0.9;
% p_bot0 = 0.95;
% Lables of all M_MN images known to the bot with p_bots wrong lables.
labels_bot = labels_truth;
error_indices = randperm(M_MN, round((1-p_bot0)*M_MN));
labels_bot(error_indices) = ~labels_bot(error_indices);

% Number of total challenges for learning
N_C = 200*1000;
% The accuracy of the bot w.r.t. the number of challenges observed.
p_bot = p_bot0*ones(1,N_C);
% The lower bound of p_t/p.
p_t_p_min = 10;
% The minimum number of observed challenges to run a k-means.
N_kmean_min = ceil(p_t_p_min);
% Boolean labels for candidate images
candidate_images = false(1,M_MN);
% Counters of candidate images.
counters = zeros(1,M_MN);
% Index of the last successfully passed challenge
last_pass_index = 0;
% The starting index of candidate images for counting
start_index = 0;
% Record trap images already detected
trap_images_detected = [];

for i=1:N_C
    % The UTS-CAPTCHA service generates a random challenge with the trap
    % image numer and neutral labels.
    [C, t, valid_labels] = generate_challenge(TI, M, MN, c, n_max, t_max, 1);
    R_truth = (C<=M); % The groud truth labels are for those small indices.

    % The bot increases counters by one for all images appearing
    % in the current challenge.
    counters(C) = counters(C)+1;

    i_start_index = i-start_index;
    if (i-start_index>=N_kmean_min)
        % The bot conducts k-means clustering if there are enough data and
        % there are some candidate images appearing at least once.
        if (any(candidate_images) && any(counters(candidate_images)>0))
            indices = find(counters>0);
            frequencies = counters(indices)'/i_start_index;
            [cluster_labels,p_pt] = kmeans(frequencies, 2, 'Start', [min(frequencies);max(frequencies)], 'EmptyAction', 'singleton');
            if (numel(p_pt)==2 && p_pt(2)/p_pt(1)>=p_t_p_min)
                % The bot detects all candidate images in the
                % high-frequency cluster as trap images.
                trap_images = indices(cluster_labels==2);
                trap_images = trap_images(candidate_images(trap_images));
                if ~isempty(trap_images)
                    % The bot reverses the classification labels of detected
                    % trap images.
                    labels_bot(trap_images) = ~labels_bot(trap_images);
                    trap_images_detected = cat(2, trap_images_detected, trap_images);
                    % The bot removes detected trap images from candidates.
                    candidate_images(trap_images) = false;
                    fprintf('Challenge %d: trap image(s) found -', i);
                    for j=1:numel(trap_images)
                        fprintf(' %d', i, trap_images(j));
                        if labels_bot(trap_images(j))~=labels_truth(trap_images(j))
                            fprintf(' (WRONG!)');
                        end
                    end
                    fprintf('\n');
                    p_bot(i:end) = sum(labels_bot==labels_truth)/M_MN;
                    success_rate = p_bot(i)^(22-n_max/2);
                    fprintf('Challenge %d: image classification accuracy = %g => CAPTCHA success rate = %g\n', i, p_bot(i), success_rate);
                end
            end
        end
        % The bot is trying to respond to the challenge following its
        % current labels.
        R_bot = labels_bot(C);
        % The UTS-CAPTCHA service checks if the responses are all correct.
        if isequal(R_bot(valid_labels), R_truth(valid_labels))
            fprintf('Challenge %d: passed!\n', i);
            % The UTS-CAPTCHA service updates TI.
            % Any mistmatches must be for neutral images.
            trap_images = C(R_bot~=R_truth);
            TI = cat(2, TI, trap_images);
            % The UTS-CAPTCHA service removes trap images if they are
            % answered successfully in this round (considered as user
            % change).
            if (remove_trap_images && t>0)
                TI = setdiff(TI, C(R_bot(1:t)==R_truth(1:t)));
            end
            % The bot sets the new start index and reset the counters for
            % all candidate images and add new candidate images.
            start_index = i;
            counters = zeros(1,M_MN);
            candidate_images(C) = true;
        end
    end
end

% Get the indices corresponding to changes of accuracy.
N_C_indices = find(diff([0 p_bot])>0);
if N_C_indices(end)<N_C
    N_C_indices = cat(2, N_C_indices, N_C);
end
p_bot_indices = p_bot(N_C_indices);

plot(N_C_indices, p_bot_indices, 'b-*');
hold on;
plot(N_C_indices, p_bot_indices.^(c-n_max/2), 'r-o');
xlabel('Number of challenges');
ylabel('Bot''s performance');
legend({'Single-image classification accuracy', 'Success rate of passing UTS-CAPTCHA'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'SouthEast');
axis tight;
axis([1 N_C 0 1]);
grid on;
