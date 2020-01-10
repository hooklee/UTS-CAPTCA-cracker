% This script simulates the sequential chi^2 tests based attack on the
% UTS-CAPTCHA scheme.
% 
% Shujun Li @ www.hooklee.com 2017

% set paratemers and variables for the simulated UTS-CAPTCHA service.
header;

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
% Record trap images already detected
trap_images_detected = [];
% Record all images seen in all passed challenges (candidate trap images)
images_C_passed = [];

% Number of total challenges for learning
N_C = 1000*1000;
% The accuracy of the bot w.r.t. the number of challenges observed.
p_bot = p_bot0*ones(1,N_C);
% The threshold used to judge an image is a trap image with the given
% significance level for a chi^2 test.
alpha = 0.01;
% The number of challenges observed for a chi^2 test to run. This needs
% large enough to let trap images to stand out of normal ones. It can
% be estimated based on an lower bound of p_t>=1/|TI| where |TI|
% can be estimated based on trap images detected (which will
% normally be smaller than |TI| due to false negatives).
counter_min = 10;
% The ratio to select large counter values (times of median count).
counter_median_ratio = 2;

% Enable or disable calculation of success rate empirically
calculate_sr = false;
if calculate_sr
    % Number of total challenges for checking the success rate
    N_sr = 2000; %#ok<UNRCH>
    % Record success rates
    success_rate0 = get_success_rate(labels_bot, N_sr, TI, M, MN, c, n_max, t_max, TI_min);
    success_rates = success_rate0*ones(1,N_C);
end

counters = zeros(1,c);
for i=1:N_C
    % Generate a random challenge with the trap image numer and neutral labels.
    [C, t, valid_labels] = generate_challenge(TI, M, MN, c, n_max, t_max, 1);
    R_truth = (C<=M); % The groud truth labels are for those small indices.

    if learning_mode
        % Learning mode, the bot will not try to respond but just passively
        % observing challenges until an existing condition is met.
        % Increase counters for images shown in the last successfully
        % passed challenge and also in the current challenge.
        C_last2_in_C = ismember(C_last2,C);
        counters(C_last2_in_C) = counters(C_last2_in_C) + 1; 
        n = i-i_last;
        % Run chi^2 test only when enough challenges are observed
        if (n>N_chi2)
            trap_images_number = 0;
            trap_images_detected2 = intersect(trap_images_detected, C_last2);
            % A true trap image should appear at least counter_min
            % times and if not it is likely to be a normal image
            % accidentally appearing that many times (which can
            % always happen but with a very small probability).
            while max(counters)>=counter_min
                E = sum(counters)/numel(counters);
                chi2 = sum((counters-E).^2/E);
                chi2_threshold = chi2inv(1-alpha,numel(counters)-1);
                if (chi2>chi2_threshold && trap_images_number<=c)
                    % Get the images with the highest occurrence
                    % counts.
                    trap_images_indices = (counters==max(counters));
                    trap_images = C_last2(trap_images_indices);
                    % Remove the identified trap images.
                    C_last2(trap_images_indices) = []; %#ok<SAGROW>
                    counters(trap_images_indices) = [];
                    % Remove known (previously detected trap images
                    % as they are double confirmed now.
                    trap_images_known = intersect(trap_images, trap_images_detected2);
                    trap_images_detected2 = setdiff(trap_images_detected2, trap_images_known);
                    trap_images = setdiff(trap_images, trap_images_known);
                    if ~isempty(trap_images)
                        trap_images_number = trap_images_number + numel(trap_images);
                        % Only need to reverse the original labels
                        % because for such images the old labels must
                        % be wrong.
                        labels_bot(trap_images) = ~labels_bot(trap_images);
                        % After reversing the labels the new ones
                        % should be the same as the groud truth.
                        fprintf('Challenge %d: trap image(s) found -', i);
                        for j=1:numel(trap_images)
                            fprintf(' %d', i, trap_images(j));
                            if labels_bot(trap_images(j))~=labels_truth(trap_images(j))
                                fprintf(' (WRONG!)');
                            end
                        end
                        fprintf('\n');
                        % Add the detected trap images into an internal
                        % list so that they can checked for false
                        % positives.
                        trap_images_detected = cat(2, trap_images_detected, trap_images);
                    end
                else
                    break; % No trap images left
                end
            end
            if ~isempty(trap_images_detected2)
                % Some previously detected trap images are not
                % double confirmed. These are likely false
                % positives so should be removed.
                % In addition, detected trap images will be eventually
                % removed by the UTS-CAPTCHA service thus should be removed
                % from the bot's list of detected trap images as well.
                trap_images_detected = setdiff(trap_images_detected, trap_images_detected2);
                fprintf('Challenge %d: trap images removed - %s\n', i, mat2str(trap_images_detected2));
            end
            % Update the current accuracy of the bot's image
            % recognition capability if any trap image is detected.
            if trap_images_number>0
                p_bot(i:end) = sum(labels_bot==labels_truth)/M_MN;
            end
            % Exit learning mode.
            if calculate_sr
                success_rate = get_success_rate(labels_bot, N_sr, TI, M, MN, c, n_max, t_max, TI_min); %#ok<UNRCH>
                success_rates(i:end) = success_rate;
            else
                success_rate = p_bot(i)^(22-n_max/2);
            end
            fprintf('Challenge %d: image classification accuracy = %g => CAPTCHA success rate = %g\n', i, p_bot(i), success_rate);
            learning_mode = false;
        end
    else
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
            % pass => The bot switches to learning mode
            learning_mode = true;
            % Need to remember the index of the last successful challenge
            % to exit the learning mode when exit condition is met.
            i_last = i;
            C_last = C;
            % The following are for debugging purposes.
            t_last = t;
            valid_labels_last = valid_labels;
            R_truth_last = R_truth;
            % Set an array of counters (number of occurrence of images) for
            % all images shown in the passed challenge.
            C_last2 = C;
            counters = zeros(1, numel(C_last2));
            % Update the number of challenges observed to detect trap
            % images assuming numel(trap_images_detected) may be slightly
            % smaller than |TI|.
            % N_chi2 = counter_min*(numel(trap_images_detected)+c/2);
            N_chi2 = counter_min*(numel(trap_images_detected)+c/4);
        else
            % Failed
            % When there are missed trap images (false negatives), the bot
            % will indeed be always locked out when at least one such
            % missed trap images appears in the challenge.
            % The following two facts can help detect such missed trap
            % images:
            % 1) once seen (and missed) they never appear in any future
            % passed challenges;
            % 2) once seen (and missed) they appears in future failed
            % challenges with the highest probability, more than any other
            % trap images and normal images.
            % Therefore, among all non-trap images which never appear in
            % any passed challenges, the most occurring image in failed
            % challenges is likely a missed trap image. Such images can be
            % tested by temporarily reversing their labels to respond to a
            % challenge in which they appear again.
            % Code to be added.... (A bit complicated)
            % This locking-out issue can also be solved by using multiple
            % parallel bots running from different machines, and their
            % results can be aggregated. In this case even if all are
            % eventually locked out, the merger of the multiple bots will
            % lead to a more accurate bot.
        end
    end
end

% Get the indices corresponding to changes of accuracy.
N_C_indices = find(diff([0 p_bot])>0);
if N_C_indices(end)<N_C
    N_C_indices = cat(2, N_C_indices, N_C);
end
p_bot_indices = p_bot(N_C_indices);
if calculate_sr
    success_rates_indices = success_rates(N_C_indices); %#ok<UNRCH>
end

plot(N_C_indices, p_bot_indices, 'b-*');
hold on;
if calculate_sr
    plot(N_C_indices, success_rates_indices, 'r-o'); %#ok<UNRCH>
else
    plot(N_C_indices, p_bot_indices.^(c-n_max/2), 'r-o');
end
xlabel('Number of challenges');
ylabel('Bot''s performance');
legend({'Single-image classification accuracy', 'Success rate of passing UTS-CAPTCHA'}, ...
    'Interpreter', 'latex', 'FontSize', 12, 'Location', 'SouthEast');
axis tight;
axis([1 N_C 0 1]);
grid on;
