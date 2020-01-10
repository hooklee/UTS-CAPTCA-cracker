function sr = get_success_rate(labels_bot, N_sr, TI, M, MN, c, n_max, t_max, TI_min)
%sr = get_success_rate(labels_bot, N_sr, TI, M, MN, c, n_max, t_max, TI_min)
% Get the success rate of a bot equipped with labels_bots as its current
% classification labels.
%
% Input arguments:
% N_sr: the number of challenegs used to estimate the success rate
% Other parameters: the same as generate_challenge() function
%
% Output argument: the success rate
% 
% Shujun Li @ www.hooklee.com 2017

sr = 0;

if nargin<1
    disp('Please provide at least one input argument!');
    return;
end

if ~exist('N_sr','var')
    N_sr = 2000;
end
if ~exist('M','var')
    M = 4033;
end
if ~exist('MN','var')
    MN = 8355;
end
if ~exist('c','var')
    c = 22;
end
if ~exist('n_max','var')
    n_max = 8;
end
if ~exist('t_max','var')
    t_max = 2;
end
if ~exist('TI_min','var')
    TI_min = 1;
end

challenges_passed = 0;
for j=1:N_sr
    [C, ~, valid_labels] = generate_challenge(TI, M, MN, c, n_max, t_max, TI_min);
    R_truth = (C<=M); % The groud truth labels are for those small indices.
    R_bot = labels_bot(C);
    if isequal(R_bot(valid_labels), R_truth(valid_labels))
        challenges_passed = challenges_passed + 1;
    end
end

sr = challenges_passed / N_sr;

end
