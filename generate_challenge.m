function [C, t, valid_labels] = generate_challenge(TI, M, MN, c, n_max, t_max, TI_min)
%[C, t, valid_labels] = generate_challenge(TI, M, MN, c, n_max, t_max, TI_min)
% This function generates a challenge and labels of neutral images given
% all the parameters needed.
% 
% Input arguments:
% See main.m for meanings and default values.
% 
% Output arguments:
% C = challenge generated;
% t = number of trap images used
% valid_labels = Boolean labels of what images are counted (= NOT neutral)
% 
% Shujun Li @ www.hooklee.com 2017

if ~exist('M','var')
    M = 4033;
end
if ~exist('MN','var')
    MN = 8355;
end
M_MN = M + MN;
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

C = zeros(1,c);
% If there are some trap images, we need to add some first.
TI_size = numel(TI);
if TI_size>=TI_min
    t = randi([1 min(t_max,TI_size)],1);
    trap_images = randperm(TI_size, t);
    C(1:t) = TI(trap_images);
else
    t = 0; % no trap image
end
% Select at least one image from M and at least one from MN
C(t+1) = randi(M,1);
C(t+2) = M + randi(MN,1);
% Remove the two selected images from all candidates
candidates = 1:M_MN;
candidates(C(t+(1:2))) = [];
% Randomly select c-(t+2) other images
C(t+3:c) = candidates(randperm(M_MN-2, c-(t+2)));
% In real world systems the images in C should also be shuffled so that
% the order does not reveal information about the images. Here, since
% we do not use the order to derive any information, we do not shuffle
% all images but keep the structure (t trap images + 1 M image + 1 MN
% image + c-(t+2) any images) to make coding simpler.
% correct labels will be C<=M (1 = in M, 0 = in MN)
% Select n random neutral images (note that the trap images should never be
% selected as neutral.
n = randi([0 n_max],1);
neutral_images = t + randperm(c-t,n);
valid_labels = true(1,c);
valid_labels(neutral_images) = false;

end
