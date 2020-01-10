% This script sets the common parameters and variables for all scripts
% in this folder.
%
% Shujun Li @ www.hooklee.com 2017

% Use parameters in Kwon and Cha's paper.
% Images in the set M \cup MN will be indexed starting from 1 (1st to M-th
% images are in the set M, and (M+1)-th to (M+MN)-th images in the set MN).
M = 4033;
MN = 8355;
M_MN = M + MN;
labels_truth = [true(1,M) false(1,MN)];
c = 22;
n_max = 8;
t_max = 2;
% The minimum size of TI to start adding trap images.
TI_min = 1;

% Set the trap image database TI to be empty at the beginning.
TI = [];
% true = remove a trap image if it is successfully answered (consider a
% different user comes)
% flase = do not remove any trap images
% remove_trap_images = false;
remove_trap_images = true;
