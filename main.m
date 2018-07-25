close all; clear all;
run('../vlfeat-0.9.20/toolbox/vl_setup');
addpath('transforms');

dataPath= "images/";
resultsPath = "images/";

%%{
data_path_ = [data_path 'APAP-railtracks/'];
exp_path = [results_path 'APAP-railtracks/'];

% find all image files in the provided data folder
data_files = [dir([data_path_ '*.jpg']);dir([data_path_ '*.JPG']);dir([data_path_ '*.png']);dir([data_path_ '*.PNG'])];

im1 = imread([data_path_ data_files(1).name]);
im2 = imread([data_path_ data_files(2).name]);

peak_thresh = 0;
edge_thresh = 500;
ransac_threshold = 0.015;
run_LFA;
%}

%%{
data_path_ = [data_path 'REW_racetracks/'];
exp_path = [results_path 'REW_racetracks/'];

% find all image files in the provided data folder
data_files = [dir([data_path_ '*.jpg']);dir([data_path_ '*.JPG']);dir([data_path_ '*.png']);dir([data_path_ '*.PNG'])];

im1 = imread([data_path_ data_files(1).name]);
im2 = imread([data_path_ data_files(2).name]);

peak_thresh = 0;
edge_thresh = 100;
ransac_threshold = 0.08;
run_LFA;
%}

%%{
data_path_ = [data_path 'REW_tower/'];
exp_path = [results_path 'REW_tower/'];

% find all image files in the provided data folder
data_files = [dir([data_path_ '*.jpg']);dir([data_path_ '*.JPG']);dir([data_path_ '*.png']);dir([data_path_ '*.PNG'])];

im1 = imread([data_path_ data_files(1).name]);
im2 = imread([data_path_ data_files(2).name]);

peak_thresh = 0;
edge_thresh = 500;
ransac_threshold = 0.04;
run_LFA;
%}

%%{
data_path_ = [data_path 'TFT_grove/'];
exp_path = [results_path 'TFT_grove/'];

% find all image files in the provided data folder
data_files = [dir([data_path_ '*.jpg']);dir([data_path_ '*.JPG']);dir([data_path_ '*.png']);dir([data_path_ '*.PNG'])];

im1 = imread([data_path_ data_files(1).name]);
im2 = imread([data_path_ data_files(2).name]);

peak_thresh = 0.007;
edge_thresh = 500;
ransac_threshold = 0.01;
run_LFA_fisheye;
%}

%%{
data_path_ = [data_path 'TFT_stairs/'];
exp_path = [results_path 'TFT_stairs/'];

% find all image files in the provided data folder
data_files = [dir([data_path_ '*.jpg']);dir([data_path_ '*.JPG']);dir([data_path_ '*.png']);dir([data_path_ '*.PNG'])];

im1 = imread([data_path_ data_files(1).name]);
im2 = imread([data_path_ data_files(2).name]);

peak_thresh = 0.007;
edge_thresh = 500;
ransac_threshold = 0.06;
run_LFA_fisheye;
%}