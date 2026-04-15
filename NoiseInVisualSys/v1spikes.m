% This script loads a data file containing a summary of an electrophysiological 
% experiment conducted by then post doc Robbe Goris in the lab of Tony Movshon (NYU). 
% A drifting grating was placed in the receptive field of a V1 cell of an 
% anestethized macaque monkey. The stimulus orientation was varied and action 
% potentials were recorded extracellularly. Each orientation was repeated 
% multiple times over the course of the experiment. S is a structure with S.condVec 
% the stimulus orientation, S.trialDur the trial duration, S.spikeCounts the total 
% number of action potentials fired during the trial, and S.spikeTimes the exact 
% timing of each spike.


%% Clean slate
clear all
close all
clc

% Load data
load('exampleCell1')

S

disp('Read summary of data format – script, line 1–9!')

