%% This example script will take you through the basics usage of the Fig class
clear; clc;

%% The Fig class
%
% The objective of this class of the GepocToolbox is to override Matlab's
% default figure() so that figures are created with (opinionated) sensible
% default values that provide a nice-looking figure.
%
% The objective is for the basic usage of the class to behave as a user
% familiar with Matlab's default plotting functions would expect.
% In fact, the idea is for the most basic usage to be as simple as just
% substituting the call to figure() for a call to Fig().
% Lets take a look at an example of this "most basic" usage.

%% Basic example
% Let us create a figure using the typical Matlab syntax.

x = 0:pi/100:2*pi;
y = sin(x);

figure(1);
plot(x, y);
title("Figure of a sine wave");
xlabel("x");
ylabel("sin(x)");

% Take a look at the appearance of the figure. It is, in our opinion, very basic.
% Lets now use the Fig class to improve its appearance.
% Simply substitute line 8, which read
%   figure(1);
% for
%   Fig(1);
% and take a look at the new result.

% The objective of this short example is to show you that Fig is designed so you
% can use it in your preexisting code, providing out-of-the-box visual improvements
% with minimal effort.
% However, there is much more you can do with Fig to help you create article-level
% figures. Take a look at example_Fig_intermediate.m to see more advanced features!

