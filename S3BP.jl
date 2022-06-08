using DifferentialEquations
using Dates

using Plots;

#use GR module
gr();

μ = 3.036e-6;
x_L4 = 0.5 - μ;
y_L4 = sqrt(3)/2;

u0 = [x_L4,y_L4,0.0,0.0];