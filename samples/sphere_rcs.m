clear; close all; clc;

path_samples = fileparts(which('rcs_unit_sphere.m')); % /src/scripts/下脚本路径
path_root = fileparts(path_samples);

addpath(genpath(path_root)); % 将根目录添加到运行路径
cd(path_root);

% 
mesh = init_mesh('sph1');


% freqs = 2*pi*linspace(0,2e1,21);
freq = 3e8;
omega = 2 * pi * freq;
fprintf( 'Solving for frequency %.2e Hz\n', freq );

% Solver options
opts = init_solvopts( omega );

k = omega * sqrt(mu0 * eps0);

% Number of the edges
nedges = size(mesh.edges,1);

fintg_fp  = @(r, robs)integ_fp(k, r, robs, 8);
fintg_p   = @(r, robs)integ_p(k, r, robs, 8);

Z1 = mkmommat(mesh, fintg_fp, 1, 1:nedges, 1:nedges)...
    *1i*omega*mu0/(4*pi);
Z2 = mkmommatgrad(mesh, fintg_p, 1, 1:nedges, 1:nedges)...
    /(4*pi*1i*omega*eps0);
Z = Z1 + Z2;

% a = mkmommatgrad(mesh, fintg_p, 1, 1:nedges, 1:nedges)...
%     /(4*pi*1i*omega*eps0);
