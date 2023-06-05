function [Nx, Nu] = fsf_matrices(system)

nx = size(system.A,1); % number of states               
ny= size(system.C,1); % number of outputs

big_A = [system.A-eye(nx) system.B
         system.C system.D];

big_Y =[ zeros(nx,ny)
         eye(ny,ny) ];

big_N = big_A\big_Y;

% full state feedback matrices
Nx = big_N(1:nx,:); 
Nu = big_N (nx+1:end,:); 

end