function A = matgen(par)
% This function generates the matrices that are used in the Fourier
% transformed linear beta-plane equatorial shallow water system for each k
% component: 
% v = eta_*Aeta_v + u_*Au_v + v_*Av_v - Q*Aqv; 
% eta = (v+v_)/2*Avv_eta + eta_*Aeta_eta + u_*Au_eta  - Q*Aqeta; 
% u = (v+v_)/2*Avu + eta_*Aeta_u + u_*Au_u- Q*Aqu;
% A is a structure variable with fields Aeta_v, Au_v and so on.

% the basic parameters
dt = par.dt;
e = par.e;
Ly = par.Ly;
Ny = par.Ny;
Lx = par.Lx;
Nx = par.Nx;
k = 2*pi/Lx*[0:Nx/2 -Nx/2+1:-1];
dy = Ly/(Ny+1);

% the basic matrices that build the transform matrices
E = (1/dt + e/2)*eye(Ny+1);
E_ = (1/dt - e/2)*eye(Ny+1);
A13 = ( [eye(Ny);zeros(1,Ny)] + [zeros(1,Ny);-eye(Ny)]  )/(2*dy);
A31 = -A13';
A23 = -( [diag(-Ny/2:Ny/2-1);zeros(1,Ny)] + [zeros(1,Ny);diag(-Ny/2+1:Ny/2)] )*dy/8;
A32 = -A23';

% initiate the transform matrices, which are functions of wavenumber k. In
% matlab, each one is expressed as a Nx long cell vector.
Aeta_v = cell(Nx,1);
Au_v = cell(Nx,1);
Av_v = cell(Nx,1);
Aqv = cell(Nx,1);

Avv_eta = cell(Nx,1);
Aeta_eta = cell(Nx,1);
Au_eta = cell(Nx,1);
Aqeta = cell(Nx,1);

Avv_u = cell(Nx,1);
Aeta_u = cell(Nx,1);
Au_u = cell(Nx,1);
Aqu = cell(Nx,1);

for ii=1:Nx;
    K = 1i*k(ii)/2*eye(Ny+1);
    Av_inv = eye(Ny)/( E(1:end-1,1:end-1) - (A31+A32)/(E+K)*(A13+A23)/2 - (A31-A32)/(E-K)*(A13-A23)/2 ).';
    
    Aeta_v{ii} = ( -A31 - (A31+A32)/(E+K)*(E_-K)/2 - (A31-A32)/(E-K)*(E_+K)/2 ).'*Av_inv;
    Au_v{ii} = ( -A32 - (A31+A32)/(E+K)*(E_-K)/2 + (A31-A32)/(E-K)*(E_+K)/2 ).'*Av_inv;
    Av_v{ii} = ( E_(1:end-1,1:end-1) + (A31+A32)/(E+K)*(A13+A23)/2 + (A31-A32)/(E-K)*(A13-A23)/2 ).'*Av_inv;
    Aqv{ii} = ( (A31+A32)/(E+K) + (A31-A32)/(E-K) ).'/2*Av_inv;
    
    Avv_eta{ii} = -( (E+K)\(A13+A23) + (E-K)\(A13-A23) ).';
    Aeta_eta{ii} = ( (E+K)\(E_-K) + (E-K)\(E_+K) ).'/2;
    Au_eta{ii} = ( (E+K)\(E_-K) - (E-K)\(E_+K) ).'/2;
    Aqeta{ii} = -( eye(Ny+1)/(E+K) + eye(Ny+1)/(E-K) ).'/2;
    
    Avv_u{ii} = -( (E+K)\(A13+A23) - (E-K)\(A13-A23) ).';
    Aeta_u{ii} = ( (E+K)\(E_-K) - (E-K)\(E_+K) ).'/2;
    Au_u{ii} = ( (E+K)\(E_-K) + (E-K)\(E_+K) ).'/2;
    Aqu{ii} = -( eye(Ny+1)/(E+K) - eye(Ny+1)/(E-K) ).'/2;
end;

% input these matrices into the structure variabe A
A.Aeta_v = Aeta_v;
A.Au_v = Au_v;
A.Av_v = Av_v;
A.Aqv = Aqv;

A.Avv_eta = Avv_eta;
A.Aeta_eta = Aeta_eta;
A.Au_eta = Au_eta;
A.Aqeta = Aqeta;

A.Avv_u = Avv_u;
A.Aeta_u = Aeta_u;
A.Au_u = Au_u;
A.Aqu = Aqu;       