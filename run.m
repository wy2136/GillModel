% This matlabe script is numerically solving the linear beta-plane
% equatorial shallow water problem for the third problem of GFD final of fall 2008
% in the Department of Earth and Environmental Sciences at Columbia University
%
% Author: Wenchang Yang
% Email: yang.wenchang@columbia.edu

clear,clc

% the basic parameters
T = 50;
Ly = 5*pi;
Lx = 12*pi;
Lz = pi;
e = 0.1;
Nt = 500;
Ny = 60;
Nx = 128;
Nz = 20;
% choose if double the resolution
header = {['Parameters: ' ...
    '   T = ' num2str(T) ...
    ';   Lx = ' num2str(Lx) ...
    ';   Ly = ' num2str(Ly) ]; ...   
    ['Current resolution: ' ...
    '   Nt = ' num2str(Nt) ...
    ';   Nx = ' num2str(Nx) ...
    ';   Ny = ' num2str(Ny)] };
s{1} = 'Keep the resolution unchanged';
s{2} = 'Double time resolution';
s{3} = 'Double zonal resolution';
s{4} = 'Double meridional resolution';
icase = menu(header,s);
switch icase
    case 2
        Nt = 2*Nt;
    case 3
        Nx = 2*Nx;
    case 4
        Ny = 2*Ny;
end;
dt = T/Nt;
dy = Ly/(Ny+1);
y = dy*( -Ny/2:Ny/2 )';
dx = Lx/Nx;
x = dx*(-Nx/2:Nx/2-1)';
dz = Lz/(Nz-1);
z = dz*(0:Nz-1)';
dv = 0.2; mv = 10; cline = [-mv:dv:-dv dv:dv:mv];%choose which contour lines to be shown

% the initial values
v0 = zeros(Nx,Ny);
u0 = zeros(Nx,Ny+1);
eta0 = zeros(Nx,Ny+1);
Q = zeros(Nx,Ny+1);
% choose a case to run
header = 'Choose a case';
ss{1} = 'Gill: symmetric heating';
ss{2} = 'Gill: asymmetric heating';
ss{3} = 'Kelvin wave';
ss{4} = 'Gaussian shape initial perturbation';
ss{5} = 'Matsuno: forced stationary motion';
ss{6} = 'Gill: symmetric heating with Gaussian shape in zonal direction';
ss{7} = 'Gill: symmetric heating + asymmetric heating';
ss{8} = 'Gill: delta heating at y = 1';
icase = menu(header,ss);
switch icase
    case 1
        x = x + 2*pi;
        L = 2;
        Fx = zeros(Nx,1);
        Fx(abs(x)<L) = cos(pi/2/L*x(abs(x)<L));
        Q = Fx*( exp(-y'.*y'/4) );
        fname = 'Gill_sym';
    case 2
        L = 2;
        Fx = zeros(Nx,1);
        Fx(abs(x)<L) = cos(pi/2/L*x(abs(x)<L));
        Q = Fx*( y'.*exp(-y'.*y'/4) );
        fname = 'Gill_asym';
    case 3
        T = 10;
        Nt = 100;
        dt = T/Nt;
        e = 0;
        eta0 = cos(4*2*pi/Lx*x)*( exp(-y'.*y'/4) );
        u0 = eta0;
        fname = 'KelvinWave';
    case 4
        e = 0;
        eta0 = exp(-x.*x/4)*(  exp(-y'.*y'/4) );
        dv = 0.05; mv = 1; cline = -mv:dv:mv;%choose which contour lines to be shown
        fname = 'GaussPert';
    case 5
        e = 0.2;
        Q = sin(0.5*x)*( exp(-y'.*y'/4) );
        fname = 'MatsunoHeating';
    case 6
        x = x + 2*pi;
        L = 2;
        Q = exp(-x.*x/2)*( exp(-y'.*y'/4) );
        fname = 'Gill_sym_2';
    case 7
        x = x + 2*pi;
        L = 2;
        Fx = zeros(Nx,1);
        Fx(abs(x)<L) = cos(pi/2/L*x(abs(x)<L));
        Q = Fx*( exp(-y'.*y'/4) + y'.*exp(-y'.*y'/4) );
        fname = 'Gill_sym_asym';    
    case 8
        x = x + 2*pi;
        L = 2;
        Fx = zeros(Nx,1);
        Fx(abs(x)<L) = cos(pi/2/L*x(abs(x)<L));
        Fy = zeros(size(y)); Fy(abs(y-1)<0.5) = 1;
        Q = ones(Nx,1)*Fy' + Fx*exp(-y'.*y'/4);
        fname = 'Gill_sym_delta';        
end;
% Fourier transform
vk0 = fft(v0);
uk0 = fft(u0);
etak0 = fft(eta0);
Qk = fft(Q);
vk = vk0;
uk = uk0;
etak = etak0;

% the transform matrices
par.dt = dt;
par.e = e;
par.Ly = Ly;
par.Ny = Ny;
par.Lx = Lx;
par.Nx = Nx;
A = matgen(par);
Aeta_v = A.Aeta_v;%%%v
Au_v = A.Au_v;
Av_v = A.Av_v;
Aqv = A.Aqv;
Avv_eta = A.Avv_eta;%%%%eta
Aeta_eta = A.Aeta_eta;
Au_eta = A.Au_eta;
Aqeta = A.Aqeta;
Avv_u = A.Avv_u;%%%%u
Aeta_u = A.Aeta_u;
Au_u = A.Au_u;
Aqu = A.Aqu;

% create an avi file named 'myavi' if iavi==1
iavi = 1;
dN = 2;
if iavi==1;
    %aviobj = avifile(fname);
    mov = VideoWriter(fname);
    open(mov)
end;

% the integration loop
for ii = 1:Nt;
    vk_ = vk;
    for jj=1:Nx;
        vk(jj,:) = etak(jj,:)*Aeta_v{jj} + uk(jj,:)*Au_v{jj} + vk_(jj,:)*Av_v{jj} + Qk(jj,:)*Aqv{jj};
        uk(jj,:) = ( vk(jj,:) + vk_(jj,:) )/2*Avv_u{jj} + etak(jj,:)*Aeta_u{jj} +uk(jj,:)*Au_u{jj} + Qk(jj,:)*Aqu{jj};
        etak(jj,:) = ( vk(jj,:) + vk_(jj,:) )/2*Avv_eta{jj} + etak(jj,:)*Aeta_eta{jj} + uk(jj,:)*Au_eta{jj} + Qk(jj,:)*Aqeta{jj};
    end;
    eta = real( ifft(etak) );
    u = real(ifft(uk));
    v = real(ifft(vk)); v = ([v zeros(Nx,1)]+[zeros(Nx,1) v])/2;
    figure(1);
    contour(x,y,eta',cline);
    axis equal;
    axis([x(1) x(end) y(1) y(end)])
    hold on;
    Nq = 5;
    quiver(x(1:Nq:end),y(1:Nq:end),u(1:Nq:end,1:Nq:end)',v(1:Nq:end,1:Nq:end)',0);
    hold off;
    title({ss{icase}; [ 't = ' num2str(ii*dt,'%4.2f')]});
    xlabel('x'); ylabel('y');
    if iavi==1
        if mod(ii,dN)==0
            %aviobj = addframe( aviobj,getframe(gcf) );
            writeVideo(mov, getframe(gcf));
        end;
        if ii==Nt
            %aviobj = close(aviobj);
            close(mov);
        end;
    end;
    drawnow
end;
print(gcf,'-dpng', '-r300', [fname, '_i', num2str(ii)]);

% plot the stream function or zonal mean flow in x-z and y-z cross section
% in case 1 and case 2 (Gill cases).
if icase==1||icase==2;
    figure(2);
    subplot(3,1,1);% stream function in x-z cross section
    uym = mean(u,2)*cos(pi/Lz*z');
    psiuw = -Lz/pi*mean(u,2)*sin(pi/Lz*z');
    contour(x,z,psiuw',-0.4:0.04:0.4);
    title({ss{icase}; [ 't = ' num2str(ii*dt,'%4.2f')]; 'u-w'});
    xlabel('x'); ylabel('z');
    set(gca,'ytick',[0 z(end)],'yticklabel',{'0','z'});

    subplot(3,1,2);%stream function in y-z cross section
    psivw = -Lz/pi*mean(v)'*sin(pi/Lz*z');
    contour(y,z,psivw',-0.1:0.004:0.1);
    title('v-w');
    xlabel('y'); ylabel('z');
    set(gca,'ytick',[0 z(end)],'yticklabel',{'0','z'});

    subplot(3,1,3);%zonal mean flow
    uxm = mean(u)'*cos(pi/Lz*z');
    contour(y,z,uxm',-0.4:0.04:0.4);
    title('u');
    xlabel('y'); ylabel('z');
    set(gca,'ytick',[0 z(end)],'yticklabel',{'0','z'});
    print(gcf,'-dpng', '-r300', [fname, '_i', num2str(ii), '_CrossSection']);
end;
