% This script calculates the error in SAR calculation due to the finite
% directivity of directional couplers. Needs an VOP-file and an S-Parameter
% matrix as inputs. (See provided examples).
% Written by Stephan Orzada @ German Cancer Research Center (DKFZ)

[file_full,path_full]=uigetfile('..\.mat','Select S-Parameters','MultiSelect','off');
load([path_full file_full],'S')

[file_full,path_full]=uigetfile('..\.mat','Select VOP-file','MultiSelect','off');
load([path_full file_full],'VOP')

%VOP=diag(ones(8,1)); %Can be used to check influence of S-Matrix only (Power Error)

N_startvec=3; %Number of starting vectors for first optimization round. Although optimization converges well, it does not hurt to start several times.
N_test=360; %Number of samples of values between 0° and 360°
phi_vec=linspace(0,2*pi,N_test); %This vector contains all phases that will be tested.


for directivity_in_dB=[20,23,25,27,30,35,40] %Several directivities for comparison.

    directivity=10^(directivity_in_dB/20); 
    
    [Nch,~,Nvops]=size(VOP);
    
    S=(S+S.')/2; %Correct S for asymmetries. (S is a symmetric complex matrix) This can be necessary, as not all simulation programs provide symmetric matrices.
    result=zeros(N_test,1);

    X_max=zeros(Nch,N_test);
    tic
    parfor a=1:N_test % Calculate Points on phi-vector

        Merror=exp(1i*phi_vec(a))*(1/directivity)*S+eye(Nch); %This is the Error Matrix used to include the measurement error in the VOPs

        Qmeas=pagemtimes(Merror, 'ctranspose', pagemtimes(VOP, Merror), 'none'); %These VOPs now include measurement error.
        for b=1:Nvops
            % Here we correct for potential asymmetries again, as due to the
            % different order of operations in the previous step, the upper
            % and the lower triangle might be slightly different.
            % (A VOP is a hermitian matrix!
            % Therefore we need the complex conjugate.)
            Qmeas(:,:,b)=(Qmeas(:,:,b)+Qmeas(:,:,b)')/2; 
        end

        [R, X_full] = rQstar(VOP, Qmeas,0,[],N_startvec);

        [result(a),Indx]=max(R);
        X_max(:,a)=X_full(:,Indx);
        disp(['Step ' num2str(a) ' of ' num2str(N_test) '. R=' num2str(max(R))])% '. Rmax=' num2str(max(result))])
    end
    figure(1)
    plot(result)

    elapsed_time=toc;
    
    save(['DiCo_error_' num2str(directivity_in_dB) 'db_' file_full],"result","X_max","elapsed_time")
end