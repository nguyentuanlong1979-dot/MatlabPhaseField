function [A] = Boundary_Neumann(A,Nx,Ny)
    A(1,:)=A(2,:);
    A(Nx,:)=A(Nx-1,:);
    A(:,1)=A(:,2);
    A(:,Ny)=A(:,Ny-1);
end